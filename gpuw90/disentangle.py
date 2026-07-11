"""
disentangle.py -- Batched Souza-Marzari-Vanderbilt disentanglement in PyTorch.

Phase 2 of the GPU Wannier tool: for entangled bands (num_bands > num_wann) select,
at every k, the num_wann-dimensional subspace that minimises the gauge-invariant
spread Omega_I, then hand a dimensionally-reduced problem (M0_reduced, U0) to the
Phase-1 wannierise() -- whose code is unchanged.

The iteration mirrors wannier90 3.1 disentangle.F90 step for step (subroutines
dis_windows / dis_project / dis_proj_froz / dis_extract / internal_find_u):

  windows   : outer window selects ndimwin(k) bands (default = all), frozen window
              marks ndimfroz(k) bands that must stay in the subspace exactly.
  project   : initial subspace V0 = polar(A_windowed) via SVD  (SMV Sec. III.D);
              with frozen states, the non-frozen part comes from the leading
              eigenvectors of Q_f P_s Q_f  (SMV Sec. III.G, dis_proj_froz).
  extract   : Jacobi-style self-consistency loop (all k updated simultaneously,
              exactly like wannier90's MPI k-distributed loop):
                Z_in^(1)   = Z[V^(0)]
                Z_in^(i)   = mix * Z_out^(i-1) + (1-mix) * Z_in^(i-1)      (i > 1)
                V^(i)      = frozen states + top (nw - nf) eigenvectors of Z_in^(i)
                Z_out^(i)  = Z[V^(i)]
              with Z[V](k) = sum_b w_b  M^{k,b} P_{k+b} M^{k,b dagger},
              P = V V^dagger, projected onto the non-frozen complement at k.
              Omega_I(i-1) ("womegai1", SMV Eq. 18) uses previous-iteration
              neighbour states; Omega_I(i) ("womegai", SMV Eq. 12) uses updated
              states; delta = Omega_I(i-1)/Omega_I(i) - 1; converged when
              |delta| < dis_conv_tol for dis_conv_window consecutive iterations.
  find_u    : U0 = polar(V^dagger A_windowed)  (unitary factor of the SVD), i.e.
              wannier90's internal_find_u Loewdin projection.

Batching strategy (the whole point of the tool): the per-k ragged dimensions
(ndimwin, ndimfroz vary with k) are flattened into fixed [nb] tensors with masks;
out-of-window and frozen directions are made inert inside ONE batched Hermitian
eigendecomposition per k by diagonal level shifts (see _Z_SHIFT notes inline).
Everything heavy is torch.matmul / torch.linalg.eigh over the k (and neighbour)
batch dimensions -- CPU/CUDA same code via `device`.

Deliberate deviation from wannier90 (documented, gauge-level only): after
convergence wannier90 rotates V inside the fixed subspace to Hamiltonian
eigenvectors (dis_extract tail) -- a pure right-rotation V -> V R that cancels
between M0_reduced and U0 (M0' = R^dag M0_red R', U0' = R^dag U0), leaving every
physical quantity (Omega_I, centres, spreads, the wannierise trajectory in gauge-
invariant terms) bit-for-bit equivalent.  We skip it; band-structure interpolation
(the reason wannier90 does it) is out of scope here.
"""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
import torch


@dataclass
class DisentangleResult:
    M0_reduced: torch.Tensor   # [nk, nnb, nw, nw] complex  -> wannierise() M0
    U0: torch.Tensor           # [nk, nw, nw]      complex  -> wannierise() U0
    V: torch.Tensor            # [nk, nb, nw]      complex  subspace basis (window rows)
    omega_I: float             # final gauge-invariant spread (Ang^2), w90 'Final Omega_I'
    history: list              # (iter, Omega_I(i-1), Omega_I(i), delta_frac)
    n_iter: int
    converged: bool
    ndimwin: np.ndarray        # [nk] bands inside outer window
    ndimfroz: np.ndarray       # [nk] frozen bands


# --------------------------------------------------------------------------- #
# Window masks (dis_windows)
# --------------------------------------------------------------------------- #
def band_windows(eig, dis_win_min, dis_win_max, froz_win_min, froz_win_max,
                 num_wann):
    """Outer/frozen window membership masks, mirroring wannier90 dis_windows.

    eig must be sorted ascending per k (checked by io_w90.read_eig), so the
    >=min & <=max selection is automatically the contiguous [imin..imax] range
    that wannier90 constructs.

    Returns (win_mask, froz_mask) as bool ndarrays [nk, nb].
    """
    win = (eig >= dis_win_min) & (eig <= dis_win_max)
    ndimwin = win.sum(axis=1)
    if (ndimwin == 0).any():
        raise ValueError("dis_windows: outer energy window contains no eigenvalues "
                         f"at k = {np.nonzero(ndimwin == 0)[0].tolist()}")
    if (ndimwin < num_wann).any():
        k_bad = int(np.argmax(ndimwin < num_wann))
        raise ValueError(
            f"dis_windows: outer window has fewer states ({int(ndimwin[k_bad])}) "
            f"than num_wann ({num_wann}) at k index {k_bad}")

    if froz_win_max is None:
        froz = np.zeros_like(win)
    else:
        froz = win & (eig >= froz_win_min) & (eig <= froz_win_max)
        ndimfroz = froz.sum(axis=1)
        if (ndimfroz > num_wann).any():
            k_bad = int(np.argmax(ndimfroz > num_wann))
            raise ValueError(
                f"dis_windows: {int(ndimfroz[k_bad])} bands inside the frozen "
                f"window but only {num_wann} target WFs at k index {k_bad}")
    return win, froz


# --------------------------------------------------------------------------- #
# Initial subspace (dis_project + dis_proj_froz)
# --------------------------------------------------------------------------- #
def _polar_unitary(B):
    """Unitary polar factor of B via SVD: polar(B) = Z @ Vh  (B = Z S Vh).

    Identical to Loewdin B (B^dag B)^{-1/2} for full-rank B, but well-defined for
    (near-)rank-deficient B -- matching wannier90's ZGESVD-based construction.
    """
    Z, _, Vh = torch.linalg.svd(B, full_matrices=False)
    return Z @ Vh


def initial_subspace(A_w, win_t, froz_t):
    """Initial V0 [nk, nb, nw]: SVD-unitarised projections, frozen-aware.

    A_w is the window-masked projection matrix (zero rows outside the window),
    so the batched SVD reproduces wannier90's per-k ndimwin x nw SVD exactly
    (the padded rows just carry zero singular components).

    Frozen handling (dis_proj_froz): diagonalise Q_f P_s Q_f and take the leading
    nw - nf eigenvectors as the non-frozen part, the frozen band unit vectors as
    the rest.  Batched via the same level-shift trick used for the Z matrix:
    QPQ has spectrum in [0, 1]; adding +2 on frozen diagonals and -2 on
    out-of-window diagonals makes the top-nw eigenvectors of ONE batched eigh be
    exactly {frozen unit vectors} + {leading complement eigenvectors}.  This also
    reproduces wannier90's "ortho-fix" (eigenvectors with ~0 eigenvalue must be
    taken orthogonal to the frozen states): our complement block is exactly
    decoupled, so its eigenvectors have zero frozen components by construction.
    """
    V0 = _polar_unitary(A_w)                     # [nk, nb, nw]
    if not bool(froz_t.any()):
        return V0

    q = (win_t & ~froz_t)                        # [nk, nb] non-frozen window rows
    qf = q.to(A_w.real.dtype)
    fz = froz_t.to(A_w.real.dtype)
    ex = (~win_t).to(A_w.real.dtype)

    Ps = V0 @ V0.conj().transpose(-1, -2)        # [nk, nb, nb] projector onto span(A)
    QPQ = Ps * qf[:, :, None] * qf[:, None, :]
    shift = torch.diag_embed(2.0 * fz - 2.0 * ex).to(QPQ.dtype)
    evals, evecs = torch.linalg.eigh(QPQ + shift)

    nw = A_w.shape[-1]
    V0 = evecs[:, :, -nw:]                       # top nw: complement leaders + frozen
    V0 = _scrub_columns(V0, evals[:, -nw:], win_t, froz_t, thresh=1.5)
    return V0


def _scrub_columns(V, evals_top, win_t, froz_t, thresh):
    """Zero the fp-noise cross-block components of the selected eigenvectors.

    The shifted matrices are exactly block-diagonal (frozen rows/cols zeroed
    before shifting), so each selected eigenvector belongs either to the frozen
    block (eigenvalue near the +shift plateau, > thresh) or to the free block.
    LAPACK/eigh leaves O(eps*||A||/gap) ~ 1e-15 cross-components; scrubbing them
    keeps frozen states exactly frozen and V exactly window-supported.  The norm
    perturbation is O(1e-30) -- invisible in float64, so no renormalisation.
    """
    free_rows = (win_t & ~froz_t)                                # [nk, nb]
    is_froz_col = (evals_top > thresh)                           # [nk, nw]
    rowmask = torch.where(is_froz_col[:, None, :],
                          froz_t[:, :, None], free_rows[:, :, None])
    return V * rowmask.to(V.dtype)


# --------------------------------------------------------------------------- #
# dis_extract (batched Jacobi self-consistency on the Z matrix)
# --------------------------------------------------------------------------- #
def disentangle(
    M0_full,
    A_full,
    eig,
    bvec,
    wb,
    kb_idx,
    dis_win_min=None,
    dis_win_max=None,
    froz_win_min=None,
    froz_win_max=None,
    num_iter=200,
    mix_ratio=0.5,
    conv_tol=1e-10,
    conv_window=3,
    device="cpu",
    dtype=torch.complex128,
    check_every=1,
    verbose=False,
):
    """Extract the optimally-connected num_wann subspace (wannier90 dis_extract).

    Parameters
    ----------
    M0_full : [nk, nnb, nb, nb] complex   full-band overlaps (straight from .mmn)
    A_full  : [nk, nb, nw]      complex   full-band projections (straight from .amn)
    eig     : [nk, nb]          real      band energies in eV (from .eig)
    bvec, wb, kb_idx                      as in wannierise()
    dis_win_min/max   : outer window (eV); default = full eigenvalue range
                        (wannier90 parameters.F90 default)
    froz_win_min/max  : frozen (inner) window; froz_win_max=None -> no frozen
                        states; froz_win_min defaults to dis_win_min when only
                        froz_win_max is given (wannier90 default)
    num_iter, mix_ratio, conv_tol, conv_window : wannier90's dis_num_iter,
                        dis_mix_ratio, dis_conv_tol, dis_conv_window
    check_every       : poll the device-side convergence flag every N iterations
                        (single scalar device->host read per poll)

    Returns DisentangleResult; feed .M0_reduced / .U0 to wannierise().
    """
    rdtype = torch.float64 if dtype == torch.complex128 else torch.float32
    M0_full = M0_full.to(device=device, dtype=dtype)
    A_full = A_full.to(device=device, dtype=dtype)
    bvec = bvec.to(device=device, dtype=rdtype)
    wb = wb.to(device=device, dtype=rdtype)
    kb_idx = kb_idx.to(device=device, dtype=torch.long)
    eig = np.asarray(eig, dtype=np.float64)

    nk, nnb, nb, nb2 = M0_full.shape
    nw = A_full.shape[-1]
    if nb != nb2 or A_full.shape[1] != nb or eig.shape != (nk, nb):
        raise ValueError(f"shape mismatch: M0 {tuple(M0_full.shape)}, "
                         f"A {tuple(A_full.shape)}, eig {eig.shape}")
    if nb < nw:
        raise ValueError(f"num_bands ({nb}) < num_wann ({nw})")

    # -- windows (wannier90 defaults when not specified) ----------------------
    if dis_win_min is None:
        dis_win_min = float(eig.min())
    if dis_win_max is None:
        dis_win_max = float(eig.max())
    if froz_win_max is not None and froz_win_min is None:
        froz_win_min = dis_win_min
    win_np, froz_np = band_windows(eig, dis_win_min, dis_win_max,
                                   froz_win_min, froz_win_max, nw)
    win_t = torch.from_numpy(win_np).to(device)
    froz_t = torch.from_numpy(froz_np).to(device)
    ndimwin = win_np.sum(axis=1)
    ndimfroz = froz_np.sum(axis=1)

    winf = win_t.to(rdtype)                     # [nk, nb] 1/0 row masks
    frozf = froz_t.to(rdtype)
    freef = (win_t & ~froz_t).to(rdtype)
    exclf = (~win_t).to(rdtype)

    wbtot = float(wb[0].sum())

    # eigenvalue-slot masks over the ASCENDING eigh spectrum of the shifted Z:
    #   [... excluded @ -2*wbtot ... | free spectrum in [0, wbtot] | frozen @ +2*wbtot]
    # top nw slots  -> the new subspace columns
    # free slots    -> positions nb-nw .. nb-nf(k)-1, whose eigenvalue sum enters
    #                  Omega_I(i-1) (ragged in nf -> per-k mask, fully batched)
    pos = torch.arange(nb, device=device)
    nf_t = torch.from_numpy(ndimfroz).to(device)
    free_slots = ((pos[None, :] >= nb - nw)
                  & (pos[None, :] < (nb - nf_t)[:, None])).to(rdtype)  # [nk, nb]

    # window-masked overlaps: rows in-window at k, cols in-window at k+b
    M0w = M0_full * winf[:, None, :, None].to(dtype) \
                  * winf[kb_idx][:, :, None, :].to(dtype)
    del M0_full

    # window-masked projections and initial subspace
    A_w = A_full * winf[:, :, None].to(dtype)
    V = initial_subspace(A_w, win_t, froz_t)

    # diagonal level shifts for the Z eigh (Z spectrum lies in [0, wbtot] since
    # ||M|| <= 1 and w_b > 0): frozen +2*wbtot, excluded -2*wbtot.  The shifts are
    # baked into BOTH Z_in and Z_out, so the (1-mix)/mix combination preserves
    # them exactly and eigenvalue-slot selection stays valid across iterations.
    z_shift = torch.diag_embed((2.0 * wbtot) * (frozf - exclf)).to(dtype)
    # zero frozen+excluded rows/cols of Z (the Q_f . Q_f projection)
    qmask2 = (freef[:, :, None] * freef[:, None, :]).to(dtype)

    wbc = wb.to(dtype)

    def neighbour_overlap(Vcur):
        """MV[k,b] = M0w[k,b] @ V[k+b]  -> [nk, nnb, nb, nw]"""
        return M0w @ Vcur[kb_idx]

    def build_z(MV):
        """Z[k] = Q_f ( sum_b w_b MV MV^dag ) Q_f + shifts   (Hermitian, [nk,nb,nb])"""
        Zk = torch.einsum("kb,kbml,kbnl->kmn", wbc, MV, MV.conj())
        return Zk * qmask2 + z_shift

    def omega_i_of(Vcur, MVcur):
        """Omega_I with all states from the same iteration (w90 'Omega_I(i)')."""
        W = Vcur.conj().transpose(-1, -2)[:, None] @ MVcur   # [nk, nnb, nw, nw]
        s = torch.einsum("kb,kb->k", wb, (W.real**2 + W.imag**2).sum(dim=(-2, -1)))
        return (nw * wbtot - s).mean(), W

    t0 = time.perf_counter()
    MV = neighbour_overlap(V)
    Z_out = build_z(MV)
    Z_in = None

    hist_w1, hist_w, hist_delta = [], [], []
    delta_hist = torch.zeros(conv_window, device=device, dtype=rdtype) + np.inf
    converged = False
    it = 0
    W_red = None

    for it in range(1, num_iter + 1):
        # -- Z_in update (mixing scheme, w90 dis_extract 'extract_1') ---------
        Z_in = Z_out if it == 1 else mix_ratio * Z_out + (1.0 - mix_ratio) * Z_in

        # -- frozen contribution to Omega_I(i-1), previous-iteration states ---
        # (w90 'extract_2': sum over frozen rows m of |[M V_prev(k+b)]_{mn}|^2;
        #  MV still holds the previous iteration's neighbour overlaps here)
        froz_term = torch.einsum(
            "kb,kbmn,km->k", wb, (MV.real**2 + MV.imag**2), frozf)

        # -- new subspace: one batched Hermitian eigendecomposition per k ------
        evals, evecs = torch.linalg.eigh(Z_in)
        V = evecs[:, :, -nw:]
        V = _scrub_columns(V, evals[:, -nw:], win_t, froz_t, thresh=1.5 * wbtot)

        # Omega_I(i-1) = mean_k [ nw*wbtot - froz_term - sum(selected eigenvalues) ]
        eigsum = (evals * free_slots).sum(dim=1)
        womegai1 = (nw * wbtot - froz_term - eigsum).mean()

        # -- Omega_I(i) with the updated states (w90 'extract_4') -------------
        MV = neighbour_overlap(V)
        womegai, W_red = omega_i_of(V, MV)

        # -- Z_out for the next iteration's mixing ----------------------------
        Z_out = build_z(MV)

        # -- convergence: |delta| < tol over the last conv_window iterations --
        delta = womegai1 / womegai - 1.0
        hist_w1.append(womegai1)
        hist_w.append(womegai)
        hist_delta.append(delta)
        delta_hist = torch.roll(delta_hist, -1)
        delta_hist[-1] = delta

        if verbose:
            print(f"  dis iter {it:5d}  Omega_I(i-1)={womegai1.item():.8f}  "
                  f"Omega_I(i)={womegai.item():.8f}  delta={delta.item():.3e}")

        if it >= conv_window and (it % check_every == 0 or it == num_iter):
            if bool((delta_hist.abs() < conv_tol).all()):
                converged = True
                break

    dis_time = time.perf_counter() - t0

    # -- outputs for wannierise() ---------------------------------------------
    # M0_reduced = V^dag M0 V' is exactly the last W_red (computed with the final V).
    # U0 = polar(V^dag A_w): wannier90 internal_find_u (SVD unitary factor).
    caa = V.conj().transpose(-1, -2) @ A_w
    U0 = _polar_unitary(caa)

    h1 = torch.stack(hist_w1).cpu().numpy()
    h = torch.stack(hist_w).cpu().numpy()
    hd = torch.stack(hist_delta).cpu().numpy()
    history = [(i + 1, float(h1[i]), float(h[i]), float(hd[i]))
               for i in range(len(h))]

    if verbose:
        print(f"  dis {'converged' if converged else 'NOT converged'} after "
              f"{it} iterations ({dis_time:.1f} s), "
              f"final Omega_I = {float(h[-1]):.8f} Ang^2")

    return DisentangleResult(
        M0_reduced=W_red,
        U0=U0,
        V=V,
        omega_I=float(h[-1]),
        history=history,
        n_iter=it,
        converged=converged,
        ndimwin=ndimwin,
        ndimfroz=ndimfroz,
    )


def disentangle_case(case, device="cpu", num_iter=None, mix_ratio=None,
                     conv_tol=None, conv_window=None, check_every=1,
                     verbose=False):
    """Convenience wrapper: run disentangle() on an io_w90.W90Case.

    Window and iteration parameters are taken from the case's .win keys when
    present, falling back to wannier90 defaults (same resolution order as
    wannier90 parameters.F90).
    """
    if case.eig is None:
        raise ValueError("case has no band energies; load_case(..., eig_path=...)")
    w = case.win or {}
    kw = dict(
        dis_win_min=w.get("dis_win_min"),
        dis_win_max=w.get("dis_win_max"),
        froz_win_min=w.get("dis_froz_min"),
        froz_win_max=w.get("dis_froz_max"),
        num_iter=num_iter if num_iter is not None else (w.get("dis_num_iter") or 200),
        mix_ratio=mix_ratio if mix_ratio is not None else (w.get("dis_mix_ratio") or 0.5),
        conv_tol=conv_tol if conv_tol is not None else (w.get("dis_conv_tol") or 1e-10),
        conv_window=conv_window if conv_window is not None else (w.get("dis_conv_window") or 3),
    )
    return disentangle(
        torch.from_numpy(case.M),
        torch.from_numpy(case.A).to(torch.complex128),
        case.eig,
        torch.from_numpy(case.bvec),
        torch.from_numpy(case.wb),
        torch.from_numpy(case.kb_idx),
        device=device,
        check_every=check_every,
        verbose=verbose,
        **kw,
    )
