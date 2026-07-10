"""
wannierise.py -- Batched Marzari-Vanderbilt spread minimisation in PyTorch.

Isolated-band (wannierise-only) stage: given the ab-initio overlaps M0^{k,b} and an
initial gauge U0^k (from Loewdin-orthonormalised .amn projections), minimise the
gauge-dependent part of the MV spread functional by steepest descent on the unitary
manifold.  The k dimension is the batch dimension, so every heavy op is a batched
matmul / batched Hermitian eigendecomposition -- the same code runs on CPU or CUDA
by flipping the `device` argument (numerics are identical; only the reduction order
differs at the ~1e-14 level).

Math (all matching wannier90 3.1 wannierise.F90 / Marzari-Vanderbilt PRB 56, 12847):

  rotated overlaps      M^{k,b} = U_k^dagger M0^{k,b} U_{k+b}
  Wannier centres       r_n = -(1/N) sum_{k,b} w_b b Im ln M^{k,b}_nn
  gradient              G_k = 4 sum_b w_b ( A[R^{k,b}] - S[T^{k,b}] )      (anti-Herm.)
     with  R_mn = M_mn conj(M_nn),   T_mn = (M_mn / M_nn) q_n,
           q_n  = Im ln M_nn + b . r_n,
           A[B] = (B - B^dagger)/2,   S[B] = (B + B^dagger)/(2i)
  step                  dW_k = -(alpha / (4 w_tot)) G_k
  update                U_k <- U_k exp(dW_k)        via i*dW Hermitian eigendecomp

Numerical-pitfall handling is documented inline at each site (ln branch cut, weights,
degenerate-eigenvector phase gauge).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch


# --------------------------------------------------------------------------- #
# Initial gauge
# --------------------------------------------------------------------------- #
def lowdin_orthonormalize(A):
    """Loewdin (symmetric) orthonormalisation of the projection matrices.

    U_k = A_k (A_k^dagger A_k)^{-1/2}.  For the isolated case A_k is [nw, nw] and the
    result is the initial unitary gauge that wannier90 reports as 'Initial State'.

    (A^dagger A)^{-1/2} is built from a Hermitian eigendecomposition; it is invariant
    under the eigenvector phase/degeneracy gauge, so no phase convention is needed.
    """
    S = A.conj().transpose(-1, -2) @ A               # [nk, nw, nw] Hermitian PD
    evals, evecs = torch.linalg.eigh(S)
    inv_sqrt = evecs @ torch.diag_embed(
        (1.0 / torch.sqrt(evals)).to(A.dtype)
    ) @ evecs.conj().transpose(-1, -2)
    return A @ inv_sqrt


# --------------------------------------------------------------------------- #
# Core quantities
# --------------------------------------------------------------------------- #
def _rotate_M(M0, U, kb_idx):
    """M^{k,b} = U_k^dagger M0^{k,b} U_{k+b}, batched over (k, b)."""
    Ul = U.conj().transpose(-1, -2)          # [nk, nw, nw]
    Ur = U[kb_idx]                           # [nk, nnb, nw, nw]
    return Ul[:, None] @ M0 @ Ur             # [nk, nnb, nw, nw]


def _centres(M, bvec, wb, nk):
    """Wannier centres r_n (Ang) and the diagonal Im ln M_nn.

    Im ln M_nn uses the principal branch atan2(Im, Re) in (-pi, pi], identical to
    wannier90's aimag(log(...)).  For well-localised insulator WFs M_nn ~ 1, so the
    phase is small and branch-safe; the same principal value is used on CPU and GPU
    so there is no cross-device branch ambiguity.
    """
    Mnn = torch.diagonal(M, dim1=-2, dim2=-1)          # [nk, nnb, nw]
    lnim = torch.atan2(Mnn.imag, Mnn.real)             # [nk, nnb, nw] principal branch
    # r_n[c] = -(1/nk) sum_{k,b} w_b b_c Im ln M_nn
    centres = -(1.0 / nk) * torch.einsum("kb,kbc,kbn->nc", wb, bvec, lnim)
    return centres, Mnn, lnim


def spread_decomposition(M, bvec, wb, nk):
    """Omega decomposition and per-WF spread, matching wannier90's report.

    Omega_I  = (1/N) sum_kb w_b (nw - sum_n |M_nn|^2)          (gauge invariant)
    Omega_OD = (1/N) sum_kb w_b sum_{m!=n} |M_mn|^2
    Omega_D  = (1/N) sum_kb w_b sum_n (Im ln M_nn + b.r_n)^2
    spread_n = <r^2>_n - |r_n|^2
    """
    nw = M.shape[-1]
    centres, Mnn, lnim = _centres(M, bvec, wb, nk)
    absM2 = (M.real**2 + M.imag**2)                    # [nk, nnb, nw, nw]
    absMnn2 = (Mnn.real**2 + Mnn.imag**2)              # [nk, nnb, nw]

    sum_all = absM2.sum(dim=(-2, -1))                  # [nk, nnb] = sum_mn |M_mn|^2
    # Omega_I uses the FULL sum over m,n (invariant part); Omega_OD is the off-diagonal
    # remainder.  (A common slip is to put sum_n|M_nn|^2 into Omega_I, which silently
    # folds Omega_OD into Omega_I and then double-counts it.)
    omega_I = (1.0 / nk) * torch.einsum("kb,kb->", wb, (nw - sum_all))
    omega_OD = (1.0 / nk) * torch.einsum("kb,kb->", wb, (sum_all - absMnn2.sum(-1)))

    bdotr = torch.einsum("kbc,nc->kbn", bvec, centres)
    q = lnim + bdotr                                    # [nk, nnb, nw]
    omega_D = (1.0 / nk) * torch.einsum("kb,kbn->", wb, q**2)

    # per-WF: <r^2>_n = (1/N) sum_kb w_b [ (1-|Mnn|^2) + (Im ln Mnn)^2 ]
    r2 = (1.0 / nk) * torch.einsum("kb,kbn->n", wb, (1.0 - absMnn2) + lnim**2)
    spread_n = r2 - (centres**2).sum(-1)

    omega = {
        "I": float(omega_I),
        "OD": float(omega_OD),
        "D": float(omega_D),
        "total": float(omega_I + omega_OD + omega_D),
    }
    return centres, spread_n, omega


def _gradient(M, bvec, wb, centres, Mnn, lnim):
    """True MV gradient g_k = dOmega/dW_k; anti-Hermitian [nk, nw, nw].

    g_k = -(1/N) * 4 sum_b w_b ( A[R^{k,b}] - S[T^{k,b}] )
    The 1/N prefactor comes from Omega's own (1/N) sum over k, and the overall sign was
    pinned by a finite-difference directional-derivative test (fd_check): with this
    definition, Re Tr[g_k^dagger X_k] equals dOmega/deps along U exp(eps X), so plain
    descent dW = -eff * g decreases Omega.
    """
    nk = M.shape[0]
    bdotr = torch.einsum("kbc,nc->kbn", bvec, centres)
    q = lnim + bdotr                                    # [nk, nnb, nw] real

    Mnn_c = Mnn.conj()
    R = M * Mnn_c[:, :, None, :]                        # R_mn = M_mn conj(M_nn)
    Tfac = (q / Mnn)                                    # [nk, nnb, nw] complex
    T = M * Tfac[:, :, None, :]                         # T_mn = (M_mn/M_nn) q_n

    Rh = R.conj().transpose(-1, -2)
    Th = T.conj().transpose(-1, -2)
    A_R = (R - Rh) * 0.5                                # anti-Hermitian part
    S_T = (T + Th) / (2.0j)                             # Hermitian part / i
    Gkb = A_R - S_T                                     # [nk, nnb, nw, nw]
    Graw = 4.0 * torch.einsum("kb,kbmn->kmn", wb.to(Gkb.dtype), Gkb)
    return -(1.0 / nk) * Graw


def _expm_antihermitian(dW):
    """exp(dW) for anti-Hermitian dW via a Hermitian eigendecomposition.

    dW anti-Hermitian  ->  H = i dW is Hermitian  ->  exp(dW) = exp(-iH) = V e^{-i l} V^d.
    This is exactly wannier90's zheev-based exponentiation.  Because exp is applied as
    a function of the operator (V f(l) V^dagger), the result is invariant to the
    eigenvector phase and to rotations inside degenerate eigenspaces -- so the known
    'degenerate eigh phase gauge' pitfall does NOT affect the update.
    """
    H = 1j * dW
    evals, evecs = torch.linalg.eigh(H)                 # H = V diag(evals) V^dagger
    phase = torch.exp(-1j * evals).to(dW.dtype)
    return evecs @ torch.diag_embed(phase) @ evecs.conj().transpose(-1, -2)


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
@dataclass
class WannieriseResult:
    U: torch.Tensor           # [nk, nw, nw] final gauge
    centres: np.ndarray       # [nw, 3] Ang
    spreads: np.ndarray       # [nw]
    omega: dict
    history: list             # list of (iter, omega_total, rms_grad)
    n_iter: int
    converged: bool


def wannierise(
    M0,
    U0,
    bvec,
    wb,
    kb_idx,
    num_iter=1000,
    step=1.0,
    conv_tol=1e-10,
    conv_window=3,
    device="cpu",
    dtype=torch.complex128,
    verbose=False,
):
    """Batched MV steepest-descent wannierisation.

    Parameters
    ----------
    M0      : [nk, nnb, nw, nw] complex  ab-initio overlaps (isolated case nb == nw)
    U0      : [nk, nw, nw]      complex  initial gauge (Loewdin-orthonormalised .amn)
    bvec    : [nk, nnb, 3]      real     Cartesian b-vectors (Ang^-1)
    wb      : [nk, nnb]         real     shell weights
    kb_idx  : [nk, nnb]         long     neighbour k index
    step    : SD step scale alpha (effective step = alpha / (4 w_tot))
    conv_tol, conv_window : stop when |dOmega| < conv_tol for `conv_window` iters
                            (matches wannier90 conv_tol / conv_window semantics)

    Returns WannieriseResult.
    """
    rdtype = torch.float64 if dtype == torch.complex128 else torch.float32
    M0 = M0.to(device=device, dtype=dtype)
    U = U0.to(device=device, dtype=dtype).clone()
    bvec = bvec.to(device=device, dtype=rdtype)
    wb = wb.to(device=device, dtype=rdtype)
    kb_idx = kb_idx.to(device=device, dtype=torch.long)

    nk = M0.shape[0]
    wbtot = float(wb[0].sum())
    # classic MV steepest-descent step: dW = -(1/(4 w_tot)) * G_raw, and since the true
    # gradient g = -(1/N) G_raw, that is dW = -(N/(4 w_tot)) * g.  step==alpha scales it.
    eff_step = step * nk / (4.0 * wbtot)

    history = []
    prev_omega = None
    n_below = 0
    converged = False
    it = 0

    for it in range(1, num_iter + 1):
        M = _rotate_M(M0, U, kb_idx)
        centres, Mnn, lnim = _centres(M, bvec, wb, nk)
        G = _gradient(M, bvec, wb, centres, Mnn, lnim)

        # anti-Hermitian gradient; RMS over all k,m,n for reporting/convergence
        rms_grad = float(torch.sqrt((G.real**2 + G.imag**2).mean()))

        # steepest-descent step on the unitary manifold (minimise -> -G direction)
        dW = -(eff_step) * G
        U = U @ _expm_antihermitian(dW)

        # spread AFTER the update (matches wannier90's per-cycle report)
        Mn = _rotate_M(M0, U, kb_idx)
        _, _, omega = spread_decomposition(Mn, bvec, wb, nk)
        omega_tot = omega["total"]
        history.append((it, omega_tot, rms_grad))
        if verbose and (it <= 3 or it % 50 == 0):
            print(f"  iter {it:5d}  Omega={omega_tot:.10f}  rms_grad={rms_grad:.3e}")

        if prev_omega is not None:
            if abs(omega_tot - prev_omega) < conv_tol:
                n_below += 1
                if n_below >= conv_window:
                    converged = True
                    prev_omega = omega_tot
                    break
            else:
                n_below = 0
        prev_omega = omega_tot

    # final state from converged U
    Mf = _rotate_M(M0, U, kb_idx)
    centres, spread_n, omega = spread_decomposition(Mf, bvec, wb, nk)
    return WannieriseResult(
        U=U,
        centres=centres.detach().cpu().numpy(),
        spreads=spread_n.detach().cpu().numpy(),
        omega=omega,
        history=history,
        n_iter=it,
        converged=converged,
    )


def initial_state(M0, U0, bvec, wb, kb_idx, device="cpu", dtype=torch.complex128):
    """Spread/centres of the initial gauge (matches wannier90 'Initial State')."""
    rdtype = torch.float64 if dtype == torch.complex128 else torch.float32
    M0 = M0.to(device=device, dtype=dtype)
    U0 = U0.to(device=device, dtype=dtype)
    bvec = bvec.to(device=device, dtype=rdtype)
    wb = wb.to(device=device, dtype=rdtype)
    kb_idx = kb_idx.to(device=device, dtype=torch.long)
    nk = M0.shape[0]
    M = _rotate_M(M0, U0, kb_idx)
    centres, spread_n, omega = spread_decomposition(M, bvec, wb, nk)
    return centres.detach().cpu().numpy(), spread_n.detach().cpu().numpy(), omega
