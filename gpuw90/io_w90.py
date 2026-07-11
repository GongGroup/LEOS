"""
io_w90.py -- Parsers for wannier90 standard files used by the GPU wannierise tool.

Reads the ab-initio overlap (.mmn) and initial-projection (.amn) matrices plus the
input file (.win) needed to reconstruct the finite-difference b-vectors and shell
weights.  Everything here is pure standard-library + numpy; the torch side lives in
wannierise.py.

Conventions (all matching wannier90 3.1):
  * .mmn stores M^{k,b}_{mn} = <u_{m,k} | u_{n,k+b}> with the band index m as the
    FAST (column-major) index inside each (k, b) block.
  * .amn stores A^{k}_{mn} = <psi_{m,k} | g_n> (band m, trial orbital n).
  * b-vectors are k+b - k in fractional reciprocal coordinates times the reciprocal
    lattice (Ang^-1, including the 2*pi factor), reconstructed from .win + the
    per-block (ik, ik_b, G) records that .mmn carries in its headers.  This means we
    do NOT need a separate .nnkp file for the isolated-band (wannierise-only) stage.
  * shell weights w_b solve the B1 completeness relation
    sum_b w_b b_a b_b = delta_ab   [Eq. (B1), PRB 56, 12847 (1997)],
    reconstructed here via the wannier90 kmesh shell least-squares rule.

Units: CODATA-2006 Bohr radius, matching the reference .wout header.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np

BOHR_TO_ANG = 0.52917720859  # CODATA 2006, matches wannier90 3.1 / reference .wout


# --------------------------------------------------------------------------- #
# Low-level file parsers
# --------------------------------------------------------------------------- #
def read_mmn(path):
    """Parse a wannier90 .mmn overlap file.

    Returns
    -------
    num_bands, num_kpts, nntot : int
    kb_idx : int ndarray [num_kpts, nntot]      0-based neighbour k index (ik_b)
    gvec   : int ndarray [num_kpts, nntot, 3]    reciprocal-lattice shift G
    M      : complex128 ndarray [num_kpts, nntot, num_bands, num_bands]
             M[k, b, m, n] = <u_{m,k} | u_{n, k+b}>
    """
    with open(path, "r") as fh:
        fh.readline()  # comment line
        num_bands, num_kpts, nntot = (int(x) for x in fh.readline().split())

        kb_idx = np.zeros((num_kpts, nntot), dtype=np.int64)
        gvec = np.zeros((num_kpts, nntot, 3), dtype=np.int64)
        M = np.zeros((num_kpts, nntot, num_bands, num_bands), dtype=np.complex128)

        # per-k running slot counter, so we tolerate any block ordering
        slot = np.zeros(num_kpts, dtype=np.int64)
        nblocks = num_kpts * nntot
        for _ in range(nblocks):
            hdr = fh.readline().split()
            ik = int(hdr[0]) - 1
            ikb = int(hdr[1]) - 1
            g = (int(hdr[2]), int(hdr[3]), int(hdr[4]))
            b = slot[ik]
            slot[ik] += 1
            kb_idx[ik, b] = ikb
            gvec[ik, b] = g
            # data block: num_bands*num_bands lines, band index m fastest (col-major)
            block = np.empty((num_bands * num_bands, 2), dtype=np.float64)
            for i in range(num_bands * num_bands):
                re_im = fh.readline().split()
                block[i, 0] = float(re_im[0])
                block[i, 1] = float(re_im[1])
            vals = block[:, 0] + 1j * block[:, 1]
            # column-major (m fast) -> reshape with Fortran order into [m, n]
            M[ik, b] = vals.reshape((num_bands, num_bands), order="F")

    return num_bands, num_kpts, nntot, kb_idx, gvec, M


def read_amn(path):
    """Parse a wannier90 .amn projection file.

    Returns
    -------
    num_bands, num_wann, num_kpts : int
    A : complex128 ndarray [num_kpts, num_bands, num_wann]
        A[k, m, n] = <psi_{m,k} | g_n>
    """
    with open(path, "r") as fh:
        fh.readline()  # comment line
        num_bands, num_kpts, num_wann = (int(x) for x in fh.readline().split())
        A = np.zeros((num_kpts, num_bands, num_wann), dtype=np.complex128)
        for line in fh:
            s = line.split()
            if len(s) < 5:
                continue
            m = int(s[0]) - 1
            n = int(s[1]) - 1
            k = int(s[2]) - 1
            A[k, m, n] = float(s[3]) + 1j * float(s[4])
    return num_bands, num_wann, num_kpts, A


def read_eig(path):
    """Parse a wannier90 .eig file (band energies for disentanglement windows).

    Format: one line per (band, k-point) pair, ``m  k  E`` with the band index m
    as the fast index (same convention as wannier90/VASP LWRITE_MMN_AMN output).

    Returns
    -------
    eig : float64 ndarray [num_kpts, num_bands]   energies in eV
    """
    data = np.loadtxt(path, dtype=np.float64)
    if data.ndim == 1:
        data = data[None, :]
    m_idx = data[:, 0].astype(np.int64)
    k_idx = data[:, 1].astype(np.int64)
    nb = int(m_idx.max())
    nk = int(k_idx.max())
    if data.shape[0] != nb * nk:
        raise ValueError(f".eig has {data.shape[0]} rows, expected nb*nk = {nb*nk}")
    eig = np.full((nk, nb), np.nan)
    eig[k_idx - 1, m_idx - 1] = data[:, 2]
    if np.isnan(eig).any():
        raise ValueError(".eig is missing (band, kpt) entries")
    # wannier90 assumes eigenvalues sorted bottom-up per k (dis_windows relies on
    # contiguity of the energy windows); enforce it here rather than mis-slice.
    if np.any(np.diff(eig, axis=1) < -1e-8):
        raise ValueError(".eig energies are not sorted ascending within a k-point")
    return eig


def read_nnkp(path):
    """Parse a wannier90 seedname.nnkp file (kpoints, lattices, neighbour list).

    Provided as an alternative source of kpoints/recip_lattice when a .win does not
    list explicit k-points, and as an independent cross-check of the b-vector neighbour
    list.  Returns a dict with keys real_lattice, recip_lattice, kpoints, nnkpts
    (the last as [num_kpts, nntot, 4] int: ik_b(0-based), g1, g2, g3).
    """
    with open(path, "r") as fh:
        lines = [_strip_comment(l) for l in fh.readlines()]
    out = {}

    def block(name):
        for i, ln in enumerate(lines):
            if ln.lower().startswith(f"begin {name}"):
                j = i + 1
                rows = []
                while j < len(lines) and not lines[j].lower().startswith(f"end {name}"):
                    if lines[j]:
                        rows.append(lines[j].split())
                    j += 1
                return rows
        return None

    rl = block("real_lattice")
    if rl:
        out["real_lattice"] = np.array([[float(x) for x in r[:3]] for r in rl])
    gl = block("recip_lattice")
    if gl:
        out["recip_lattice"] = np.array([[float(x) for x in r[:3]] for r in gl])
    kp = block("kpoints")
    if kp:
        nk = int(kp[0][0])
        out["kpoints"] = np.array([[float(x) for x in kp[1 + i][:3]]
                                   for i in range(nk)])
    nn = block("nnkpts")
    if nn:
        nntot = int(nn[0][0])
        nk = (len(nn) - 1) // nntot
        arr = np.zeros((nk, nntot, 4), dtype=np.int64)
        slot = np.zeros(nk, dtype=np.int64)
        for r in nn[1:]:
            ik = int(r[0]) - 1
            b = slot[ik]
            slot[ik] += 1
            arr[ik, b] = [int(r[1]) - 1, int(r[2]), int(r[3]), int(r[4])]
        out["nnkpts"] = arr
    return out


def _strip_comment(line: str) -> str:
    for c in ("!", "#"):
        idx = line.find(c)
        if idx >= 0:
            line = line[:idx]
    return line.strip()


def read_win(path):
    """Parse the subset of a wannier90 .win file we need.

    Returns a dict with keys:
      num_wann, num_bands, mp_grid (list[3]),
      real_lattice  [3,3] rows a_i in Ang,
      recip_lattice [3,3] rows b_i in Ang^-1 (includes 2*pi),
      kpoints       [num_kpts,3] fractional.
    """
    with open(path, "r") as fh:
        raw = fh.readlines()
    lines = [_strip_comment(l) for l in raw]

    out = {"num_wann": None, "num_bands": None, "mp_grid": None}
    # disentanglement keys (None = not present; defaults are wannier90's and are
    # resolved inside disentangle(): win = full eigval range, froz_min = win_min,
    # num_iter 200, mix 0.5, conv_tol 1e-10, conv_window 3)
    for key in ("dis_win_min", "dis_win_max", "dis_froz_min", "dis_froz_max",
                "dis_mix_ratio", "dis_conv_tol", "dis_num_iter", "dis_conv_window"):
        out[key] = None
    cell = None
    cell_unit = "ang"
    kpts = None

    i = 0
    while i < len(lines):
        line = lines[i]
        if not line:
            i += 1
            continue
        low = line.lower()

        m = re.match(r"num_wann\s*[=:]\s*(\d+)", low)
        if m:
            out["num_wann"] = int(m.group(1))
            i += 1
            continue
        m = re.match(r"num_bands\s*[=:]\s*(\d+)", low)
        if m:
            out["num_bands"] = int(m.group(1))
            i += 1
            continue
        m = re.match(r"mp_grid\s*[=:]\s*(\d+)\s+(\d+)\s+(\d+)", low)
        if m:
            out["mp_grid"] = [int(m.group(1)), int(m.group(2)), int(m.group(3))]
            i += 1
            continue
        m = re.match(r"(dis_win_min|dis_win_max|dis_froz_min|dis_froz_max|"
                     r"dis_mix_ratio|dis_conv_tol)\s*[=:]\s*([-+0-9.ed]+)", low)
        if m:
            out[m.group(1)] = float(m.group(2).replace("d", "e"))
            i += 1
            continue
        m = re.match(r"(dis_num_iter|dis_conv_window)\s*[=:]\s*(\d+)", low)
        if m:
            out[m.group(1)] = int(m.group(2))
            i += 1
            continue

        if low.startswith("begin unit_cell_cart"):
            j = i + 1
            rows = []
            unit = "ang"
            while j < len(lines) and not lines[j].lower().startswith("end unit_cell_cart"):
                tok = lines[j].split()
                if len(tok) == 1 and tok[0].lower() in ("bohr", "ang", "angstrom"):
                    unit = "bohr" if tok[0].lower() == "bohr" else "ang"
                elif len(tok) >= 3:
                    rows.append([float(tok[0]), float(tok[1]), float(tok[2])])
                j += 1
            cell = np.array(rows, dtype=np.float64)
            cell_unit = unit
            i = j + 1
            continue

        if low.startswith("begin kpoints"):
            j = i + 1
            rows = []
            while j < len(lines) and not lines[j].lower().startswith("end kpoints"):
                tok = lines[j].split()
                if len(tok) >= 3:
                    rows.append([float(tok[0]), float(tok[1]), float(tok[2])])
                j += 1
            kpts = np.array(rows, dtype=np.float64)
            i = j + 1
            continue

        i += 1

    if out["num_bands"] is None:
        out["num_bands"] = out["num_wann"]
    if cell is None:
        raise ValueError("unit_cell_cart block not found in .win")
    if cell_unit == "bohr":
        cell = cell * BOHR_TO_ANG
    out["real_lattice"] = cell
    out["recip_lattice"] = 2.0 * np.pi * np.linalg.inv(cell).T
    out["kpoints"] = kpts
    return out


# --------------------------------------------------------------------------- #
# b-vectors and shell weights (wannier90 kmesh rule)
# --------------------------------------------------------------------------- #
def compute_bvectors(kpoints, recip_lattice, kb_idx, gvec):
    """Cartesian b-vectors (Ang^-1) for every (k, neighbour) pair.

    b_cart[k,b] = ((kpt[ik_b] + G) - kpt[ik]) @ recip_lattice
    """
    nk, nnb = kb_idx.shape
    frac = kpoints[kb_idx] + gvec - kpoints[:, None, :]   # [nk, nnb, 3] fractional
    bvec = frac @ recip_lattice                            # [nk, nnb, 3] Ang^-1
    return bvec


def compute_shell_weights(bvec, tol=1e-5):
    """Shell weights w_b satisfying the B1 completeness relation.

    Groups the b-vectors around k=0 into shells by |b|, then solves the least-squares
    system  sum_shell w_s * (sum_{b in shell} b_a b_b) = delta_ab  for the six
    independent symmetric components.  All k share the same shell structure on a
    uniform mesh, so weights are assigned per (k,b) by matching |b|.

    Returns
    -------
    wb  : float ndarray [nk, nnb]     weight for each (k,b)
    info: dict with shell radii, per-shell weights and the B1 residual (diagnostics)
    """
    nk, nnb, _ = bvec.shape
    b0 = bvec[0]                       # [nnb, 3]
    blen = np.linalg.norm(b0, axis=1)  # [nnb]

    # group into shells by |b|
    order = np.argsort(blen)
    shell_radii = []
    shell_of = np.full(nnb, -1, dtype=np.int64)
    for idx in order:
        placed = False
        for s, r in enumerate(shell_radii):
            if abs(blen[idx] - r) < tol:
                shell_of[idx] = s
                placed = True
                break
        if not placed:
            shell_radii.append(blen[idx])
            shell_of[idx] = len(shell_radii) - 1
    nshell = len(shell_radii)

    # least-squares matrix: rows = (xx, yy, zz, xy, xz, yz), cols = shells
    def sym6(vsum):
        return np.array([vsum[0, 0], vsum[1, 1], vsum[2, 2],
                         vsum[0, 1], vsum[0, 2], vsum[1, 2]])

    A = np.zeros((6, nshell))
    for s in range(nshell):
        vs = b0[shell_of == s]
        outer = np.einsum("bi,bj->ij", vs, vs)   # sum_{b in shell} b_i b_j
        A[:, s] = sym6(outer)
    q = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    w_shell, *_ = np.linalg.lstsq(A, q, rcond=None)

    # assign per (k,b): each k has the same shell structure; match by |b|
    wb = np.zeros((nk, nnb))
    for k in range(nk):
        bl = np.linalg.norm(bvec[k], axis=1)
        for b in range(nnb):
            s = int(np.argmin([abs(bl[b] - r) for r in shell_radii]))
            wb[k, b] = w_shell[s]

    # B1 residual on k=0 as a diagnostic
    B1 = np.einsum("b,bi,bj->ij", wb[0], b0, b0)
    residual = float(np.max(np.abs(B1 - np.eye(3))))

    info = {
        "shell_radii": [float(r) for r in shell_radii],
        "w_shell": [float(w) for w in w_shell],
        "b1_residual": residual,
        "wbtot": float(wb[0].sum()),
    }
    return wb, info


# --------------------------------------------------------------------------- #
# Bundled loader
# --------------------------------------------------------------------------- #
@dataclass
class W90Case:
    seedname: str
    num_bands: int
    num_wann: int
    num_kpts: int
    nntot: int
    M: np.ndarray          # [nk, nnb, nb, nb] complex128
    A: np.ndarray          # [nk, nb, nw]      complex128
    kb_idx: np.ndarray     # [nk, nnb] int
    bvec: np.ndarray       # [nk, nnb, 3] float64 (Ang^-1)
    wb: np.ndarray         # [nk, nnb] float64
    kpoints: np.ndarray    # [nk, 3] fractional
    real_lattice: np.ndarray   # [3,3] Ang
    recip_lattice: np.ndarray  # [3,3] Ang^-1
    kmesh_info: dict
    win: dict = None       # parsed .win keys (incl. dis_* when present)
    eig: np.ndarray = None  # [nk, nb] band energies (eV), only when eig_path given


def load_case(win_path, mmn_path, amn_path, seedname="case", nnkp_path=None,
              eig_path=None, allow_entangled=False):
    """Read a full wannier90 case and reconstruct b-vectors + shell weights.

    b-vectors need kpoints (fractional) + recip_lattice.  Both are taken from the .win
    (self-contained, production LEOS path).  If the .win has no explicit `kpoints`
    block, pass `nnkp_path` to source them (and the reciprocal lattice) from a
    seedname.nnkp instead.

    Entangled guard: if num_bands > num_wann the raw M/A must NOT be fed straight
    into wannierise() (it would silently produce wrong centres).  Loading such a
    case requires `allow_entangled=True` plus an `eig_path`, and the caller must
    run disentangle.disentangle() first.
    """
    win = read_win(win_path)
    nb_m, nk_m, nnb, kb_idx, gvec, M = read_mmn(mmn_path)
    nb_a, nw_a, nk_a, A = read_amn(amn_path)

    assert nb_m == nb_a, f"num_bands mismatch mmn/amn: {nb_m} vs {nb_a}"
    assert nk_m == nk_a, f"num_kpts mismatch mmn/amn: {nk_m} vs {nk_a}"

    if nb_m > nw_a and not allow_entangled:
        raise ValueError(
            f"num_bands ({nb_m}) > num_wann ({nw_a}): entangled bands. "
            "wannierise() alone would return silently wrong centres. "
            "Load with allow_entangled=True (and eig_path=...) and run "
            "disentangle.disentangle() first; feed its M0_reduced/U0 to wannierise()."
        )
    if nb_m > nw_a and eig_path is None:
        raise ValueError(
            "entangled case needs band energies to place the outer/frozen windows: "
            "pass eig_path='seedname.eig'"
        )

    eig = None
    if eig_path is not None:
        eig = read_eig(eig_path)
        if eig.shape != (nk_m, nb_m):
            raise ValueError(
                f".eig shape {eig.shape} does not match .mmn (nk={nk_m}, nb={nb_m})"
            )

    kpoints = win["kpoints"]
    recip = win["recip_lattice"]
    real = win["real_lattice"]
    if kpoints is None:
        if nnkp_path is None:
            raise ValueError("no explicit kpoints in .win; provide nnkp_path")
        nnkp = read_nnkp(nnkp_path)
        kpoints = nnkp["kpoints"]
        recip = nnkp.get("recip_lattice", recip)
        real = nnkp.get("real_lattice", real)

    bvec = compute_bvectors(kpoints, recip, kb_idx, gvec)
    wb, kmesh_info = compute_shell_weights(bvec)

    return W90Case(
        seedname=seedname,
        num_bands=nb_m,
        num_wann=nw_a,
        num_kpts=nk_m,
        nntot=nnb,
        M=M,
        A=A,
        kb_idx=kb_idx,
        bvec=bvec,
        wb=wb,
        kpoints=kpoints,
        real_lattice=real,
        recip_lattice=recip,
        kmesh_info=kmesh_info,
        win=win,
        eig=eig,
    )


# --------------------------------------------------------------------------- #
# Reference .wout parser (for the correctness check)
# --------------------------------------------------------------------------- #
def parse_wout_final_state(path):
    """Extract the 'Final State' WF centres/spreads and Omega decomposition.

    Returns dict: centres [nw,3], spreads [nw], omega {I,D,OD,total}.
    Also returns the 'Initial State' block if present (same structure) under 'initial'.
    """
    with open(path, "r") as fh:
        text = fh.read()
    lines = text.splitlines()

    def grab_block(tag):
        centres, spreads = [], []
        for i, ln in enumerate(lines):
            if ln.strip() == tag:
                j = i + 1
                while j < len(lines):
                    m = re.search(
                        r"WF centre and spread\s+\d+\s+\(\s*([-\d.]+),\s*([-\d.]+),"
                        r"\s*([-\d.]+)\s*\)\s+([-\d.]+)",
                        lines[j],
                    )
                    if m:
                        centres.append([float(m.group(1)), float(m.group(2)),
                                        float(m.group(3))])
                        spreads.append(float(m.group(4)))
                        j += 1
                    elif "Sum of centres" in lines[j] or lines[j].strip() == "":
                        break
                    else:
                        j += 1
                break
        return np.array(centres), np.array(spreads)

    centres, spreads = grab_block("Final State")
    init_c, init_s = grab_block("Initial State")

    omega = {}
    for key, pat in (("I", r"Omega I\s*=\s*([-\d.]+)"),
                     ("D", r"Omega D\s*=\s*([-\d.]+)"),
                     ("OD", r"Omega OD\s*=\s*([-\d.]+)"),
                     ("total", r"Omega Total\s*=\s*([-\d.]+)")):
        m = re.search(pat, text)
        omega[key] = float(m.group(1)) if m else None

    out = {"centres": centres, "spreads": spreads, "omega": omega}
    if len(init_c):
        out["initial"] = {"centres": init_c, "spreads": init_s}
    return out


def parse_wout_dis_info(path):
    """Extract disentanglement results from a reference .wout.

    Returns dict: omega_I (the 'Final Omega_I' after dis_extract, Ang^2),
    converged (bool), n_iter (number of DIS iterations run).
    """
    with open(path, "r") as fh:
        text = fh.read()
    m = re.search(r"Final Omega_I\s+([-\d.]+)", text)
    omega_i = float(m.group(1)) if m else None
    converged = "Disentanglement convergence criteria satisfied" in text
    iters = re.findall(r"^\s+(\d+)\s+[-\d.]+\s+[-\d.]+\s+[-+\dEe.]+\s+[-\d.]+\s+<-- DIS",
                       text, flags=re.M)
    n_iter = int(iters[-1]) if iters else None
    return {"omega_I": omega_i, "converged": converged, "n_iter": n_iter}


if __name__ == "__main__":
    import sys
    d = sys.argv[1] if len(sys.argv) > 1 else "testdata/gaas"
    case = load_case(f"{d}/gaas.win", f"{d}/gaas.mmn", f"{d}/gaas.amn", "gaas")
    print(f"num_bands={case.num_bands} num_wann={case.num_wann} "
          f"num_kpts={case.num_kpts} nntot={case.nntot}")
    print("kmesh:", case.kmesh_info)
    print("b0[0..2]:\n", case.bvec[0, :3])
