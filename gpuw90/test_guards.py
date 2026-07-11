"""
test_guards.py -- Unit tests for the entangled-input guards (Phase 2, issue #14).

Covers the failure mode that used to be silent: feeding an entangled case
(num_bands > num_wann) straight into wannierise() without disentanglement.

Run either way:
    python test_guards.py
    python -m pytest test_guards.py -q
"""

from __future__ import annotations

import os
import tempfile

import numpy as np
import torch

import io_w90
import wannierise as wz

NB, NW, NK, NNTOT = 3, 2, 1, 1


def _write_tiny_case(d):
    """Minimal entangled case (nb=3 > nw=2) in wannier90 file formats."""
    with open(os.path.join(d, "t.win"), "w") as fh:
        fh.write("num_wann = 2\nnum_bands = 3\nmp_grid = 1 1 1\n"
                 "begin unit_cell_cart\n1 0 0\n0 1 0\n0 0 1\nend unit_cell_cart\n"
                 "begin kpoints\n0.0 0.0 0.0\nend kpoints\n")
    with open(os.path.join(d, "t.mmn"), "w") as fh:
        fh.write(f"tiny\n{NB} {NK} {NNTOT}\n1 1 0 0 1\n")
        eye = np.eye(NB)
        for n in range(NB):          # column-major: m fast
            for m in range(NB):
                fh.write(f"{eye[m, n]:.1f} 0.0\n")
    with open(os.path.join(d, "t.amn"), "w") as fh:
        fh.write(f"tiny\n{NB} {NK} {NW}\n")
        for n in range(NW):
            for m in range(NB):
                v = 1.0 if m == n else 0.0
                fh.write(f"{m+1} {n+1} 1 {v:.1f} 0.0\n")
    with open(os.path.join(d, "t.eig"), "w") as fh:
        for k in range(NK):
            for m in range(NB):
                fh.write(f"{m+1} {k+1} {float(m):.6f}\n")


def test_load_case_rejects_entangled_by_default():
    with tempfile.TemporaryDirectory() as d:
        _write_tiny_case(d)
        try:
            io_w90.load_case(f"{d}/t.win", f"{d}/t.mmn", f"{d}/t.amn", "t")
        except ValueError as e:
            msg = str(e)
            assert "num_bands" in msg and "disentangle" in msg, msg
        else:
            raise AssertionError("load_case accepted nb > nw without "
                                 "allow_entangled -- silent-wrong-result guard "
                                 "is missing")


def test_load_case_entangled_needs_eig():
    with tempfile.TemporaryDirectory() as d:
        _write_tiny_case(d)
        try:
            io_w90.load_case(f"{d}/t.win", f"{d}/t.mmn", f"{d}/t.amn", "t",
                             allow_entangled=True)
        except ValueError as e:
            assert "eig_path" in str(e), str(e)
        else:
            raise AssertionError("entangled load without eig_path must fail")


def test_load_case_entangled_happy_path():
    with tempfile.TemporaryDirectory() as d:
        _write_tiny_case(d)
        case = io_w90.load_case(f"{d}/t.win", f"{d}/t.mmn", f"{d}/t.amn", "t",
                                eig_path=f"{d}/t.eig", allow_entangled=True)
        assert case.num_bands == NB and case.num_wann == NW
        assert case.eig.shape == (NK, NB)
        assert case.win["dis_win_min"] is None  # no dis keys in tiny .win


def test_wannierise_rejects_nonsquare_inputs():
    nk, nb, nw, nnb = 2, 3, 2, 1
    M0 = torch.eye(nb, dtype=torch.complex128).expand(nk, nnb, nb, nb)
    A = torch.randn(nk, nb, nw, dtype=torch.complex128)
    U0 = wz.lowdin_orthonormalize(A)             # [nk, nb, nw]: NOT square
    bvec = torch.zeros(nk, nnb, 3, dtype=torch.float64)
    wb = torch.ones(nk, nnb, dtype=torch.float64)
    kb = torch.zeros(nk, nnb, dtype=torch.long)
    try:
        wz.wannierise(M0, U0, bvec, wb, kb, num_iter=1)
    except ValueError as e:
        assert "disentangle" in str(e), str(e)
    else:
        raise AssertionError("wannierise accepted full-band M0 with a "
                             "rectangular U0 -- this used to return silently "
                             "wrong centres")


def test_wannierise_rejects_mismatched_nw():
    nk, nnb = 2, 1
    M0 = torch.eye(3, dtype=torch.complex128).expand(nk, nnb, 3, 3)
    U0 = torch.eye(2, dtype=torch.complex128).expand(nk, 2, 2)
    bvec = torch.zeros(nk, nnb, 3, dtype=torch.float64)
    wb = torch.ones(nk, nnb, dtype=torch.float64)
    kb = torch.zeros(nk, nnb, dtype=torch.long)
    try:
        wz.wannierise(M0, U0, bvec, wb, kb, num_iter=1)
    except ValueError as e:
        assert "mismatch" in str(e), str(e)
    else:
        raise AssertionError("wannierise accepted M0 nw=3 with U0 nw=2")


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"{len(fns)} guard tests passed")
