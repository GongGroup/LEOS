"""Finite-difference self-check of the MV gradient sign / prefactor.

Compares the analytic directional derivative Re Tr[g_k^dagger X_k] against the
central finite difference [Omega(U e^{eps X}) - Omega(U e^{-eps X})] / (2 eps) for a
random anti-Hermitian direction X.  Expected ratio ~ +1.0000 across eps: it confirms
that wannierise._gradient returns the true dOmega/dW, so plain descent dW = -eff*g
decreases Omega.  Run this after any change to the gradient.

    python fd_check.py [testdata_dir] [seed]
"""

import os
import sys

import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import io_w90
import wannierise as wz

HERE = os.path.dirname(os.path.abspath(__file__))
datadir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "testdata/gaas")
seed = sys.argv[2] if len(sys.argv) > 2 else "gaas"

torch.manual_seed(0)
case = io_w90.load_case(f"{datadir}/{seed}.win", f"{datadir}/{seed}.mmn",
                        f"{datadir}/{seed}.amn", seed)
M0 = torch.from_numpy(case.M)
A = torch.from_numpy(case.A).to(torch.complex128)
bvec = torch.from_numpy(case.bvec)
wb = torch.from_numpy(case.wb)
kb = torch.from_numpy(case.kb_idx)
nk = M0.shape[0]
nw = case.num_wann

U = wz.lowdin_orthonormalize(A)


def omega_total(U):
    M = wz._rotate_M(M0, U, kb)
    _, _, om = wz.spread_decomposition(M, bvec, wb, nk)
    return om["total"]


# analytic gradient at U
M = wz._rotate_M(M0, U, kb)
centres, Mnn, lnim = wz._centres(M, bvec, wb, nk)
G = wz._gradient(M, bvec, wb, centres, Mnn, lnim)

# random anti-Hermitian direction
Xr = torch.randn(nk, nw, nw, dtype=torch.complex128)
X = (Xr - Xr.conj().transpose(-1, -2)) / 2

analytic = sum(
    torch.trace(G[k].conj().transpose(-1, -2) @ X[k]).real.item() for k in range(nk)
)

print(f"case {seed}: nk={nk} nw={nw}   Omega_init={omega_total(U):.8f}")
for eps in (1e-4, 1e-5, 1e-6):
    Up = U @ wz._expm_antihermitian(eps * X)
    Um = U @ wz._expm_antihermitian(-eps * X)
    fd = (omega_total(Up) - omega_total(Um)) / (2 * eps)
    print(f"  eps={eps:.0e}  FD={fd:+.8e}  analytic={analytic:+.8e}  "
          f"ratio={fd / analytic:+.4f}  (expect +1.0000)")
