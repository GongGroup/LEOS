"""
check_against_w90.py -- Correctness gate against a wannier90 reference .wout.

Runs the batched PyTorch wannierisation on a gold-standard case and compares Wannier
centres / spreads against wannier90's 'Final State'.  Acceptance (set by the project):
physical equivalence -- centre deviation <= 1e-4 Ang, with near-degenerate WF swaps
allowed (WFs are matched by nearest centre, not by index).

Usage:
    python check_against_w90.py [testdata_dir] [--device cpu|cuda] [--num-iter N]
                                [--step S] [--verbose]
"""

from __future__ import annotations

import argparse

import numpy as np
import torch

import io_w90
import wannierise as wz

TOL_CENTRE = 1e-4  # Ang, project acceptance criterion


def min_image_delta(c_a, c_b, real_lattice):
    """Minimum-image Cartesian distance between two centres (periodic cell)."""
    inv = np.linalg.inv(real_lattice)
    d = c_a - c_b
    f = d @ inv
    f -= np.round(f)                 # wrap to [-0.5, 0.5)
    return f @ real_lattice


def match_wfs(centres_mine, centres_ref, real_lattice):
    """Greedy nearest-centre assignment (allows degenerate WF swaps).

    Returns perm such that centres_mine[perm[j]] matches centres_ref[j].
    """
    nw = len(centres_ref)
    used = [False] * nw
    perm = [-1] * nw
    for j in range(nw):
        best, best_d = -1, 1e30
        for i in range(nw):
            if used[i]:
                continue
            d = np.linalg.norm(min_image_delta(centres_mine[i], centres_ref[j],
                                               real_lattice))
            if d < best_d:
                best_d, best = d, i
        perm[j] = best
        used[best] = True
    return perm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("datadir", nargs="?", default="testdata/gaas")
    ap.add_argument("--seed", default="gaas")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--num-iter", type=int, default=2000)
    ap.add_argument("--step", type=float, default=1.0)
    ap.add_argument("--conv-tol", type=float, default=1e-10)
    ap.add_argument("--conv-window", type=int, default=3)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    d = args.datadir
    s = args.seed
    case = io_w90.load_case(f"{d}/{s}.win", f"{d}/{s}.mmn", f"{d}/{s}.amn", s)
    ref = io_w90.parse_wout_final_state(f"{d}/{s}_ref.wout")

    print("=" * 78)
    print(f"gpuw90 correctness check :: {d}/{s}   device={args.device}")
    print("=" * 78)
    print(f"num_bands={case.num_bands} num_wann={case.num_wann} "
          f"num_kpts={case.num_kpts} nntot={case.nntot}")
    ki = case.kmesh_info
    print(f"kmesh: shells={[f'{r:.6f}' for r in ki['shell_radii']]}  "
          f"w_shell={[f'{w:.6f}' for w in ki['w_shell']]}  "
          f"B1_residual={ki['b1_residual']:.2e}")
    assert case.num_bands == case.num_wann, \
        "this stage handles isolated bands only (nb == nw); disentangle is next stage"

    # torch inputs
    M0 = torch.from_numpy(case.M)
    A = torch.from_numpy(case.A)
    bvec = torch.from_numpy(case.bvec)
    wb = torch.from_numpy(case.wb)
    kb_idx = torch.from_numpy(case.kb_idx)

    # initial gauge: Loewdin-orthonormalised projections
    U0 = wz.lowdin_orthonormalize(A.to(torch.complex128))

    # --- Initial State sanity check (validates io + weights + centre formula) ---
    if "initial" in ref:
        ci, si, oi = wz.initial_state(M0, U0, bvec, wb, kb_idx, device=args.device)
        perm = match_wfs(ci, ref["initial"]["centres"], case.real_lattice)
        maxdev_i = max(
            np.linalg.norm(min_image_delta(ci[perm[j]], ref["initial"]["centres"][j],
                                           case.real_lattice))
            for j in range(case.num_wann)
        )
        print(f"\nInitial State: Omega_total mine={oi['total']:.8f} "
              f"ref={ref['initial']['spreads'].sum():.8f}  "
              f"max centre dev={maxdev_i:.2e} Ang")

    # --- Run wannierise ---
    print(f"\nrunning wannierise (SD, step={args.step}, max_iter={args.num_iter}) ...")
    res = wz.wannierise(
        M0, U0, bvec, wb, kb_idx,
        num_iter=args.num_iter, step=args.step,
        conv_tol=args.conv_tol, conv_window=args.conv_window,
        device=args.device, verbose=args.verbose,
    )
    print(f"converged={res.converged} after {res.n_iter} iters   "
          f"Omega_total={res.omega['total']:.9f}")

    # --- Compare Final State ---
    perm = match_wfs(res.centres, ref["centres"], case.real_lattice)
    print("\nFinal State: per-WF comparison (matched by nearest centre)")
    print("-" * 78)
    print(f"{'WF':>3} {'centre (mine, Ang)':>34} {'|d_centre|':>12} "
          f"{'spread_mine':>12} {'d_spread':>11}")
    print("-" * 78)
    max_cdev = 0.0
    max_sdev = 0.0
    for j in range(case.num_wann):
        i = perm[j]
        cm = res.centres[i]
        cr = ref["centres"][j]
        dcen = np.linalg.norm(min_image_delta(cm, cr, case.real_lattice))
        dspr = abs(res.spreads[i] - ref["spreads"][j])
        max_cdev = max(max_cdev, dcen)
        max_sdev = max(max_sdev, dspr)
        print(f"{j+1:>3} ({cm[0]:9.6f},{cm[1]:9.6f},{cm[2]:9.6f}) "
              f"{dcen:12.2e} {res.spreads[i]:12.6f} {dspr:11.2e}")
    print("-" * 78)

    om = res.omega
    ro = ref["omega"]
    print(f"Omega   I: mine={om['I']:.9f} ref={ro['I']:.9f} "
          f"d={abs(om['I']-ro['I']):.2e}")
    print(f"Omega   D: mine={om['D']:.9f} ref={ro['D']:.9f} "
          f"d={abs(om['D']-ro['D']):.2e}")
    print(f"Omega  OD: mine={om['OD']:.9f} ref={ro['OD']:.9f} "
          f"d={abs(om['OD']-ro['OD']):.2e}")
    print(f"Omega TOT: mine={om['total']:.9f} ref={ro['total']:.9f} "
          f"d={abs(om['total']-ro['total']):.2e}")

    print("\n" + "=" * 78)
    print(f"MAX CENTRE DEVIATION = {max_cdev:.3e} Ang   (tolerance {TOL_CENTRE:.0e})")
    print(f"MAX SPREAD DEVIATION = {max_sdev:.3e} Ang^2")
    verdict = "PASS" if max_cdev <= TOL_CENTRE else "FAIL"
    print(f"VERDICT: {verdict}")
    print("=" * 78)
    return 0 if max_cdev <= TOL_CENTRE else 1


if __name__ == "__main__":
    raise SystemExit(main())
