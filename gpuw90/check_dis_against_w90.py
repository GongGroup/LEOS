"""
check_dis_against_w90.py -- Correctness gate for the disentanglement stage.

Runs disentangle() + (optionally) wannierise() on an entangled gold case and
compares against a serial wannier90.x reference .wout produced from the SAME
.mmn/.amn/.eig/.win inputs.

Acceptance gates (protocol-independent quantities):
  1. Final Omega_I after disentanglement vs the reference 'Final Omega_I'
     (gauge-invariant output of the subspace selection), |d| <= --omega-tol.
  2. Post-disentanglement 'Initial State' centres vs the reference Initial State
     (the Loewdin-projected gauge inside the extracted subspace -- fully
     deterministic, no spread minimisation involved), <= --tol-centre (1e-4 Ang),
     WFs matched by nearest centre (degenerate swaps allowed), minimum-image.

Informational only (NOT gated): the Final State comparison.  The site reference
runs use wannier90's default wannierise protocol -- CG, num_iter=100, NO
convergence check (conv_window=-1) -- and none of the gold runs is converged at
iteration 100 (deltas still 1e-2..1e-1 Ang^2/iter).  Comparing our fully
converged SD minimum against a truncated CG trajectory point is a protocol
difference, not a correctness signal; we print the deviation table for the
record.  Pass --wannierise to include this stage (it dominates the runtime for
the 288-WF case).

Usage:
    python check_dis_against_w90.py testdata/vircef  --seed wannier90.1
    python check_dis_against_w90.py testdata/li2mno3 --seed wannier90.1 --wannierise
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import torch

import io_w90
import wannierise as wz
import disentangle as dz
from check_against_w90 import match_wfs, min_image_delta

TOL_CENTRE = 1e-4   # Ang, project acceptance criterion
TOL_OMEGA_I = 1e-4  # Ang^2, |Omega_I - ref| gate


def compare_states(centres, spreads, ref_c, ref_s, real_lattice, label):
    perm = match_wfs(centres, ref_c, real_lattice)
    nw = len(ref_c)
    max_cdev = 0.0
    max_sdev = 0.0
    for j in range(nw):
        i = perm[j]
        max_cdev = max(max_cdev, float(np.linalg.norm(
            min_image_delta(centres[i], ref_c[j], real_lattice))))
        max_sdev = max(max_sdev, abs(float(spreads[i]) - float(ref_s[j])))
    print(f"{label}: max centre dev = {max_cdev:.3e} Ang, "
          f"max spread dev = {max_sdev:.3e} Ang^2")
    return max_cdev, max_sdev


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("datadir", nargs="?", default="testdata/vircef")
    ap.add_argument("--seed", default="wannier90.1")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--dis-num-iter", type=int, default=None,
                    help="override dis_num_iter (default: .win value or 200)")
    ap.add_argument("--wannierise", action="store_true",
                    help="also run wannierise to convergence (informational "
                         "Final State comparison; slow for 288 WFs on CPU)")
    ap.add_argument("--num-iter", type=int, default=3000,
                    help="wannierise max iterations")
    ap.add_argument("--step", type=float, default=1.0)
    ap.add_argument("--tol-centre", type=float, default=TOL_CENTRE)
    ap.add_argument("--omega-tol", type=float, default=TOL_OMEGA_I)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    d, s = args.datadir, args.seed
    print("=" * 78)
    print(f"gpuw90 disentangle gold check :: {d}/{s}   device={args.device}")
    print("=" * 78)

    t0 = time.perf_counter()
    case = io_w90.load_case(f"{d}/{s}.win", f"{d}/{s}.mmn", f"{d}/{s}.amn", s,
                            eig_path=f"{d}/{s}.eig", allow_entangled=True)
    t_load = time.perf_counter() - t0
    ref = io_w90.parse_wout_final_state(f"{d}/{s}_ref.wout")
    dis_ref = io_w90.parse_wout_dis_info(f"{d}/{s}_ref.wout")

    ki = case.kmesh_info
    print(f"num_bands={case.num_bands} num_wann={case.num_wann} "
          f"num_kpts={case.num_kpts} nntot={case.nntot}   (load {t_load:.1f} s)")
    print(f"kmesh: shells={[f'{r:.6f}' for r in ki['shell_radii']]}")
    print(f"       w_shell={[f'{w:.6f}' for w in ki['w_shell']]}  "
          f"B1_residual={ki['b1_residual']:.2e}")
    w = case.win
    print(f"windows: outer=[{w.get('dis_win_min')}, {w.get('dis_win_max')}] "
          f"frozen=[{w.get('dis_froz_min')}, {w.get('dis_froz_max')}] "
          f"(None -> wannier90 default)")

    # ---- disentangle ---------------------------------------------------------
    t0 = time.perf_counter()
    res = dz.disentangle_case(case, device=args.device,
                              num_iter=args.dis_num_iter, verbose=args.verbose)
    t_dis = time.perf_counter() - t0
    print(f"\ndisentangle: {'converged' if res.converged else 'NOT converged'} "
          f"after {res.n_iter} iters, {t_dis:.1f} s "
          f"({t_dis / max(res.n_iter, 1):.2f} s/iter)")
    print(f"ndimwin: min={res.ndimwin.min()} max={res.ndimwin.max()}   "
          f"ndimfroz: min={res.ndimfroz.min()} max={res.ndimfroz.max()}")

    ref_conv = "converged" if dis_ref["converged"] else "NOT converged"
    print(f"Omega_I: mine={res.omega_I:.8f}  ref={dis_ref['omega_I']:.8f} "
          f"(ref {ref_conv}, {dis_ref['n_iter']} iters)")
    d_omega = abs(res.omega_I - dis_ref["omega_I"])
    print(f"|d Omega_I| = {d_omega:.3e} Ang^2   (gate {args.omega_tol:.0e})")

    # ---- Initial State (post-dis Loewdin gauge; protocol-independent) --------
    bvec = torch.from_numpy(case.bvec)
    wb = torch.from_numpy(case.wb)
    kb_idx = torch.from_numpy(case.kb_idx)
    ci, si, oi = wz.initial_state(res.M0_reduced, res.U0, bvec, wb, kb_idx,
                                  device=args.device)
    print(f"\nInitial State (post-dis): Omega_I={oi['I']:.8f} "
          f"Omega_OD={oi['OD']:.8f} Omega_D={oi['D']:.8f} "
          f"total={oi['total']:.8f}")
    print(f"  sanity |Omega_I(initial) - Omega_I(dis)| = "
          f"{abs(oi['I'] - res.omega_I):.2e} (gauge-invariant, must be ~0)")
    max_cdev = max_sdev = None
    if "initial" in ref:
        max_cdev, max_sdev = compare_states(
            ci, si, ref["initial"]["centres"], ref["initial"]["spreads"],
            case.real_lattice, "Initial State vs ref")
    else:
        print("  (reference .wout has no Initial State block -- gate skipped)")

    # ---- optional wannierise + Final State (informational) -------------------
    if args.wannierise:
        print(f"\nrunning wannierise (SD, step={args.step}, "
              f"max_iter={args.num_iter}) ...")
        t0 = time.perf_counter()
        wres = wz.wannierise(res.M0_reduced, res.U0, bvec, wb, kb_idx,
                             num_iter=args.num_iter, step=args.step,
                             device=args.device, verbose=args.verbose)
        t_wann = time.perf_counter() - t0
        print(f"wannierise: {'converged' if wres.converged else 'NOT converged'} "
              f"after {wres.n_iter} iters, {t_wann:.1f} s "
              f"({t_wann / max(wres.n_iter, 1):.3f} s/iter)")
        om, ro = wres.omega, ref["omega"]
        print(f"Omega   I: mine={om['I']:.8f} ref={ro['I']:.8f}")
        print(f"Omega   D: mine={om['D']:.8f} ref={ro['D']:.8f}")
        print(f"Omega  OD: mine={om['OD']:.8f} ref={ro['OD']:.8f}")
        print(f"Omega TOT: mine={om['total']:.8f} ref={ro['total']:.8f}")
        compare_states(wres.centres, wres.spreads, ref["centres"],
                       ref["spreads"], case.real_lattice,
                       "Final State vs ref [informational: ref = CG truncated "
                       "at 100 iters, unconverged]")

    # ---- verdict --------------------------------------------------------------
    print("\n" + "=" * 78)
    ok = d_omega <= args.omega_tol
    print(f"GATE 1  |d Omega_I| = {d_omega:.3e} <= {args.omega_tol:.0e}  "
          f"{'PASS' if ok else 'FAIL'}")
    if max_cdev is not None:
        ok_c = max_cdev <= args.tol_centre
        print(f"GATE 2  initial-state max centre dev = {max_cdev:.3e} <= "
              f"{args.tol_centre:.0e}  {'PASS' if ok_c else 'FAIL'}")
        ok = ok and ok_c
    print(f"VERDICT: {'PASS' if ok else 'FAIL'}")
    print("=" * 78)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
