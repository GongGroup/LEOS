"""
bench_gpu.py -- dis + wannierise wall-clock benchmark (CPU baseline vs CUDA).

Runs the full Phase-2 pipeline on one spin channel with FIXED iteration counts
(no convergence early-exit) so s/iter is comparable across devices, then prints
a compact table.  Correctness is NOT judged here -- run check_dis_against_w90.py
for the gold gates first.

Local CPU baseline (recorded 2026-07-11, Apple M-series, torch 2.8 CPU):
    python bench_gpu.py testdata/li2mno3 --seed wannier90.1

sai-pod V100 (16V100 partition needs >=4 GPUs/job; we use 1 of 4 -- see
run_bench_v100.sbatch next to this file):
    python bench_gpu.py <datadir> --seed wannier90.1 --device cuda
"""

from __future__ import annotations

import argparse
import time

import torch

import io_w90
import wannierise as wz
import disentangle as dz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("datadir")
    ap.add_argument("--seed", default="wannier90.1")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--dis-iters", type=int, default=50,
                    help="fixed disentanglement iterations (conv check disabled)")
    ap.add_argument("--wann-iters", type=int, default=100,
                    help="fixed wannierise iterations (conv check disabled)")
    ap.add_argument("--check-every", type=int, default=10,
                    help="host poll stride inside the hot loops")
    args = ap.parse_args()

    d, s, dev = args.datadir, args.seed, args.device
    if dev.startswith("cuda") and not torch.cuda.is_available():
        raise SystemExit("CUDA requested but not available (check sm_70 wheel: "
                         "V100 needs a torch<=2.6/cu124-era build; the site "
                         "mace0.3.15 env ships cu128 without sm_70)")

    t0 = time.perf_counter()
    case = io_w90.load_case(f"{d}/{s}.win", f"{d}/{s}.mmn", f"{d}/{s}.amn", s,
                            eig_path=f"{d}/{s}.eig", allow_entangled=True)
    t_load = time.perf_counter() - t0
    print(f"case {d}/{s}: nb={case.num_bands} nw={case.num_wann} "
          f"nk={case.num_kpts} nnb={case.nntot}  device={dev}  load {t_load:.1f} s")

    # conv_tol=0 -> the |delta| < tol window can never trigger: fixed-length runs
    t0 = time.perf_counter()
    res = dz.disentangle_case(case, device=dev, num_iter=args.dis_iters,
                              conv_tol=0.0, check_every=args.check_every)
    if dev.startswith("cuda"):
        torch.cuda.synchronize()
    t_dis = time.perf_counter() - t0

    bvec = torch.from_numpy(case.bvec)
    wb = torch.from_numpy(case.wb)
    kb_idx = torch.from_numpy(case.kb_idx)
    t0 = time.perf_counter()
    wres = wz.wannierise(res.M0_reduced, res.U0, bvec, wb, kb_idx,
                         num_iter=args.wann_iters, conv_tol=0.0,
                         device=dev, check_every=args.check_every)
    if dev.startswith("cuda"):
        torch.cuda.synchronize()
    t_wann = time.perf_counter() - t0

    print(f"\n{'stage':<14}{'iters':>7}{'total s':>10}{'s/iter':>10}")
    print(f"{'disentangle':<14}{res.n_iter:>7}{t_dis:>10.2f}"
          f"{t_dis / max(res.n_iter, 1):>10.3f}")
    print(f"{'wannierise':<14}{wres.n_iter:>7}{t_wann:>10.2f}"
          f"{t_wann / max(wres.n_iter, 1):>10.3f}")
    print(f"\nOmega_I after {res.n_iter} dis iters: {res.omega_I:.8f} Ang^2  "
          f"(fixed-iteration run; NOT the converged value)")
    print(f"Omega_total after {wres.n_iter} wann iters: "
          f"{wres.omega['total']:.8f} Ang^2")


if __name__ == "__main__":
    main()
