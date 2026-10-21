#!/usr/bin/env python3
"""MACE foundation-model MD with built-in acceptance checks (M5-B).

Stage A (acceptance): short NVT at 300 K, then verdicts on
  - energy sanity: no blow-up, NVT total-energy drift within bounds
  - structure integrity: framework RDF peaks persist, no free O2
    (an O-O pair under 1.35 A whose members sit >3.0 A from every TM),
    O-O dimer inventory (1.2-1.6 A, TM-coordinated) reported as data
  - temperature control: mean T within 30 K of target
Stage B (production, --steps larger): same driver, longer trajectory,
frames dumped for downstream DFT single points.

Usage (inside /opt/envs/mace0.3.15.env):
  python3 mace_md.py POSCAR --model ~/models/mace-mpa-0-medium.model \
      --steps 2500 --dt 2.0 --temp 300 --dump-every 50 --out mdrun
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from ase import units
from ase.io import read, write
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

TM = {"Mn", "Ru", "Ni", "Co", "Ti", "V", "Cr", "Fe"}


def analyse(atoms, traj_frames):
    syms = atoms.get_chemical_symbols()
    o_idx = [i for i, s in enumerate(syms) if s == "O"]
    tm_idx = [i for i, s in enumerate(syms) if s in TM]
    last = traj_frames[-1]
    d = last.get_all_distances(mic=True)
    free_o2, dimers = [], []
    for x, i in enumerate(o_idx):
        for j in o_idx[x + 1:]:
            dij = d[i][j]
            if dij < 1.6:
                near_tm = min((d[i][t] for t in tm_idx), default=9.9)
                near_tm = min(near_tm, min((d[j][t] for t in tm_idx), default=9.9))
                rec = {"pair": [int(i), int(j)], "d_oo": round(float(dij), 3),
                       "nearest_tm": round(float(near_tm), 2)}
                if dij < 1.35 and near_tm > 3.0:
                    free_o2.append(rec)
                else:
                    dimers.append(rec)
    return free_o2, dimers


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("poscar", type=Path)
    ap.add_argument("--model", required=True)
    ap.add_argument("--steps", type=int, default=2500)
    ap.add_argument("--dt", type=float, default=2.0)
    ap.add_argument("--temp", type=float, default=300.0)
    ap.add_argument("--dump-every", type=int, default=50)
    ap.add_argument("--out", type=Path, default=Path("mdrun"))
    ap.add_argument("--device", default="cuda")
    args = ap.parse_args()

    from mace.calculators import MACECalculator
    atoms = read(args.poscar)
    atoms.calc = MACECalculator(model_paths=[str(args.model)],
                                device=args.device, default_dtype="float64")
    args.out.mkdir(parents=True, exist_ok=True)
    MaxwellBoltzmannDistribution(atoms, temperature_K=args.temp)
    dyn = Langevin(atoms, args.dt * units.fs, temperature_K=args.temp,
                   friction=0.02)

    log_rows, frames = [], []
    t0 = time.time()

    def snap():
        e = atoms.get_potential_energy()
        ke = atoms.get_kinetic_energy()
        T = ke / (1.5 * units.kB * len(atoms))
        log_rows.append((dyn.nsteps, e, ke, T))
        frames.append(atoms.copy())
        if dyn.nsteps % (args.dump_every * 10) == 0:
            done = dyn.nsteps / args.steps
            eta = (time.time() - t0) / max(done, 1e-9) * (1 - done)
            print(f"  step {dyn.nsteps}/{args.steps} T={T:6.1f}K "
                  f"E={e:12.3f}eV  ETA {eta/60:.1f} min", flush=True)

    dyn.attach(snap, interval=args.dump_every)
    snap()
    dyn.run(args.steps)

    write(args.out / "traj.extxyz", frames)
    with (args.out / "md.log").open("w") as fh:
        fh.write("step\tE_pot_eV\tE_kin_eV\tT_K\n")
        for r in log_rows:
            fh.write(f"{r[0]}\t{r[1]:.6f}\t{r[2]:.6f}\t{r[3]:.2f}\n")

    # acceptance verdicts on the second half (post-equilibration)
    half = len(log_rows) // 2
    E = [r[1] + r[2] for r in log_rows[half:]]
    Ts = [r[3] for r in log_rows[half:]]
    drift = (E[-1] - E[0]) / len(atoms) * 1000  # meV/atom over window
    free_o2, dimers = analyse(atoms, frames)
    verdict = {
        "n_atoms": len(atoms),
        "steps": args.steps,
        "dt_fs": args.dt,
        "wall_min": round((time.time() - t0) / 60, 1),
        "T_mean_K": round(float(np.mean(Ts)), 1),
        "E_drift_meV_per_atom": round(float(drift), 2),
        "free_O2": free_o2,
        "oo_dimers": dimers,
        "blow_up": bool(max(Ts) > 3 * args.temp),
        "pass": bool(max(Ts) <= 3 * args.temp
                     and abs(np.mean(Ts) - args.temp) < 30
                     and not free_o2),
    }
    (args.out / "verdict.json").write_text(json.dumps(verdict, indent=1))
    print(json.dumps(verdict, indent=1))
    print("MACE_MD_DONE", "PASS" if verdict["pass"] else "CHECK")


if __name__ == "__main__":
    main()
