#!/usr/bin/env python3
"""Build the M5-B three-system Li-rich cathode structure set.

Systems (paper's applicability-boundary scan):
  li2mno3  - ordered 3d baseline. C2/m, Strobel & Lambert-Andron 1988
             (a=4.937, b=8.532, c=5.030, beta=109.46), 2x1x2 supercell.
  li2ruo3  - 4d touchstone: same framework, Ru lattice parameters
             (Miura 2007-ish; the starting geometry is relaxed downstream).
  nmc      - Li-rich disordered Li1.2Ni0.13Mn0.54Co0.13O2 on an R-3m
             layered framework; TM-layer occupation randomised with a
             poor-man's SQS (fixed seed, best-of-N nearest-neighbour
             correlation match to the random limit).

Delithiation points x = 0 / 0.25 / 0.50 (fraction of Li removed), Li picked
to keep the remaining Li sublattice maximally spread (best-of-N seeded).

Run on sai-pod inside the mace env (needs ase):
  source /opt/envs/mace0.3.15.env && python3 build_m5b.py --out m5b_structures
"""
from __future__ import annotations

import argparse
import itertools
import random
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import make_supercell
from ase.io import write
from ase.spacegroup import crystal

DELITH = (0.0, 0.25, 0.50)
SEED = 42
NTRIAL = 200


def li2mno3(tm="Mn"):
    a, b, c, beta = 4.937, 8.532, 5.030, 109.46
    if tm == "Ru":
        a, b, c, beta = 4.92, 8.77, 5.11, 108.7
    prim = crystal(
        symbols=["Li", "Li", "Li", tm, "O", "O"],
        basis=[(0.0, 0.5, 0.0),        # Li 2b
               (0.0, 0.0, 0.5),        # Li 2c
               (0.0, 0.6606, 0.5),     # Li 4h
               (0.0, 0.16708, 0.0),    # TM 4g
               (0.2189, 0.0, 0.2273),  # O 4i
               (0.2540, 0.32119, 0.2233)],  # O 8j
        spacegroup=12, cellpar=[a, b, c, 90.0, beta, 90.0])
    sup = make_supercell(prim, np.diag([2, 1, 2]))
    return sup


def nmc_lirich():
    """Layered R-3m framework, 4x3x1 of the rhombohedral-hex cell in its
    hexagonal setting: alternating Li layer / O / TM layer / O stacking."""
    a, c = 2.86, 14.23
    prim = crystal(symbols=["Li", "Co", "O"],
                   basis=[(0, 0, 0), (0, 0, 0.5), (0, 0, 0.2599)],
                   spacegroup=166, cellpar=[a, a, c, 90, 90, 120])
    sup = make_supercell(prim, np.diag([4, 2, 1]))   # 24 Li + 24 TM + 48 O
    # TM-layer target occupation for Li1.2(Ni.13Mn.54Co.13)O2 on 24 sites:
    # 5 excess Li + 3 Ni + 13 Mn + 3 Co  -> Li1.208 TM0.792 O2
    fill = ["Li"] * 5 + ["Ni"] * 3 + ["Mn"] * 13 + ["Co"] * 3
    tm_idx = [i for i, s in enumerate(sup.get_chemical_symbols()) if s == "Co"]
    assert len(tm_idx) == len(fill)
    rng = random.Random(SEED)
    d = sup.get_all_distances(mic=True)
    nn = {i: sorted(tm_idx, key=lambda j: d[i][j])[1:7] for i in tm_idx}  # 6 in-plane NN
    frac = {s: fill.count(s) / len(fill) for s in set(fill)}

    def score(assign):
        # deviation of same/cross-pair counts from the random limit
        pair = {}
        for i in tm_idx:
            for j in nn[i]:
                if j > i:
                    key = tuple(sorted((assign[i], assign[j])))
                    pair[key] = pair.get(key, 0) + 1
        tot = sum(pair.values())
        dev = 0.0
        for x, y in itertools.combinations_with_replacement(sorted(frac), 2):
            expect = (frac[x] ** 2 if x == y else 2 * frac[x] * frac[y])
            dev += (pair.get((x, y), 0) / tot - expect) ** 2
        return dev

    best, best_dev = None, 1e9
    for _ in range(NTRIAL):
        perm = fill[:]
        rng.shuffle(perm)
        assign = dict(zip(tm_idx, perm))
        s = score(assign)
        if s < best_dev:
            best, best_dev = assign, s
    syms = sup.get_chemical_symbols()
    for i, s in best.items():
        syms[i] = s
    sup.set_chemical_symbols(syms)
    print(f"  nmc SQS pair-correlation deviation: {best_dev:.5f}")
    return sup


def delithiate(atoms: Atoms, x: float) -> Atoms:
    """Remove round(x * nLi) Li keeping the survivors maximally spread."""
    li = [i for i, s in enumerate(atoms.get_chemical_symbols()) if s == "Li"]
    n_rm = round(x * len(li))
    if n_rm == 0:
        return atoms.copy()
    rng = random.Random(SEED + int(100 * x))
    d = atoms.get_all_distances(mic=True)
    best, best_min = None, -1.0
    for _ in range(NTRIAL):
        keep = rng.sample(li, len(li) - n_rm)
        m = min(d[i][j] for i in keep for j in keep if j > i)
        if m > best_min:
            best, best_min = set(keep), m
    out = atoms.copy()
    del out[[i for i in li if i not in best]]
    print(f"  x={x}: removed {n_rm} Li, min Li-Li of survivors {best_min:.2f} A")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path("m5b_structures"))
    args = ap.parse_args()
    systems = {
        "li2mno3": li2mno3("Mn"),
        "li2ruo3": li2mno3("Ru"),
        "nmc": nmc_lirich(),
    }
    for name, atoms in systems.items():
        print(f"[{name}] {atoms.get_chemical_formula()} "
              f"({len(atoms)} atoms, cell {atoms.cell.lengths().round(2)})")
        for x in DELITH:
            d = args.out / name / f"x{x:.2f}"
            d.mkdir(parents=True, exist_ok=True)
            sub = delithiate(atoms, x)
            write(d / "POSCAR", sub, format="vasp", direct=True, sort=True)
    print("BUILD_M5B_DONE")


if __name__ == "__main__":
    main()
