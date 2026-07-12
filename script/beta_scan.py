#!/usr/bin/env python3
"""Beta plateau scan: the 74-case benchmark at several beta values.

Pins the electronegativity-weight plateau to data instead of prose (the
audited "beta 10-15 is flat" sentence was falsified: 15 flips Cr(CO)6,
20 flips tmao).  Every beta is scored with the same element-mean
dev <= 0.25 criterion as full_bench.py; cases without a reference stay
out of the counts (denominator 72).

Output: one row per beta with pass counts and the cases whose status
differs from the beta=10 default.

Usage: python3 script/beta_scan.py [--out FILE]
"""
from __future__ import annotations

import argparse
import functools
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "script"))

from full_bench import collect, element_means          # noqa: E402
from pyleos.parsers import load_case                   # noqa: E402
from pyleos.transport import assign_case               # noqa: E402

print = functools.partial(print, flush=True)

BETAS = [3.0, 5.0, 8.0, 9.0, 10.0, 11.0, 12.0, 15.0, 20.0]
REF_BETA = 10.0


def score(case, ref_map, beta):
    try:
        res = assign_case(case, beta=beta)
    except Exception:
        return "ERROR"
    means = element_means(case.structure.elements, res.os_int)
    dev = max((abs(means[el] - t) for el, t in ref_map.items() if el in means),
              default=None)
    if dev is None:
        return "run"
    return "OK" if dev <= 0.25 else "~" if dev <= 1.0 else "FAIL"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "data/v3/beta_scan_20260712.tsv")
    args = ap.parse_args()

    # key by (batch, name): case names are not unique across batches (tio2)
    cases = [(batch, name, load_case(d), ref)
             for batch, name, d, ref in collect()]
    status = {}
    for i, beta in enumerate(BETAS):
        t0 = time.time()
        for batch, name, case, ref in cases:
            status[(beta, batch, name)] = score(case, ref, beta)
        counts = {s: sum(1 for (b, *_), v in status.items()
                         if b == beta and v == s)
                  for s in ("OK", "~", "FAIL")}
        print(f"[{i+1}/{len(BETAS)}] beta={beta:4.1f}  OK={counts['OK']}"
              f"  ~={counts['~']}  FAIL={counts['FAIL']}"
              f"  ({time.time()-t0:.0f}s)")

    with args.out.open("w") as fh:
        fh.write("beta\tn_OK\tn_mid\tn_FAIL\tchanged_vs_beta10\n")
        for beta in BETAS:
            def n_with(s, b=beta):
                return sum(1 for bt, n, _, _ in cases
                           if status[(b, bt, n)] == s)
            changed = ",".join(
                f"{n}:{status[(REF_BETA, bt, n)]}->{status[(beta, bt, n)]}"
                for bt, n, _, _ in cases
                if status[(beta, bt, n)] != status[(REF_BETA, bt, n)]) or "-"
            fh.write(f"{beta:g}\t{n_with('OK')}\t{n_with('~')}"
                     f"\t{n_with('FAIL')}\t{changed}\n")
    print(f"-> {args.out}")


if __name__ == "__main__":
    main()
