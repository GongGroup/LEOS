#!/usr/bin/env python3
"""Collapse a frames_to_dft output tree into an oxidation-state time series.

Scans <root>/frame_*/EOSCAR, runs the transport assignment (tier-2 arbiter
default) plus reliability on each frame, and writes one TSV row per frame:
per-element mean OS (transport and greedy), R_elem, quality flags. This is
the M5-B production deliverable: metal vs oxygen ledger along the MACE
trajectory, with the greedy column doubling as the large-system drift
comparison.

Usage: python3 tools/timeseries.py m5b_dft/li2mno3/x0.25 --dt-ps 0.5 \
           [--out timeseries.tsv]
"""
from __future__ import annotations

import argparse
import functools
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pyleos.parsers import load_case                              # noqa: E402
from pyleos.transport import assign_case, reliability_margins     # noqa: E402

print = functools.partial(print, flush=True)


def element_means(elements, values):
    by = defaultdict(list)
    for el, v in zip(elements, values):
        if v is not None:
            by[el].append(v)
    return {el: sum(v) / len(v) for el, v in by.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root", type=Path)
    ap.add_argument("--dt-ps", type=float, default=0.5,
                    help="time between dumped frames in ps")
    ap.add_argument("--beta", type=float, default=10.0)
    ap.add_argument("--tier2", default="arbiter")
    ap.add_argument("--reliability", action="store_true")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    frames = sorted(args.root.glob("frame_*"))
    done = [d for d in frames if (d / "EOSCAR").is_file()]
    print(f"{len(done)}/{len(frames)} frames have EOSCAR")
    if not done:
        return
    rows, elements_seen = [], []
    for d in done:
        k = int(d.name.split("_")[1])
        case = load_case(d)
        res = assign_case(case, beta=args.beta, tier2=args.tier2)
        if args.reliability:
            res = reliability_margins(case, beta=args.beta, res=res,
                                      tier2=args.tier2)
        els = case.structure.elements
        tmean = element_means(els, res.os_int)
        gmean = element_means(els, case.eoscar.greedy_os) if case.eoscar else {}
        for el in tmean:
            if el not in elements_seen:
                elements_seen.append(el)
        rows.append((k, k * args.dt_ps, tmean, gmean,
                     res.case_conf_elem if args.reliability else None,
                     len(res.stray), res.ghost))
        print(f"frame {k:3d}  " + "  ".join(
            f"{el}:{tmean[el]:+.3f}" for el in sorted(tmean)))

    out = args.out or (args.root / "timeseries.tsv")
    with out.open("w") as fh:
        head = ["frame", "t_ps"]
        for el in elements_seen:
            head += [f"{el}_flow", f"{el}_greedy"]
        head += ["R_elem", "n_stray", "n_ghost"]
        fh.write("\t".join(head) + "\n")
        for k, t, tmean, gmean, r, ns, ng in rows:
            row = [str(k), f"{t:.2f}"]
            for el in elements_seen:
                row.append(f"{tmean.get(el, float('nan')):.4f}")
                row.append(f"{gmean.get(el, float('nan')):.4f}"
                           if el in gmean else "-")
            row += [f"{r:.1f}" if r is not None else "-", str(ns), str(ng)]
            fh.write("\t".join(row) + "\n")
    print(f"wrote {out}")

    # stability summary: max excursion of each element's transport mean
    print("\nstability (transport): max |dev from trajectory mean|")
    for el in elements_seen:
        vs = [r[2][el] for r in rows if el in r[2]]
        mu = sum(vs) / len(vs)
        print(f"  {el}: mean {mu:+.3f}, max excursion "
              f"{max(abs(v - mu) for v in vs):.3f}")


if __name__ == "__main__":
    main()
