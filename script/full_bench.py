#!/usr/bin/env python3
"""Full 74-case benchmark: transport OS + two-level reliability R.

Covers all four batches (v1 dataset, v1 oxides, M4 wave1, M4 wave2) against
pyleos/tests/reference_os.json. Emits a TSV (default data/v3/reliability_74.tsv)
with per-case dev and the R_elem / R_site reliability scores, plus a summary
of how R_elem separates passing from failing cases (the indicator-validity
check for the paper's R(%) comparison table).

Usage: python3 script/full_bench.py [--out FILE] [--beta B] [--gamma G]
       [--no-reliability]
"""
from __future__ import annotations

import argparse
import functools
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pyleos.parsers import load_case                      # noqa: E402
from pyleos.transport import assign_case, reliability_margins   # noqa: E402

print = functools.partial(print, flush=True)
REF = json.loads((ROOT / "pyleos/tests/reference_os.json").read_text())


def find_oxide_dir(name: str):
    for d in sorted((ROOT / "data/v1/oxides").glob(f"leos-oxide-{name}-*")):
        if (d / "EOSCAR").is_file():
            return d
    return None


def collect():
    """Yield (batch, name, dir, ref_map)."""
    for name, meta in REF["dataset"].items():
        d = ROOT / "data/v1/eos_dataset" / name
        if d.is_dir() and (d / "EOSCAR").is_file():
            yield "v1_mol", name, d, meta.get("ref") or meta.get("ref_mean") or {}
    for name, meta in REF["oxides"].items():
        d = find_oxide_dir(name)
        if d is not None:
            yield "v1_oxide", name, d, meta.get("ref") or meta.get("ref_mean") or {}
    for batch, sub in (("wave1", "m4_wave1"), ("wave2", "m4_wave2")):
        for name, meta in REF[batch].items():
            if name.startswith("_"):
                continue
            d = ROOT / "data/v3" / sub / name
            if d.is_dir() and (d / "EOSCAR").is_file():
                yield batch, name, d, (meta.get("ref") or meta.get("ref_mean")
                                       or {}) if isinstance(meta, dict) else {}


def element_means(elements, values):
    by = defaultdict(list)
    for el, v in zip(elements, values):
        by[el].append(v)
    return {el: sum(v) / len(v) for el, v in by.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=ROOT / "data/v3/reliability_74.tsv")
    ap.add_argument("--beta", type=float, default=10.0)
    ap.add_argument("--gamma", type=float, default=0.0)
    ap.add_argument("--tier2", default="linear",
                    choices=["linear", "arbiter", "baseline"])
    ap.add_argument("--only-wp", action="store_true",
                    help="restrict to cases with a #WP block (tier-2 bench)")
    ap.add_argument("--no-reliability", action="store_true")
    ap.add_argument("--greedy", action="store_true",
                    help="score the in-situ greedy OS column of EOSCAR "
                         "instead of the transport assignment (same cases, "
                         "same dev criterion; R columns are meaningless "
                         "here and emitted as '-')")
    args = ap.parse_args()

    rows = []
    t_start = time.time()
    cases = list(collect())
    for i, (batch, name, d, ref_map) in enumerate(cases):
        t0 = time.time()
        print(f"[{i+1:2d}/{len(cases)}] {name} ...", end=" ")
        try:
            case = load_case(d)
            if args.only_wp and not (case.eoscar and case.eoscar.wp_detail):
                print("skip (no #WP)")
                continue
            if args.greedy:
                res = None
            else:
                res = assign_case(case, beta=args.beta, gamma=args.gamma,
                                  tier2=args.tier2)
                if not args.no_reliability:
                    res = reliability_margins(case, beta=args.beta,
                                              gamma=args.gamma, res=res,
                                              tier2=args.tier2)
        except Exception as exc:
            print(f"ERROR {str(exc)[:50]}")
            rows.append((batch, name, "ERROR", "", "", "", "", ""))
            continue
        os_vals = case.eoscar.greedy_os if args.greedy else res.os_int
        means = element_means(case.structure.elements, os_vals)
        dev = max((abs(means[el] - t) for el, t in ref_map.items() if el in means),
                  default=None)
        status = ("run" if dev is None else
                  "OK" if dev <= 0.25 else "~" if dev <= 1.0 else "FAIL")
        has_r = res is not None and res.atom_margin
        rows.append((batch, name, status,
                     f"{dev:.2f}" if dev is not None else "-",
                     f"{res.case_conf_elem:.1f}" if has_r else "-",
                     f"{res.case_conf:.1f}" if has_r else "-",
                     len(res.stray) if res is not None else "-",
                     res.ghost if res is not None else "-"))
        print(f"{status} dev={rows[-1][3]} R_elem={rows[-1][4]} "
              f"({time.time()-t0:.1f}s)")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        fh.write("batch\tcase\tstatus\tdev\tR_elem\tR_site\tn_stray\tn_ghost\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    # Indicator validity: does low R_elem flag the failing cases?
    def rvals(pred):
        return [float(r[4]) for r in rows if r[4] not in ("-", "") and pred(r[2])]
    ok_r, bad_r = rvals(lambda s: s == "OK"), rvals(lambda s: s in ("~", "FAIL"))
    print(f"\n{len(rows)} cases in {time.time()-t_start:.0f}s -> {args.out}")
    if ok_r and bad_r:
        med = lambda v: sorted(v)[len(v) // 2]
        print(f"R_elem median: OK cases {med(ok_r):.0f} "
              f"vs deviating cases {med(bad_r):.0f}")


if __name__ == "__main__":
    main()
