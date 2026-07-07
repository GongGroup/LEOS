#!/usr/bin/env python3
"""pyleos CLI: transport OS assignment for single cases and regression runs.

Usage:
  python3 pyleos/cli.py case <case_dir> [--beta B] [--eneg allen|pauling] [-v]
  python3 pyleos/cli.py regress [--root leos_work] [--beta B] [--eneg SCALE]
"""
from __future__ import annotations

import argparse
import functools
import json
import sys
from collections import defaultdict
from pathlib import Path

print = functools.partial(print, flush=True)   # progress must survive redirection

if __package__ is None or __package__ == "":          # direct-script support
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from pyleos.parsers import load_case                # noqa: E402
    from pyleos.transport import assign_case            # noqa: E402
else:
    from .parsers import load_case
    from .transport import assign_case

REF_PATH = Path(__file__).resolve().parent / "tests" / "reference_os.json"


def element_summary(elements, values):
    by_el = defaultdict(list)
    for el, v in zip(elements, values):
        by_el[el].append(v)
    return {el: (min(v), max(v), sum(v) / len(v)) for el, v in sorted(by_el.items())}


def run_case(case_dir: Path, beta: float, eneg_scale: str, verbose: bool = False,
             drop_stray: bool = False, gamma: float = 0.0):
    case = load_case(case_dir)
    res = assign_case(case, beta=beta, eneg_scale=eneg_scale, drop_stray=drop_stray,
                      gamma=gamma)
    st = case.structure
    print(f"\n=== {case.name}  (ISPIN={case.ispin}, WFs={len(res.wfs)}, "
          f"electrons={res.supply_total}) ===")
    if case.eoscar and case.eoscar.dne is not None:
        print(f"greedy header: NELECT(Wannier)={case.eoscar.nelect_wann}  "
              f"|dNe|={case.eoscar.dne:.4f}")
    if res.ghost:
        print(f"note: {res.ghost} ghost electrons excluded "
              f"(num_wann exceeds occupied manifold; flow capped at NELECT)")
    if res.infeasible:
        print(f"WARNING: {res.infeasible} electrons could not be routed "
              f"(capacity infeasibility)")
    if res.stray:
        worst = max(res.stray, key=lambda w: w.spread)
        print(f"WARNING: {len(res.stray)} stray/unlocalised WFs "
              f"(worst spread {worst.spread:.1f} A^2) - Wannierisation suspect")
    if verbose:
        print(f"{'ion':>4} {'el':>3} {'ne':>4} {'OS_flow':>8} {'OS_greedy':>10}")
        for a, el in enumerate(st.elements):
            g = case.eoscar.greedy_os[a] if case.eoscar else None
            print(f"{a+1:>4} {el:>3} {res.ne[a]:>4} {res.os_int[a]:>+8} "
                  f"{g if g is None else format(g, '+d'):>10}")
    flow_sum = element_summary(st.elements, res.os_int)
    greedy_sum = (element_summary(st.elements, case.eoscar.greedy_os)
                  if case.eoscar else {})
    print(f"{'el':>4} {'flow min..max (mean)':>24} {'greedy min..max (mean)':>26}")
    for el, (lo, hi, mean) in flow_sum.items():
        gtxt = ""
        if el in greedy_sum:
            glo, ghi, gmean = greedy_sum[el]
            gtxt = f"{glo:+d}..{ghi:+d} ({gmean:+.2f})"
        print(f"{el:>4} {f'{lo:+d}..{hi:+d} ({mean:+.2f})':>24} {gtxt:>26}")
    return case, res


def find_oxide_dir(root: Path, name: str):
    for d in sorted(root.glob(f"leos-oxide-{name}-*")):
        if (d / "EOSCAR").is_file() and (
                (d / "wannier90.wout").is_file() or (d / "wannier90.1.wout").is_file()):
            return d
    return None


def regress(root: Path, beta: float, eneg_scale: str):
    ref = json.loads(REF_PATH.read_text())
    rows = []
    targets = []
    ds_root = root / "leos_eos_dataset"
    for name, meta in ref["dataset"].items():
        d = ds_root / name
        if d.is_dir() and (d / "EOSCAR").is_file():
            targets.append((name, d, meta))
    for name, meta in ref["oxides"].items():
        d = find_oxide_dir(root, name)
        if d is not None:
            targets.append((name, d, meta))

    print(f"regression: {len(targets)} cases, beta={beta}, eneg={eneg_scale}")
    import time
    for name, d, meta in targets:
        t0 = time.time()
        print(f"[{name}] ...", end=" ")
        try:
            case = load_case(d)
            res = assign_case(case, beta=beta, eneg_scale=eneg_scale)
            print(f"done in {time.time()-t0:.1f}s")
        except Exception as exc:                        # keep the sweep going
            print(f"ERROR after {time.time()-t0:.1f}s: {str(exc)[:60]}")
            rows.append((name, "ERROR", str(exc)[:60], "", ""))
            continue
        st = case.structure
        flow_sum = element_summary(st.elements, res.os_int)
        greedy_sum = (element_summary(st.elements, case.eoscar.greedy_os)
                      if case.eoscar else {})
        ref_map = meta.get("ref") or meta.get("ref_mean") or {}
        dev_flow = dev_greedy = 0.0
        for el, target in ref_map.items():
            if el in flow_sum:
                dev_flow = max(dev_flow, abs(flow_sum[el][2] - target))
            if el in greedy_sum:
                dev_greedy = max(dev_greedy, abs(greedy_sum[el][2] - target))
        status = "OK" if ref_map and dev_flow <= 0.25 else (
            "~" if ref_map and dev_flow <= 1.0 else ("FAIL" if ref_map else "run"))
        flags = []
        if res.stray:
            flags.append(f"stray:{len(res.stray)}")
        if res.infeasible:
            flags.append(f"infeas:{res.infeasible}")
        rows.append((name, status, f"{dev_flow:.2f}" if ref_map else "-",
                     f"{dev_greedy:.2f}" if ref_map else "-", ",".join(flags)))

    print(f"\n{'case':<26} {'flow':<6} {'dev_flow':>8} {'dev_greedy':>10}  flags")
    n_ok = n_better = n_worse = 0
    for name, status, df, dg, flags in rows:
        print(f"{name:<26} {status:<6} {df:>8} {dg:>10}  {flags}")
        if status == "OK":
            n_ok += 1
        try:
            dfv, dgv = float(df), float(dg)
            if dfv < dgv - 0.05:
                n_better += 1
            elif dfv > dgv + 0.05:
                n_worse += 1
        except ValueError:
            pass
    print(f"\nsummary: {n_ok}/{len(rows)} within 0.25 of reference; "
          f"transport better than greedy in {n_better}, worse in {n_worse}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("case")
    p1.add_argument("dir", type=Path)
    p1.add_argument("--beta", type=float, default=10.0)
    p1.add_argument("--eneg", default="allen", choices=["allen", "pauling"])
    p1.add_argument("-v", "--verbose", action="store_true")
    p1.add_argument("--drop-stray", action="store_true",
                    help="exclude unlocalised WFs (attribution diagnostic)")
    p1.add_argument("--gamma", type=float, default=0.0,
                    help="tier-2 projector-composition weight (needs #WP data)")
    p2 = sub.add_parser("regress")
    p2.add_argument("--root", type=Path, default=Path("leos_work"))
    p2.add_argument("--beta", type=float, default=10.0)
    p2.add_argument("--eneg", default="allen", choices=["allen", "pauling"])
    args = ap.parse_args()
    if args.cmd == "case":
        run_case(args.dir, args.beta, args.eneg, args.verbose, args.drop_stray,
                 args.gamma)
    else:
        regress(args.root, args.beta, args.eneg)


if __name__ == "__main__":
    main()
