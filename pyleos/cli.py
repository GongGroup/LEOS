#!/usr/bin/env python3
"""pyleos CLI: the single entry point (python3 -m pyleos, or python3 pyleos/cli.py).

Usage:
  python3 -m pyleos case <case_dir> [--beta B] [--eneg allen|pauling] [-v]
  python3 -m pyleos regress [--beta B] [--eneg SCALE] [--reliability]
  python3 -m pyleos bench full|clean|test ...     (see pyleos/bench.py)
  python3 -m pyleos batch stage|submit|status|... (see pyleos/batch.py)
  python3 -m pyleos check                         (cases/INDEX.tsv vs STRUCTURES.tsv vs disk)

regress runs the 36 v1 cases (bench sets v1_mol and v1_oxide of
cases/INDEX.tsv) and compares transport with the in-situ greedy column.
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
    from pyleos.cases import bench74_dirs, check_registry           # noqa: E402
    from pyleos.parsers import load_case                            # noqa: E402
    from pyleos.transport import assign_case, reliability_margins   # noqa: E402
else:
    from .cases import bench74_dirs, check_registry
    from .parsers import load_case
    from .transport import assign_case, reliability_margins

REF_PATH = Path(__file__).resolve().parent / "sets" / "reference_os.json"


def element_summary(elements, values):
    by_el = defaultdict(list)
    for el, v in zip(elements, values):
        by_el[el].append(v)
    return {el: (min(v), max(v), sum(v) / len(v)) for el, v in sorted(by_el.items())}


def run_case(case_dir: Path, beta: float, eneg_scale: str, verbose: bool = False,
             drop_stray: bool = False, gamma: float = 0.0,
             fragments=None):
    case = load_case(case_dir)
    res = assign_case(case, beta=beta, eneg_scale=eneg_scale, drop_stray=drop_stray,
                      gamma=gamma, fragments=fragments)
    res = reliability_margins(case, beta=beta, eneg_scale=eneg_scale,
                              drop_stray=drop_stray, gamma=gamma, res=res,
                              fragments=fragments)
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
    if res.atom_margin:
        print(f"reliability R_elem = {res.case_conf_elem:.0f}/100 "
              f"(element ledger)   R_site = {res.case_conf:.0f}/100 "
              f"(site level; low values flag intra-element degeneracy "
              f"such as homonuclear bonds)")
    for i, members in enumerate(res.fragments):
        names = "+".join(f"{st.elements[a]}{a+1}" for a in members)
        print(f"fragment {i+1} [{names}]: ne={res.frag_ne[i]} "
              f"total OS {res.frag_os[i]:+d}")
    if not res.fragments:
        # Fragment suggestion: electron pairs with a near-tie alternative
        # host are the fingerprint of delocalised multi-centre units.
        tied = sorted((m, w) for w, m in enumerate(res.wf_margin) if m < 0.5)
        if tied:
            ws = ", ".join(f"WF{w+1}(margin {m:.2f})" for m, w in tied[:3])
            print(f"note: {len(tied)} electron pair(s) sit at a near-tie "
                  f"between hosts ({ws}) - if these span a delocalised "
                  f"multi-centre unit, consider --fragment")
    if verbose:
        print(f"{'ion':>4} {'el':>3} {'ne':>4} {'OS_flow':>8} {'OS_greedy':>10} {'R':>5}")
        for a, el in enumerate(st.elements):
            g = case.eoscar.greedy_os[a] if case.eoscar else None
            osa = res.os_int[a]
            print(f"{a+1:>4} {el:>3} {res.ne[a]:>4} "
                  f"{'frag' if osa is None else format(osa, '+d'):>8} "
                  f"{g if g is None else format(g, '+d'):>10} "
                  f"{res.atom_conf[a]:>5.0f}")
    kept = [(el, v) for el, v in zip(st.elements, res.os_int) if v is not None]
    flow_sum = element_summary([e for e, _ in kept], [v for _, v in kept])
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


def regress(beta: float, eneg_scale: str, reliability: bool = False):
    ref = json.loads(REF_PATH.read_text())
    dirs = bench74_dirs()
    rows = []
    targets = []
    for bench_set, sub in (("v1_mol", "dataset"), ("v1_oxide", "oxides")):
        for name, meta in ref[sub].items():
            d = dirs.get((bench_set, name))
            if d is not None and (d / "EOSCAR").is_file():
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
        if reliability:
            res = reliability_margins(case, beta=beta, eneg_scale=eneg_scale,
                                      res=res)
        flags = []
        if res.stray:
            flags.append(f"stray:{len(res.stray)}")
        if res.infeasible:
            flags.append(f"infeas:{res.infeasible}")
        if res.atom_margin:
            flags.append(f"R:{res.case_conf:.0f}")
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


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv[:1] == ["bench"]:
        from pyleos import bench
        return bench.main(argv[1:])
    if argv[:1] == ["batch"]:
        from pyleos import batch
        sys.argv = ["pyleos batch"] + argv[1:]
        return batch.main()
    ap = argparse.ArgumentParser(prog="pyleos", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("bench", help="benchmarks: python3 -m pyleos bench -h")
    sub.add_parser("batch", help="cluster batch driver: python3 -m pyleos batch -h")
    p1 = sub.add_parser("case")
    p1.add_argument("dir", type=Path)
    p1.add_argument("--beta", type=float, default=10.0)
    p1.add_argument("--eneg", default="allen", choices=["allen", "pauling"])
    p1.add_argument("-v", "--verbose", action="store_true")
    p1.add_argument("--drop-stray", action="store_true",
                    help="exclude unlocalised WFs (attribution diagnostic)")
    p1.add_argument("--gamma", type=float, default=0.0,
                    help="tier-2 projector-composition weight (needs #WP data)")
    p1.add_argument("--fragment", action="append", default=[],
                    help="comma-separated 1-based atom indices to merge into "
                         "one delocalised unit (repeatable), e.g. --fragment 2,3,4")
    sub.add_parser("check", help="consistency of cases/INDEX.tsv, STRUCTURES.tsv and cases/")
    p2 = sub.add_parser("regress")
    p2.add_argument("--beta", type=float, default=10.0)
    p2.add_argument("--eneg", default="allen", choices=["allen", "pauling"])
    p2.add_argument("--reliability", action="store_true",
                    help="compute the +-1 re-solve reliability R per case")
    args = ap.parse_args(argv)
    if args.cmd == "case":
        frags = [[int(t) - 1 for t in f.split(",")] for f in args.fragment]
        run_case(args.dir, args.beta, args.eneg, args.verbose, args.drop_stray,
                 args.gamma, fragments=frags or None)
    elif args.cmd == "check":
        problems = check_registry()
        for line in problems:
            print(line)
        print(f"{len(problems)} problem(s)")
        sys.exit(1 if problems else 0)
    else:
        regress(args.beta, args.eneg, args.reliability)


if __name__ == "__main__":
    main()
