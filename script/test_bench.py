#!/usr/bin/env python3
"""Score assignment plugins on the pre-registered held-out test sets.

The sets and their labels are defined in pyleos/tests/test_set.json and were
committed before any test calculation existed:

  regression    the 12 clean-set cases whose archived Wannier functions had
                stray WFs, recomputed with gpuw90. Element references come
                from reference_os.json and site specs from clean_set.json.
                Cases whose recomputed Wannier set has an unconverged spin
                channel are scored advisory-only (reported, not counted).
  homonuclear   13 new compounds with homonuclear bonds; references and site
                specs are in test_set.json itself.

Each case is read from <root>/<case>/ (POSCAR, EOSCAR and the wannier90 wout
file(s)), or from the case's registered "dir" when --root is not given. The
plugin sees exactly what the dev evaluator shows it (script/clean_bench.py:
same static check, same table trackers, same scrub). For regression cases the
EOSCAR #WF block is dropped as well: its per-WF occupations belong to the
archived gauge, not to the recomputed WFs, so the plugin gets the recomputed
centres and spreads and NELECT(Wannier) only.

Usage: python3 script/test_bench.py --set regression|homonuclear|both
       --method PLUGIN.py [--root DIR] [--json FILE]
       python3 script/test_bench.py --selftest     (checks the scoring path
       against clean_bench on dev cases; reads no test data)
"""
import argparse
import json
import sys
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "script"))

import clean_bench as cb                                        # noqa: E402
from pyleos.parsers import load_case                            # noqa: E402

TEST = json.loads((ROOT / "pyleos/tests/test_set.json").read_text())


def cases(which):
    """Yield (set, batch, name, ref, site_spec, registered_dir, advisory)."""
    if which in ("regression", "both"):
        refs = {(b, n): ref for b, n, _d, ref in cb.collect()}
        status = TEST["test_regression"].get("wannier_status", {})
        for name, meta in TEST["test_regression"]["cases"].items():
            b = meta["batch"]
            spec = cb.CLEAN["cases"][b][name]
            yield ("regression", b, name, refs[(b, name)], spec, meta.get("dir"),
                   status.get(name, "converged") != "converged")
    if which in ("homonuclear", "both"):
        for name, meta in TEST["test_homonuclear"]["cases"].items():
            yield ("homonuclear", "test", name, meta["ref"], meta["site"],
                   meta.get("dir"), False)


def score_case(plugin, case, name, ref, spec, drop_wf_detail):
    st = case.structure
    labels = cb.site_labels(name, st, ref, spec)
    n_lab = sum(lab is not None for lab in labels)
    row = dict(status="ERROR", dev=None, atoms_ok=0, atoms_labelled=n_lab,
               wrong="", seconds=0.0, error="")
    case = cb.scrub(case)
    if drop_wf_detail and case.eoscar is not None:
        case.eoscar.wf_detail = []
    t0 = time.time()
    try:
        os_vals = [int(v) for v in plugin.assign(case)]
        if len(os_vals) != len(st.elements):
            raise ValueError(f"returned {len(os_vals)} values for {len(st.elements)} atoms")
    except Exception as exc:                           # scored as a failure
        row.update(error=f"{type(exc).__name__}: {exc}"[:200], seconds=round(time.time() - t0, 2))
        return row, None
    row["seconds"] = round(time.time() - t0, 2)
    means = cb.element_means(st.elements, os_vals)
    dev = max(abs(means[el] - t) for el, t in ref.items() if el in means)
    wrong = Counter((el, lab, o) for el, lab, o in zip(st.elements, labels, os_vals)
                    if lab is not None and o != lab)
    row.update(status="OK" if dev <= 0.25 else "~" if dev <= 1.0 else "FAIL",
               dev=round(dev, 4), atoms_ok=n_lab - sum(wrong.values()),
               wrong=",".join(f"{el}{lab:+d}->{o:+d}x{k}" for (el, lab, o), k in sorted(wrong.items())))
    return row, os_vals


def selftest():
    """Score three dev cases through score_case and compare with clean_bench."""
    plugin = cb.load_plugin(cb.DEFAULT_METHOD)
    probe = {("wave2", "h2s2"), ("wave2", "rb2cucl4"), ("v1_mol", "CH2F2")}
    ok = True
    for b, n, d, ref in cb.collect():
        if (b, n) not in probe:
            continue
        row, _ = score_case(plugin, load_case(d), n, ref, cb.CLEAN["cases"][b][n], False)
        print(f"  {b}/{n}: {row['status']} atoms {row['atoms_ok']}/{row['atoms_labelled']} {row['wrong']}")
        ok &= row["status"] != "ERROR"
    print("selftest", "passed" if ok else "FAILED", "(compare with: python3 script/clean_bench.py --split dev)")
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--set", choices=["regression", "homonuclear", "both"], default="both")
    ap.add_argument("--method", type=Path, default=cb.DEFAULT_METHOD)
    ap.add_argument("--root", type=Path, help="directory holding one subdirectory per case")
    ap.add_argument("--json", type=Path)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        sys.exit(0 if selftest() else 1)

    violations, literals = cb.static_check(args.method.read_text())
    if violations:
        sys.exit(f"{args.method}: forbidden constructs {violations}")
    cb.install_trackers()
    plugin = cb.load_plugin(args.method)

    rows = []
    for tset, batch, name, ref, spec, reg_dir, advisory in cases(args.set):
        d = args.root / name if args.root else (ROOT / reg_dir if reg_dir else None)
        if d is None or not (d / "POSCAR").is_file():
            sys.exit(f"{tset}/{name}: no data (dir {d}); a case is scored only once its data has landed")
        row, os_vals = score_case(plugin, load_case(d), name, ref, spec, tset == "regression")
        row.update(set=tset, batch=batch, case=name, advisory=advisory, os=os_vals)
        rows.append(row)

    for tset in ("regression", "homonuclear"):
        part = [r for r in rows if r["set"] == tset]
        if not part:
            continue
        print(f"== {tset}")
        for r in part:
            dev = "-" if r["dev"] is None else f"{r['dev']:.2f}"
            print(f"  {r['case']:8s} {'advisory' if r['advisory'] else 'scored  '} {r['status']:5s} dev={dev} "
                  f"atoms {r['atoms_ok']}/{r['atoms_labelled']} {r['wrong'] or r['error']}")
        gate = [r for r in part if not r["advisory"]]
        cb.print_summary(cb.summary(gate), f"{tset} (scored cases)")
        if len(gate) < len(part):
            cb.print_summary(cb.summary(part), f"{tset} (all cases incl. advisory)")
    flags = {"tables_used": sorted(cb.TABLES_USED), "element_literals": literals,
             "meta": getattr(plugin, "META", None)}
    print(f"tables used: {flags['tables_used'] or 'none'}")
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        out = {"method": str(args.method), "set": args.set, "flags": flags, "cases": rows,
               "summary": {t: {"scored": cb.summary([r for r in rows if r["set"] == t and not r["advisory"]]),
                               "all": cb.summary([r for r in rows if r["set"] == t])}
                           for t in ("regression", "homonuclear") if any(r["set"] == t for r in rows)}}
        args.json.write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
