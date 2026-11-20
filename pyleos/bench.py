#!/usr/bin/env python3
"""Benchmarks: the 74-case table, the clean set and the held-out test sets.

  python3 -m pyleos bench full  --out FILE [--beta B] [--gamma G]
                                [--tier2 T] [--only-wp] [--no-reliability] [--greedy] [--corrected]
  python3 -m pyleos bench clean [--split all|dev] [--method PLUGIN.py]
                                [--json FILE] [--compare BASELINE.json] [--out FILE] [--greedy]
  python3 -m pyleos bench test  --set regression|homonuclear|both
                                --method PLUGIN.py [--root DIR] [--json FILE]
  python3 -m pyleos bench test  --selftest

full   Transport OS + two-level reliability R on the 74 cases (four sets:
       v1 dataset, v1 oxides, wave1, wave2) against sets/reference_os.json.
       Writes a TSV with per-case dev and the R_elem / R_site reliability
       scores, plus a summary of how R_elem separates passing from failing
       cases. results/bench74_golden.tsv is this output; CI diffs columns 1-4
       against it. --out is required so a routine run cannot overwrite it.
       --corrected scores the cases of sets/corrected.json on their
       structure-corrected runs (results/bench74_corrected.tsv).

clean  Element-mean dev plus per-atom accuracy on the cases selected in
       sets/clean_set.json. Element values come from reference_os.json,
       per-atom labels from clean_set.json; the labels are checked to
       reproduce the reference element means before scoring. Writes no file
       unless --out / --json is given.

test   Scores a plugin on the pre-registered held-out test sets of
       sets/test_set.json, committed before any test calculation existed:
         regression   the 12 clean-set cases whose archived Wannier functions
                      had stray WFs, recomputed with gpuw90. Element
                      references come from reference_os.json and site specs
                      from clean_set.json. Cases whose recomputed Wannier set
                      has an unconverged spin channel are scored advisory-only
                      (reported, not counted).
         homonuclear  13 new compounds with homonuclear bonds; references and
                      site specs are in test_set.json itself.
       Each case is read from <root>/<case>/ (POSCAR, EOSCAR and the wannier90
       wout file(s)), or from the case's registered "dir" when --root is not
       given. The plugin sees exactly what the clean evaluator shows it (same
       static check, same table trackers, same scrub). For regression cases the
       EOSCAR #WF block is dropped as well: its per-WF occupations belong to the
       archived gauge, not to the recomputed WFs, so the plugin gets the
       recomputed centres and spreads and NELECT(Wannier) only.

Assignment methods (clean, test) are plugins: a .py file defining
assign(case) -> list of integer oxidation states in POSCAR atom order
(optional META dict). The default plugin pyleos/methods/baseline.py is the
shipped transport. Before a plugin sees a case, outputs of other methods
(EOSCAR greedy/flow columns, greedy host fields of #WF) and the #WP block
(invalid in archived data, see patch/README.md) are removed, and the case
path is anonymised. The run records whether the plugin touched the
electronegativity or the oxidation-state bound tables.
"""
from __future__ import annotations

import argparse
import dataclasses
import functools
import importlib.util
import itertools
import json
import math
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import pyleos.capacity as capacity                              # noqa: E402
import pyleos.transport as transport                            # noqa: E402
from pyleos.cases import bench74_dirs                           # noqa: E402
from pyleos.parsers import load_case                            # noqa: E402
from pyleos.transport import assign_case, reliability_margins   # noqa: E402

print = functools.partial(print, flush=True)
SETS = ROOT / "pyleos" / "sets"
REF = json.loads((SETS / "reference_os.json").read_text())
CLEAN = json.loads((SETS / "clean_set.json").read_text())
TEST = json.loads((SETS / "test_set.json").read_text())
DEFAULT_METHOD = ROOT / "pyleos" / "methods" / "baseline.py"


# ================================================================ full (74 cases)
def collect(root: Path = ROOT):
    """Yield (batch, name, dir, ref_map) for the 74 cases, in reference order.
    root: alternate tree mirroring the repo layout (cases/INDEX.tsv + cases/)."""
    dirs = bench74_dirs(root)
    for batch, sub in (("v1_mol", "dataset"), ("v1_oxide", "oxides"),
                       ("wave1", "wave1"), ("wave2", "wave2")):
        for name, meta in REF[sub].items():
            if batch.startswith("wave") and name.startswith("_"):
                continue                            # wave-level _comment keys
            d = dirs.get((batch, name))
            if d is not None and d.is_dir() and (d / "EOSCAR").is_file():
                ref = (meta.get("ref") or meta.get("ref_mean") or {}) \
                    if isinstance(meta, dict) else {}
                yield batch, name, d, ref


def element_means(elements, values):
    by = defaultdict(list)
    for el, v in zip(elements, values):
        by[el].append(v)
    return {el: sum(v) / len(v) for el, v in by.items()}


def main_full(argv=None):
    ap = argparse.ArgumentParser(prog="pyleos bench full",
                                 description="74-case benchmark (see pyleos/bench.py)")
    ap.add_argument("--out", type=Path, required=True,
                    help="output TSV (results/bench74_golden.tsv only when "
                         "deliberately regenerating the golden table)")
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
    ap.add_argument("--corrected", action="store_true",
                    help="score the cases listed in sets/corrected.json on their "
                         "structure-corrected runs instead (adds a run column; "
                         "never compare this output with the golden table)")
    args = ap.parse_args(argv)

    rows = []
    t_start = time.time()
    cases = list(collect())
    if args.corrected:
        fix = json.loads((SETS / "corrected.json").read_text())["cases"]
        cases = [(b, n, ROOT / "cases" / fix[f"{b}/{n}"]["run"] if f"{b}/{n}" in fix else d, ref)
                 for b, n, d, ref in cases]
    run_of = {(b, n): str(d.relative_to(ROOT / "cases")) for b, n, d, _ in cases}
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
            rows.append((batch, name, "ERROR", "", "", "", "", "") + ((run_of[(batch, name)],) if args.corrected else ()))
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
                     res.ghost if res is not None else "-")
                    + ((run_of[(batch, name)],) if args.corrected else ()))
        print(f"{status} dev={rows[-1][3]} R_elem={rows[-1][4]} "
              f"({time.time()-t0:.1f}s)")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        fh.write("batch\tcase\tstatus\tdev\tR_elem\tR_site\tn_stray\tn_ghost"
                 + ("\trun" if args.corrected else "") + "\n")
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


# ================================================================ clean set
def collect_clean():
    """collect() with the data directories that clean_set.json overrides
    (reruns that fix a clean-set case without touching the 74-case
    benchmark data)."""
    over = {(b, n): ROOT / d for b, cases in CLEAN.get("overrides", {}).items()
            if not b.startswith("_") for n, d in cases.items()}
    for batch, name, d, ref in collect():
        yield batch, name, over.get((batch, name), d), ref

# ---------------------------------------------------------------- tracking
TABLES_USED = set()


class _TrackedTable(dict):
    def __init__(self, data, tag):
        super().__init__(data)
        self._tag = tag

    def _hit(self):
        TABLES_USED.add(self._tag)

    def __getitem__(self, k):
        self._hit(); return super().__getitem__(k)

    def get(self, k, default=None):
        self._hit(); return super().get(k, default)

    def __contains__(self, k):
        self._hit(); return super().__contains__(k)

    def __iter__(self):
        self._hit(); return super().__iter__()

    def keys(self):
        self._hit(); return super().keys()

    def values(self):
        self._hit(); return super().values()

    def items(self):
        self._hit(); return super().items()


def install_trackers():
    capacity.ENEG_ALLEN = _TrackedTable(capacity.ENEG_ALLEN, "electronegativity")
    capacity.ENEG_PAULING = _TrackedTable(capacity.ENEG_PAULING, "electronegativity")
    capacity.OS_BOUNDS = _TrackedTable(capacity.OS_BOUNDS, "os_bounds")
    orig = capacity.capacity_bounds

    @functools.wraps(orig)
    def tracked_bounds(*a, **k):
        TABLES_USED.add("os_bounds")
        return orig(*a, **k)
    capacity.capacity_bounds = tracked_bounds
    transport.capacity_bounds = tracked_bounds


# ---------------------------------------------------------------- plugins
STATIC_FORBIDDEN = ["open(", "Path(", "import os", "from os", "subprocess",
                    "reference_os", "clean_set", "test_set", "__file__", "eval(",
                    "exec(", "importlib", "glob", "json.load", "sys.path",
                    "socket", "urllib", "pickle"]
ELEMENTS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn "
            "Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag "
            "Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm "
            "Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Th U").split()
ELEMENT_LITERAL = re.compile(r"""['"](%s)['"]""" % "|".join(ELEMENTS))


def static_check(src: str):
    code = "\n".join(line.split("#", 1)[0] for line in src.splitlines())
    violations = [tok for tok in STATIC_FORBIDDEN if tok in code]
    literals = sorted(set(ELEMENT_LITERAL.findall(code)))
    return violations, literals


def load_plugin(path: Path):
    spec = importlib.util.spec_from_file_location("leos_candidate", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    if not callable(getattr(mod, "assign", None)):
        sys.exit(f"{path}: plugin must define assign(case)")
    return mod


def scrub(case):
    """Remove other methods' answers and invalid data; anonymise the case."""
    case.path = Path("case")
    e = case.eoscar
    if e is not None:
        e.greedy_ne, e.greedy_os, e.flow_os = [], [], None
        e.wp_detail = {}
        e.wf_detail = [dataclasses.replace(w, at1=0, e1=0.0, at2=0, e2=0.0)
                       for w in e.wf_detail]
    return case


# ---------------------------------------------------------------- labels
def distance(st, i: int, j: int) -> float:
    df = [st.frac[i][k] - st.frac[j][k] for k in range(3)]
    df = [x - round(x) for x in df]
    best = math.inf
    for n in itertools.product((-1, 0, 1), repeat=3):
        f = [df[k] + n[k] for k in range(3)]
        cart = [sum(f[k] * st.lattice[k][c] for k in range(3)) for c in range(3)]
        best = min(best, math.sqrt(sum(x * x for x in cart)))
    return best


def site_labels(name, st, ref, spec):
    """Per-atom reference OS (None where only the element mean is defined)."""
    labels = []
    for i, el in enumerate(st.elements):
        s = spec.get(el, ref.get(el))
        if s is None or s == "mean":
            labels.append(None)
        elif isinstance(s, dict):
            n = sum(1 for j, e in enumerate(st.elements)
                    if j != i and e == s["count"] and distance(st, i, j) < s["within"])
            if str(n) not in s["map"]:
                raise ValueError(f"{name}: {el}{i + 1} has {n} {s['count']} "
                                 f"neighbours within {s['within']} A, not in map")
            labels.append(s["map"][str(n)])
        elif float(s).is_integer():
            labels.append(int(s))
        else:
            raise ValueError(f"{name}: {el} reference {s} is not an integer; "
                             f"give a site rule or \"mean\" in clean_set.json")
    by = defaultdict(list)
    for el, lab in zip(st.elements, labels):
        if lab is not None:
            by[el].append(lab)
    for el, labs in by.items():
        if el in ref and abs(sum(labs) / len(labs) - ref[el]) > 1e-3:
            raise ValueError(f"{name}: {el} site labels average "
                             f"{sum(labs) / len(labs):.4f}, reference {ref[el]}")
    return labels


# ---------------------------------------------------------------- scoring
def summary(rows):
    st = Counter(r["status"] for r in rows)
    ok = sum(r["atoms_ok"] for r in rows)
    n = sum(r["atoms_labelled"] for r in rows)
    return {"cases": len(rows), "elem_ok": st["OK"], "elem_mid": st["~"],
            "elem_fail": st["FAIL"] + st["ERROR"],
            "all_atoms_correct_cases": sum(r["atoms_ok"] == r["atoms_labelled"]
                                           and r["status"] != "ERROR" for r in rows),
            "atoms_ok": ok, "atoms_labelled": n,
            "atom_accuracy": round(ok / n, 4) if n else None}


def print_summary(s, title):
    if s["cases"]:
        print(f"{title}: {s['cases']} cases | element mean OK {s['elem_ok']}, "
              f"~ {s['elem_mid']}, FAIL {s['elem_fail']} | all labelled atoms "
              f"correct {s['all_atoms_correct_cases']}/{s['cases']} | atoms "
              f"{s['atoms_ok']}/{s['atoms_labelled']} ({100 * s['atom_accuracy']:.1f}%)")


def main_clean(argv=None):
    ap = argparse.ArgumentParser(prog="pyleos bench clean",
                                 description="clean-set benchmark (see pyleos/bench.py)")
    ap.add_argument("--split", choices=["all", "dev"], default="all")
    ap.add_argument("--method", type=Path, default=DEFAULT_METHOD,
                    help="assignment plugin (.py defining assign(case))")
    ap.add_argument("--greedy", action="store_true",
                    help="score the in-situ greedy OS column of EOSCAR instead")
    ap.add_argument("--json", type=Path, help="write machine-readable results")
    ap.add_argument("--compare", type=Path,
                    help="results JSON of a reference run; report per-case regressions")
    ap.add_argument("--out", type=Path, help="write a per-case TSV")
    args = ap.parse_args(argv)

    wanted = {(b, n): spec for b, cases in CLEAN["cases"].items()
              for n, spec in cases.items()}
    excluded = {(e["batch"], e["case"]) for e in CLEAN["excluded"]}
    if set(wanted) & excluded:
        sys.exit(f"cases both selected and excluded: {sorted(set(wanted) & excluded)}")
    held_out = {tuple(x) for x in CLEAN["splits"]["test_regression"]}
    if args.split == "dev":
        wanted = {k: v for k, v in wanted.items() if k not in held_out}

    src = args.method.read_text()
    violations, literals = static_check(src)
    if not args.greedy and violations:
        sys.exit(f"{args.method}: forbidden constructs {violations} "
                 f"(plugins get the case object only; no file or system access)")
    install_trackers()
    plugin = None if args.greedy else load_plugin(args.method)

    rows, seen = [], set()
    for batch, name, d, ref in collect_clean():
        if (batch, name) not in wanted:
            if (batch, name) not in excluded and (batch, name) not in held_out:
                sys.exit(f"{batch}/{name} is neither selected nor excluded")
            continue
        seen.add((batch, name))
        case = load_case(d)
        st = case.structure
        labels = site_labels(name, st, ref, wanted[(batch, name)])
        n_lab = sum(lab is not None for lab in labels)
        row = dict(batch=batch, case=name, status="ERROR", dev=None, atoms_ok=0,
                   atoms_labelled=n_lab, wrong="", seconds=0.0, error="")
        t0 = time.time()
        try:
            if args.greedy:
                os_vals = list(case.eoscar.greedy_os)
            else:
                os_vals = plugin.assign(scrub(case))
            os_vals = [int(v) for v in os_vals]
            if len(os_vals) != len(st.elements):
                raise ValueError(f"returned {len(os_vals)} values for "
                                 f"{len(st.elements)} atoms")
        except Exception as exc:                       # scored as a failure
            row["error"] = f"{type(exc).__name__}: {exc}"[:200]
            row["seconds"] = round(time.time() - t0, 2)
            rows.append(row)
            continue
        row["seconds"] = round(time.time() - t0, 2)
        means = element_means(st.elements, os_vals)
        dev = max(abs(means[el] - t) for el, t in ref.items() if el in means)
        wrong = Counter((el, lab, o) for el, lab, o in zip(st.elements, labels, os_vals)
                        if lab is not None and o != lab)
        row.update(status="OK" if dev <= 0.25 else "~" if dev <= 1.0 else "FAIL",
                   dev=round(dev, 4), atoms_ok=n_lab - sum(wrong.values()),
                   wrong=",".join(f"{el}{lab:+d}->{o:+d}x{k}"
                                  for (el, lab, o), k in sorted(wrong.items())))
        rows.append(row)
    missing = set(wanted) - seen
    if missing:
        sys.exit(f"selected cases without data: {sorted(missing)}")

    for r in rows:
        if r["status"] != "OK" or r["atoms_ok"] != r["atoms_labelled"]:
            dev = "-" if r["dev"] is None else f"{r['dev']:.2f}"
            print(f"  {r['batch']:8s} {r['case']:12s} {r['status']:5s} dev={dev} "
                  f"atoms {r['atoms_ok']}/{r['atoms_labelled']} {r['wrong'] or r['error']}")
    s = summary(rows)
    print_summary(s, f"clean set ({args.split})")
    flags = {"tables_used": sorted(TABLES_USED), "element_literals": literals,
             "meta": getattr(plugin, "META", None) if plugin else None}
    print(f"tables used: {flags['tables_used'] or 'none'}"
          + (f" | element literals in source: {literals}" if literals else ""))

    regressions = []
    if args.compare:
        ref_rows = {(r["batch"], r["case"]): r
                    for r in json.loads(args.compare.read_text())["cases"]}
        rank = {"OK": 0, "~": 1, "FAIL": 2, "ERROR": 3}
        for r in rows:
            b = ref_rows.get((r["batch"], r["case"]))
            if b and (r["atoms_ok"] < b["atoms_ok"] or rank[r["status"]] > rank[b["status"]]):
                regressions.append({"batch": r["batch"], "case": r["case"],
                                    "atoms_ok": [b["atoms_ok"], r["atoms_ok"]],
                                    "status": [b["status"], r["status"]]})
        print(f"regressions vs {args.compare.name}: {len(regressions)}"
              + "".join(f"\n  {g['case']}: atoms {g['atoms_ok'][0]}->{g['atoms_ok'][1]}, "
                        f"{g['status'][0]}->{g['status'][1]}" for g in regressions))

    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(
            {"method": "greedy" if args.greedy else str(args.method), "split": args.split,
             "summary": s, "flags": flags, "regressions": regressions, "cases": rows},
            indent=1) + "\n")
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open("w") as fh:
            fh.write("batch\tcase\tstatus\tdev\tatoms_ok\tatoms_labelled\twrong_sites\n")
            for r in rows:
                dev = "-" if r["dev"] is None else f"{r['dev']:.2f}"
                fh.write(f"{r['batch']}\t{r['case']}\t{r['status']}\t{dev}\t"
                         f"{r['atoms_ok']}\t{r['atoms_labelled']}\t{r['wrong'] or '-'}\n")


# ================================================================ held-out test sets
def test_cases(which):
    """Yield (set, batch, name, ref, site_spec, registered_dir, advisory)."""
    if which in ("regression", "both"):
        refs = {(b, n): ref for b, n, _d, ref in collect_clean()}
        status = TEST["test_regression"].get("wannier_status", {})
        for name, meta in TEST["test_regression"]["cases"].items():
            b = meta["batch"]
            spec = CLEAN["cases"][b][name]
            yield ("regression", b, name, refs[(b, name)], spec, meta.get("dir"),
                   status.get(name, "converged") != "converged")
    if which in ("homonuclear", "both"):
        for name, meta in TEST["test_homonuclear"]["cases"].items():
            yield ("homonuclear", "test", name, meta["ref"], meta["site"],
                   meta.get("dir"), False)


def score_case(plugin, case, name, ref, spec, drop_wf_detail):
    st = case.structure
    labels = site_labels(name, st, ref, spec)
    n_lab = sum(lab is not None for lab in labels)
    row = dict(status="ERROR", dev=None, atoms_ok=0, atoms_labelled=n_lab,
               wrong="", seconds=0.0, error="")
    case = scrub(case)
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
    means = element_means(st.elements, os_vals)
    dev = max(abs(means[el] - t) for el, t in ref.items() if el in means)
    wrong = Counter((el, lab, o) for el, lab, o in zip(st.elements, labels, os_vals)
                    if lab is not None and o != lab)
    row.update(status="OK" if dev <= 0.25 else "~" if dev <= 1.0 else "FAIL",
               dev=round(dev, 4), atoms_ok=n_lab - sum(wrong.values()),
               wrong=",".join(f"{el}{lab:+d}->{o:+d}x{k}" for (el, lab, o), k in sorted(wrong.items())))
    return row, os_vals


def selftest():
    """Score three dev cases through score_case and compare with the clean bench."""
    plugin = load_plugin(DEFAULT_METHOD)
    probe = {("wave2", "h2s2"), ("wave2", "rb2cucl4"), ("v1_mol", "CH2F2")}
    ok = True
    for b, n, d, ref in collect_clean():
        if (b, n) not in probe:
            continue
        row, _ = score_case(plugin, load_case(d), n, ref, CLEAN["cases"][b][n], False)
        print(f"  {b}/{n}: {row['status']} atoms {row['atoms_ok']}/{row['atoms_labelled']} {row['wrong']}")
        ok &= row["status"] != "ERROR"
    print("selftest", "passed" if ok else "FAILED", "(compare with: python3 -m pyleos bench clean --split dev)")
    return ok


def main_test(argv=None):
    ap = argparse.ArgumentParser(prog="pyleos bench test",
                                 description="held-out test sets (see pyleos/bench.py)")
    ap.add_argument("--set", choices=["regression", "homonuclear", "both"], default="both")
    ap.add_argument("--method", type=Path, default=DEFAULT_METHOD)
    ap.add_argument("--root", type=Path, help="directory holding one subdirectory per case")
    ap.add_argument("--json", type=Path)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args(argv)
    if args.selftest:
        sys.exit(0 if selftest() else 1)

    violations, literals = static_check(args.method.read_text())
    if violations:
        sys.exit(f"{args.method}: forbidden constructs {violations}")
    install_trackers()
    plugin = load_plugin(args.method)

    rows = []
    for tset, batch, name, ref, spec, reg_dir, advisory in test_cases(args.set):
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
        print_summary(summary(gate), f"{tset} (scored cases)")
        if len(gate) < len(part):
            print_summary(summary(part), f"{tset} (all cases incl. advisory)")
    flags = {"tables_used": sorted(TABLES_USED), "element_literals": literals,
             "meta": getattr(plugin, "META", None)}
    print(f"tables used: {flags['tables_used'] or 'none'}")
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        out = {"method": str(args.method), "set": args.set, "flags": flags, "cases": rows,
               "summary": {t: {"scored": summary([r for r in rows if r["set"] == t and not r["advisory"]]),
                               "all": summary([r for r in rows if r["set"] == t])}
                           for t in ("regression", "homonuclear") if any(r["set"] == t for r in rows)}}
        args.json.write_text(json.dumps(out, indent=1) + "\n")


COMMANDS = {"full": main_full, "clean": main_clean, "test": main_test}


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv or argv[0] not in COMMANDS:
        sys.exit(__doc__)
    COMMANDS[argv[0]](argv[1:])


if __name__ == "__main__":
    main()
