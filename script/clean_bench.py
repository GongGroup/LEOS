#!/usr/bin/env python3
"""Clean-set benchmark: element-mean dev plus per-atom accuracy.

Scores the cases selected in pyleos/tests/clean_set.json. Element values come
from pyleos/tests/reference_os.json, per-atom labels from clean_set.json; the
labels are checked to reproduce the reference element means before scoring.
Unlike full_bench.py this writes no file unless --out / --json is given.

Assignment methods are plugins: a .py file defining assign(case) -> list of
integer oxidation states in POSCAR atom order (optional META dict). The
default plugin pyleos/search/baseline.py is the shipped transport. Before a
plugin sees a case, outputs of other methods (EOSCAR greedy/flow columns,
greedy host fields of #WF) and the #WP block (invalid in archived data, see
patch/README.md) are removed, and the case path is anonymised. The run
records whether the plugin touched the electronegativity or the oxidation-
state bound tables.

Usage: python3 script/clean_bench.py [--split all|dev] [--method PLUGIN.py]
       [--json FILE] [--compare BASELINE.json] [--out FILE] [--greedy]
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
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "script"))

import pyleos.capacity as capacity                              # noqa: E402
import pyleos.transport as transport                            # noqa: E402
from full_bench import collect as collect_all, element_means    # noqa: E402
from pyleos.parsers import load_case                            # noqa: E402

print = functools.partial(print, flush=True)
CLEAN = json.loads((ROOT / "pyleos/tests/clean_set.json").read_text())
DEFAULT_METHOD = ROOT / "pyleos/search/baseline.py"


def collect():
    """full_bench.collect() with the data directories that clean_set.json
    overrides (reruns that fix a clean-set case without touching the 74-case
    benchmark data)."""
    over = {(b, n): ROOT / d for b, cases in CLEAN.get("overrides", {}).items()
            if not b.startswith("_") for n, d in cases.items()}
    for batch, name, d, ref in collect_all():
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


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--split", choices=["all", "dev"], default="all")
    ap.add_argument("--method", type=Path, default=DEFAULT_METHOD,
                    help="assignment plugin (.py defining assign(case))")
    ap.add_argument("--greedy", action="store_true",
                    help="score the in-situ greedy OS column of EOSCAR instead")
    ap.add_argument("--json", type=Path, help="write machine-readable results")
    ap.add_argument("--compare", type=Path,
                    help="results JSON of a reference run; report per-case regressions")
    ap.add_argument("--out", type=Path, help="write a per-case TSV")
    args = ap.parse_args()

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
    for batch, name, d, ref in collect():
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


if __name__ == "__main__":
    main()
