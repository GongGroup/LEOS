#!/usr/bin/env python3
"""Benchmark figure of the in-situ electronegativity assignment on clean38.

Counterpart of fig2_benchmark.py for the frozen plugin
pyleos/methods/frozen/r3-T-B-cap-repair.py ("r3-T-B", this work).  clean38 is
the development split of the clean benchmark: 38 cases, 793 atoms, 769 of
them with a per-atom reference (`python3 -m pyleos bench clean --split dev`).

Upper block - accuracy on clean38
(a) draw_accuracy   atoms correct, cases fully correct and element-mean pass
                    of five assignment methods.
(b) draw_cases      incorrect atoms of this work, the min-cost flow and the
                    per-WF greedy rule in every case that at least one of the
                    three does not get fully right.
(c) draw_elements   periodic-table map: atoms correct / labelled per element.
Lower block - where the polarity comes from
(d) draw_chi        in-situ chi = Z_eff / rho^2 of the elements that carry
                    lone pairs against the Allen electronegativity.
(e) draw_constants  atoms correct as each of the four constants is varied.
(f) draw_margins    winner / loser ratio of the deciding quantity for every
                    contested heteronuclear WF, by bond class and rule.

Run from this directory:
    python3 r3tb_fig2_benchmark.py              # uses the caches in _cache/r3tb/
    python3 r3tb_fig2_benchmark.py --recompute  # rebuild both caches first
                                                # (about 75 s on a laptop)
    python3 r3tb_fig2_benchmark.py --panel margins --scale 3 --out FILE.png
                                                # one panel at poster scale

Writes r3tb_fig2_benchmark.pdf / .png here, the two halves used on slides
(_cache/r3tb/fig2_top.png, panels a-c; _cache/r3tb/fig2_bottom.png, panels
d-f; 300 dpi) and the log _cache/r3tb/fig2_benchmark.run.log.

Nothing plotted is typed in by hand.  (a), (b): results/clean38_*.tsv.
(c), (d), (f): the frozen plugin re-run in process on the evaluator's own
case list, labels and scrub (cache _cache/r3tb/fig2_analysis.json; checked
case by case against results/clean38_r3-T-B.tsv).  (e): the frozen evaluator
run on copies of the plugin that differ from it in one constant (cache
_cache/r3tb/fig2_scan.json; the copies live in a temporary directory outside
the repository and the frozen file is never edited).  The Allen
electronegativities of pyleos.capacity are used as the abscissa of (d) and to
flag one row of (f); the assignment itself reads no table.

Every panel is a function draw_<name>(ax) that fetches its own data, so that
a poster can redraw single panels: all offsets, marker sizes and line widths
follow rcParams["font.size"] (7.5 pt = journal scale); `scaled(k)` is a
context manager that multiplies the journal style by k.

Needs numpy + matplotlib + the standard library, and the repository (pyleos).
"""
import argparse
import contextlib
import csv
import functools
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path

sys.dont_write_bytecode = True      # leave no __pycache__ in the repository

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedLocator, NullFormatter

import style
from style import C, MARKER, REPO, HERE, formula

sys.path.insert(0, str(REPO))

# ------------------------------------------------------------------ inputs
RESULTS = REPO / "results"
CACHE = HERE / "_cache" / "r3tb"
FROZEN = REPO / "pyleos/methods/frozen/r3-T-B-cap-repair.py"
ANALYSIS_JSON = CACHE / "fig2_analysis.json"
SCAN_JSON = CACHE / "fig2_scan.json"
LOG = CACHE / "fig2_benchmark.run.log"
NAME_OUT = "r3tb_fig2_benchmark"

# key, style key (colour + marker of the team brief), label, per-case table,
# how the table was produced (`python3 -m pyleos bench clean --split dev ...`)
METHODS = (
    ("r3tb", "aux1", "In-situ electronegativity (this work)",
     "clean38_r3-T-B.tsv", "--method pyleos/methods/frozen/r3-T-B-cap-repair.py"),
    ("flow", "flow", "Min-cost flow",
     "clean38_flow.tsv", "--method pyleos/methods/baseline.py"),
    ("greedy", "nearest", "Per-WF greedy rule",
     "clean38_greedy.tsv", "--greedy"),
    ("geo", "geo", r"Min-cost flow, $\beta$ = 0",
     "clean38_geo.tsv", "--method pyleos/methods/baseline_geo.py"),
    ("lewis", "aux2", "Lewis-first with Allen table",
     "clean38_r2-W-A.tsv",
     "--method pyleos/methods/frozen/r2-W-A-lewis-first-hardened.py"),
)
STYLE = {m[0]: m[1] for m in METHODS}
MLABEL = {m[0]: m[2] for m in METHODS}
THIS = "r3tb"
TRIO = ("r3tb", "flow", "greedy")           # methods of panel (b)

# Panel (a): (numerator, denominator, axis label)
METRICS = (("atoms_ok", "atoms_labelled", "Atoms correct (%)"),
           ("full", "cases", "Cases fully correct (%)"),
           ("elem_ok", "cases", "Element-mean pass (%)"))

# Panel (e): the four constants of the plugin.  `line` is the start of the
# source line that is replaced in a copy; `values` is the scan grid.
def _grid(lo, hi, step, extra=()):
    n = int(round((hi - lo) / step))
    return sorted({round(lo + i * step, 10) for i in range(n + 1)} | set(extra))


SCAN = {
    "PAIR_SCREEN": dict(values=_grid(0.0, 1.0, 0.025, (5.0 / 16.0,)),
                        label="Pair screening\n(electrons)", default_text="5/16",
                        xlim=(-0.05, 1.05), ticks=(0, 0.5, 1), minor=0.25),
    "PARITY": dict(values=_grid(0.0, 0.5, 0.0125),
                   label="Dipole-free limit\n(rms radii)", default_text="1/8",
                   xlim=(-0.025, 0.525), ticks=(0, 0.25, 0.5), minor=0.125),
    "REACH_SIGMA": dict(values=_grid(0.5, 2.0, 0.05),
                        label="WF extent\n(rms radii)", default_text="1",
                        xlim=(0.425, 2.075), ticks=(0.5, 1, 1.5, 2), minor=0.25),
    "FLOAT_EPS": dict(values=[float(f"1e{p}") for p in range(-12, -1)],
                      label="Numerical\ntolerance", default_text="10$^{-6}$",
                      xlim=(10 ** -12.7, 10 ** -1.3), ticks=(1e-10, 1e-4),
                      log=True, default_side="left"),
}
SCAN_WIDTH = (1.0, 1.0, 1.0, 1.0)      # relative widths of the sub-axes of (e)

# Display formulas of the case directory names (labels only, no data).
NAME = {
    "as4s4": "As4S4", "cro5": "CrO5", "h2s2": "H2S2", "h2s3": "H2S3",
    "h2s5": "H2S5", "k2s2o3": "K2S2O3", "kmgf3": "KMgF3", "mocl2": "MoCl2",
    "nan3": "NaN3", "rb2cucl4": "Rb2CuCl4", "re2": "Re2", "ref6": "ReF6",
    "sega": "GaSe", "srbe_oh_4": "SrBe(OH)4",
}

# Periodic table, periods 1-6 (the f row holds Ce-Lu); symbol -> (group, row)
_PERIODS = ("H He", "Li Be B C N O F Ne", "Na Mg Al Si P S Cl Ar",
            "K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr",
            "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe",
            "Cs Ba La Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn")
_F_ROW = "Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu"
F_ROW_Y = 6.3                       # row coordinate of the lanthanides


def _table():
    pos = {}
    for row, line in enumerate(_PERIODS):
        syms = line.split()
        if row == 0:
            groups = (1, 18)
        elif row in (1, 2):
            groups = (1, 2) + tuple(range(13, 19))
        else:
            groups = tuple(range(1, 19))
        for sym, g in zip(syms, groups):
            pos[sym] = (g, float(row))
    for i, sym in enumerate(_F_ROW.split()):
        pos[sym] = (4 + i, F_ROW_Y)
    return pos


TABLE = _table()


def block_of(sym):
    g, row = TABLE[sym]
    if row == F_ROW_Y:
        return "f"
    return "s" if g <= 2 else "d" if g <= 12 else "p"


# ---------------------------------------------------------------- utilities
def rel(path):
    try:
        return str(Path(path).relative_to(REPO))
    except ValueError:
        return str(path)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_rows(path):
    with Path(path).open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def upright(s):
    """Plain text as mathtext, so that it is set like the formulas."""
    return rf"$\mathrm{{{s}}}$"


def pretty(case):
    s = formula(NAME.get(case.lower(), case))
    return s if "$" in s else upright(s)


def ranks(v):
    """Average ranks (ties share the mean rank)."""
    v = np.asarray(v, float)
    order = np.argsort(v, kind="stable")
    r = np.empty(len(v))
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
            j += 1
        r[order[i:j + 1]] = (i + j) / 2.0
        i = j + 1
    return r


def spearman(x, y):
    rx, ry = ranks(x), ranks(y)
    return float(np.corrcoef(rx, ry)[0, 1])


# ------------------------------------------------------------ data: (a, b)
def clean38_keys():
    """(batch, case) of the development split, in the evaluator's order."""
    import pyleos.bench as cb
    held = {tuple(x) for x in cb.CLEAN["splits"]["test_regression"]}
    return [(b, n) for b, cases in cb.CLEAN["cases"].items() for n in cases
            if (b, n) not in held]


def read_bench():
    """Per-case rows and totals of the five methods (results/clean38_*.tsv)."""
    want = set(clean38_keys())
    out = {}
    for key, _sty, _lab, fname, _cmd in METHODS:
        rows = read_rows(RESULTS / fname)
        got = [(r["batch"], r["case"]) for r in rows]
        if set(got) != want or len(got) != len(want):
            raise ValueError(f"{fname}: cases differ from the development "
                             f"split of pyleos/sets/clean_set.json")
        for r in rows:
            r["atoms_ok"] = int(r["atoms_ok"])
            r["atoms_labelled"] = int(r["atoms_labelled"])
        out[key] = {
            "file": RESULTS / fname, "rows": {(r["batch"], r["case"]): r for r in rows},
            "summary": {
                "cases": len(rows),
                "atoms_ok": sum(r["atoms_ok"] for r in rows),
                "atoms_labelled": sum(r["atoms_labelled"] for r in rows),
                "full": sum(r["atoms_ok"] == r["atoms_labelled"]
                            and r["status"] != "ERROR" for r in rows),
                "elem_ok": sum(r["status"] == "OK" for r in rows)}}
    ref = out[THIS]["rows"]
    for key in out:
        for k, r in out[key]["rows"].items():
            if r["atoms_labelled"] != ref[k]["atoms_labelled"]:
                raise ValueError(f"{key}: {k} has another number of labelled atoms")
    return out


# --------------------------------------------------------- data: (c, d, f)
def _ends(plugin, case):
    """The two ends (winner, runner-up) of every contested WF that has a
    heteronuclear partner, with Z_eff and rho^2 of both ends.

    `_analyse` returns the geometry of every WF and, in `notes`, the
    elements, Z_eff and distances of the two ends, but not the atoms or the
    lone-pair radii.  They are rebuilt here with the plugin's own shell model
    (`_Shells.seen_from`, `_in_bond_region`, PAIR_SCREEN) in the order of
    `_analyse`, and every record is required to equal the plugin's note."""
    an = plugin._analyse(case)
    geo, els, zval, na = an["geo"], an["els"], an["zval"], an["na"]
    own = [[] for _ in range(na)]
    for g in geo:
        if g["own"]:
            a = g["near"]
            own[a].append((math.sqrt(g["d2"][a]), g["rad"], g["e"], g["pp"]))
    shells = [plugin._Shells(rows) for rows in own]
    out = []
    for g in geo:
        if not g["contested"]:
            continue
        vecs, d2 = g["vecs"], g["d2"]
        zeff, rho2, weak = [], [], []
        for a in range(na):
            inner, r2 = shells[a].seen_from(math.sqrt(d2[a]))
            z = zval[els[a]] - inner - plugin.PAIR_SCREEN
            zeff.append(z)
            rho2.append(r2)
            weak.append(d2[a] / z if z > 0.0 else math.inf)
        order = sorted(range(na), key=lambda a: (weak[a], d2[a], a))
        top = order[0]

        def partners_of(t):
            return [b for b in order if els[b] != els[t]
                    and plugin._in_bond_region(vecs[t], vecs[b])]

        partners = partners_of(top)
        if rho2[top] is not None:
            best = top
            for b in partners:
                if rho2[b] is not None and zeff[b] / rho2[b] > zeff[best] / rho2[best]:
                    best = b
            if best != top:
                top, partners = best, partners_of(best)
        if partners:
            b = partners[0]
            out.append({"top": top, "b": b, "zeff": (zeff[top], zeff[b]),
                        "rho2": (rho2[top], rho2[b]),
                        "d": (math.sqrt(d2[top]), math.sqrt(d2[b])),
                        "hosts": len(g["hosts"])})
    notes = an["notes"]
    if len(out) != len(notes):
        raise RuntimeError("ends: number of records differs from the plugin's notes")
    for o, n in zip(out, notes):
        both = o["rho2"][0] is not None and o["rho2"][1] is not None
        if not (els[o["top"]] == n["winner"] and els[o["b"]] == n["loser"]
                and o["zeff"] == n["zeff"] and o["d"] == n["d"]
                and ("chi" if both else "force") == n["rule"]):
            raise RuntimeError("ends: a record differs from the plugin's note")
        if both and n["chi_margin"] is not None:
            mine = (o["zeff"][0] / o["rho2"][0]) / (o["zeff"][1] / o["rho2"][1])
            if mine != n["chi_margin"]:
                raise RuntimeError("ends: chi ratio differs from the plugin's note")
        o["rule"] = n["rule"]
        o["margin"] = n["chi_margin"] if n["rule"] == "chi" else n["force_margin"]
    return an, out


def _check_pair_margins(plugin, case, elements, ends):
    """The records of one case, grouped by (winner, loser), must reproduce the
    plugin's public diagnostic pair_margins(case)."""
    mine = {}
    for o in ends:
        rec = mine.setdefault((elements[o["top"]], elements[o["b"]]),
                              {"n": 0, "by_chi": 0, "by_force": 0,
                               "min_chi_margin": None, "min_force_margin": None})
        rec["n"] += 1
        rec["by_" + o["rule"]] += 1
        key = f"min_{o['rule']}_margin"
        if o["margin"] is not None and (rec[key] is None or o["margin"] < rec[key]):
            rec[key] = o["margin"]
    if mine != plugin.pair_margins(case):
        raise RuntimeError("records differ from the plugin's pair_margins()")


def compute_analysis():
    """Re-run the frozen plugin on clean38 (evaluator's cases, labels, scrub)."""
    import pyleos.bench as cb
    from pyleos.parsers import load_case
    plugin = cb.load_plugin(FROZEN)
    spec = {(b, n): s for b, cases in cb.CLEAN["cases"].items()
            for n, s in cases.items()}
    want = set(clean38_keys())
    tsv = {(r["batch"], r["case"]): r
           for r in read_rows(RESULTS / "clean38_r3-T-B.tsv")}
    cases, wfs, chi_atoms = [], [], defaultdict(list)
    el_atoms, el_lab, el_ok, el_cases = Counter(), Counter(), Counter(), defaultdict(set)
    t_assign = 0.0
    for batch, name, d, ref in cb.collect_clean():
        if (batch, name) not in want:
            continue
        case = cb.scrub(load_case(d))
        st = case.structure
        labels = cb.site_labels(name, st, ref, spec[(batch, name)])
        t0 = time.time()
        os_vals = [int(v) for v in plugin.assign(case)]
        t_assign += time.time() - t0
        wrong = Counter((el, lab, o) for el, lab, o in zip(st.elements, labels, os_vals)
                        if lab is not None and o != lab)
        wrong_s = ",".join(f"{el}{lab:+d}->{o:+d}x{k}"
                           for (el, lab, o), k in sorted(wrong.items())) or "-"
        n_lab = sum(lab is not None for lab in labels)
        n_ok = n_lab - sum(wrong.values())
        row = tsv[(batch, name)]
        if (int(row["atoms_ok"]), int(row["atoms_labelled"]), row["wrong_sites"]) \
                != (n_ok, n_lab, wrong_s):
            raise RuntimeError(f"{batch}/{name}: re-run differs from "
                               f"results/clean38_r3-T-B.tsv")
        for el, lab, o in zip(st.elements, labels, os_vals):
            el_atoms[el] += 1
            el_cases[el].add(f"{batch}/{name}")
            if lab is not None:
                el_lab[el] += 1
                el_ok[el] += (o == lab)
        an, ends = _ends(plugin, case)
        _check_pair_margins(plugin, case, st.elements, ends)
        per_atom = defaultdict(list)
        for o in ends:
            chi = [None, None]
            for k, a in enumerate((o["top"], o["b"])):
                if o["rho2"][k] is not None:
                    chi[k] = o["zeff"][k] / o["rho2"][k]
                    per_atom[a].append(chi[k])
            m = o["margin"]
            wfs.append({"batch": batch, "case": name,
                        "winner": st.elements[o["top"]], "loser": st.elements[o["b"]],
                        "rule": o["rule"], "hosts": o["hosts"],
                        "margin": "inf" if m == math.inf else m,
                        "chi": chi, "zeff": list(o["zeff"]), "d": list(o["d"])})
        for a, v in sorted(per_atom.items()):
            chi_atoms[st.elements[a]].append(
                {"case": f"{batch}/{name}", "atom": a + 1, "n": len(v),
                 "median": float(np.median(v)), "min": min(v), "max": max(v)})
        n_wf = sum(len(ch.centres) for ch in case.channels)
        cases.append({"batch": batch, "case": name, "dir": rel(d),
                      "atoms": len(st.elements), "labelled": n_lab, "ok": n_ok,
                      "wrong": wrong_s, "wfs": n_wf,
                      "contested": sum(g["contested"] for g in an["geo"]),
                      "hetero": len(ends),
                      "chi": sum(o["rule"] == "chi" for o in ends),
                      "force": sum(o["rule"] == "force" for o in ends)})
    if {(c["batch"], c["case"]) for c in cases} != want:
        raise RuntimeError("clean38 cases without data")
    return {
        "produced": datetime.now().isoformat(timespec="seconds"),
        "plugin": rel(FROZEN), "plugin_sha256": sha256(FROZEN),
        "constants": {k: getattr(plugin, k) for k in SCAN},
        "assign_seconds": round(t_assign, 2),
        "cases": cases,
        "elements": {el: {"atoms": el_atoms[el], "labelled": el_lab[el],
                          "ok": el_ok[el], "cases": sorted(el_cases[el])}
                     for el in sorted(el_atoms)},
        "wfs": wfs,
        "chi_atoms": dict(chi_atoms),
    }


def analysis(recompute=False):
    if ANALYSIS_JSON.is_file() and not recompute:
        res = json.loads(ANALYSIS_JSON.read_text())
        if res.get("plugin_sha256") == sha256(FROZEN):
            return res
        print(f"{rel(ANALYSIS_JSON)} belongs to another plugin file: recomputing",
              flush=True)
    print("running the frozen plugin on clean38 ...", flush=True)
    res = compute_analysis()
    CACHE.mkdir(parents=True, exist_ok=True)
    ANALYSIS_JSON.write_text(json.dumps(res, indent=1) + "\n")
    return res


# --------------------------------------------------------------- data: (e)
def _scan_job(job):
    name, value, src, folder = job
    new, n = re.subn(rf"(?m)^{name} = [^#\n]*", f"{name} = {value!r}   ", src,
                     count=1)
    if n != 1:
        raise RuntimeError(f"{name}: assignment line not found in {FROZEN.name}")
    copy = folder / f"r3-T-B_{name}_{value!r}.py"
    copy.write_text(new)
    out = copy.with_suffix(".json")
    cmd = [sys.executable, "-m", "pyleos", "bench", "clean", "--split", "dev",
           "--method", str(copy), "--json", str(out)]
    r = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True,
                       env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
    if r.returncode != 0:
        raise RuntimeError(f"{name} = {value!r}: evaluator failed\n{r.stderr[-400:]}")
    res = json.loads(out.read_text())
    if res["split"] != "dev":
        raise RuntimeError("scan: not the development split")
    return name, value, res


def compute_scan(scan_dir=None):
    """Frozen evaluator on copies of the plugin with one constant changed.

    The copies are written to scan_dir (default: a temporary directory that
    is removed afterwards), never into the repository; each differs from the
    frozen file in the one line that assigns the constant."""
    import pyleos.bench as cb
    src = FROZEN.read_text()
    plugin = cb.load_plugin(FROZEN)
    defaults = {k: getattr(plugin, k) for k in SCAN}
    folder = Path(scan_dir) if scan_dir else Path(tempfile.mkdtemp(prefix="r3tb_scan_"))
    folder.mkdir(parents=True, exist_ok=True)
    if REPO in folder.resolve().parents or folder.resolve() == REPO:
        raise RuntimeError("the scan copies must stay outside the repository")
    jobs = []
    for name, spec in SCAN.items():
        if defaults[name] not in spec["values"]:
            raise RuntimeError(f"{name}: the default {defaults[name]} is not on the grid")
        jobs += [(name, v, src, folder) for v in spec["values"]]
    print(f"scan: {len(jobs)} evaluator runs on copies in {folder} ...", flush=True)
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=min(6, os.cpu_count() or 2)) as ex:
        done = list(ex.map(_scan_job, jobs))
    print(f"scan: done in {time.time() - t0:.0f} s", flush=True)
    if not scan_dir:
        shutil.rmtree(folder)

    tsv = {(r["batch"], r["case"]): (int(r["atoms_ok"]), r["status"])
           for r in read_rows(RESULTS / "clean38_r3-T-B.tsv")}
    out = {}
    for name in SCAN:
        runs = [(v, res) for n, v, res in done if n == name]
        base = next(res for v, res in runs if v == defaults[name])
        base_rows = {(c["batch"], c["case"]): (c["atoms_ok"], c["status"])
                     for c in base["cases"]}
        if base_rows != tsv:
            raise RuntimeError(f"{name}: the copy with the default value does not "
                               f"reproduce results/clean38_r3-T-B.tsv")
        points = []
        for v, res in runs:
            s = res["summary"]
            changed = [{"case": f"{c['batch']}/{c['case']}", "atoms_ok": c["atoms_ok"],
                        "atoms_labelled": c["atoms_labelled"], "status": c["status"],
                        "error": c["error"]}
                       for c in res["cases"]
                       if (c["atoms_ok"], c["status"]) != base_rows[(c["batch"], c["case"])]]
            points.append({"value": v, "atoms_ok": s["atoms_ok"],
                           "atoms_labelled": s["atoms_labelled"],
                           "full": s["all_atoms_correct_cases"],
                           "elem_ok": s["elem_ok"], "cases": s["cases"],
                           "changed": changed})
        out[name] = {"default": defaults[name], "points": points}
    return {
        "produced": datetime.now().isoformat(timespec="seconds"),
        "plugin": rel(FROZEN), "plugin_sha256": sha256(FROZEN),
        "command": "python3 -m pyleos bench clean --split dev --method <copy> "
                   "--json <copy>.json",
        "copy": "frozen plugin with the line '<NAME> = ...' replaced by "
                "'<NAME> = <value>'; one constant per copy, the other three at "
                "their defaults; copies written outside the repository",
        "copies_dir": str(folder) if scan_dir else "temporary directory (removed)",
        "scan": out,
    }


def scan(recompute=False, scan_dir=None):
    if SCAN_JSON.is_file() and not recompute:
        res = json.loads(SCAN_JSON.read_text())
        grid_ok = all([p["value"] for p in res["scan"].get(k, {}).get("points", [])]
                      == SCAN[k]["values"] for k in SCAN)
        if res.get("plugin_sha256") == sha256(FROZEN) and grid_ok:
            return res
        print(f"{rel(SCAN_JSON)} belongs to another plugin file or grid: "
              f"recomputing", flush=True)
    res = compute_scan(scan_dir)
    CACHE.mkdir(parents=True, exist_ok=True)
    SCAN_JSON.write_text(json.dumps(res, indent=1) + "\n")
    return res


_OPTS = {"recompute": False, "scan_dir": None}


@functools.lru_cache(maxsize=None)
def data():
    """Everything the panels plot (read once)."""
    return {"bench": read_bench(),
            "analysis": analysis(_OPTS["recompute"]),
            "scan": scan(_OPTS["recompute"], _OPTS["scan_dir"])}


def plateau(points, default):
    """Largest interval of scan values around the default over which every
    case has the result of the default run."""
    vals = [p["value"] for p in points]
    i0 = vals.index(default)
    lo = hi = i0
    while lo > 0 and not points[lo - 1]["changed"]:
        lo -= 1
    while hi < len(vals) - 1 and not points[hi + 1]["changed"]:
        hi += 1
    return vals[lo], vals[hi]


def chi_elements():
    """Per element: in-situ chi (median, quartiles and range over its atoms;
    the chi of an atom is the median over the contested WFs it is an end of)
    and the Allen electronegativity."""
    from pyleos.capacity import ENEG_ALLEN
    out = []
    for el, atoms in data()["analysis"]["chi_atoms"].items():
        v = np.array([a["median"] for a in atoms])
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        out.append({"el": el, "allen": float(ENEG_ALLEN[el]), "chi": float(med),
                    "q1": float(q1), "q3": float(q3),
                    "lo": float(v.min()), "hi": float(v.max()),
                    "atoms": len(atoms), "ends": sum(a["n"] for a in atoms),
                    "end_lo": min(a["min"] for a in atoms),
                    "end_hi": max(a["max"] for a in atoms),
                    "cases": sorted({a["case"] for a in atoms}),
                    "block": block_of(el)})
    return sorted(out, key=lambda r: r["allen"])


def margin_classes():
    """Bond classes of panel (f): contested WFs with a heteronuclear partner,
    grouped by (winner, loser) and by whether one atom or several atoms of the
    winning element hold the WF."""
    from pyleos.capacity import ENEG_ALLEN
    groups = defaultdict(list)
    for w in data()["analysis"]["wfs"]:
        groups[(w["winner"], w["loser"], w["hosts"] == 1)].append(w)
    out = []
    for (win, los, single), ws in groups.items():
        m = [math.inf if w["margin"] == "inf" else w["margin"] for w in ws]
        hosts = sorted({w["hosts"] for w in ws})
        out.append({"winner": win, "loser": los, "single": single, "hosts": hosts,
                    "n": len(ws), "median": float(np.median(m)),
                    "force": [x for x, w in zip(m, ws) if w["rule"] == "force"],
                    "chi": [x for x, w in zip(m, ws) if w["rule"] == "chi"],
                    "allen": ENEG_ALLEN[win] - ENEG_ALLEN[los],
                    "cases": sorted({w["case"] for w in ws})})
    key = lambda c: (c["median"], min(c["force"] + c["chi"]), c["winner"], c["loser"])
    return (sorted((c for c in out if c["single"]), key=key),
            sorted((c for c in out if not c["single"]), key=key))


# ----------------------------------------------------------------- drawing
BASE_FONT = 7.5         # font.size of the journal style


def K():
    """Scale of the current style relative to the journal style."""
    return plt.rcParams["font.size"] / BASE_FONT


_SCALED_KEYS = ("font.size", "axes.labelsize", "axes.titlesize", "xtick.labelsize",
                "ytick.labelsize", "legend.fontsize", "axes.linewidth",
                "xtick.major.size", "ytick.major.size", "xtick.major.width",
                "ytick.major.width", "xtick.minor.size", "ytick.minor.size",
                "xtick.minor.width", "ytick.minor.width", "lines.linewidth",
                "lines.markersize", "patch.linewidth", "hatch.linewidth",
                "xtick.major.pad", "ytick.major.pad", "axes.labelpad")


@contextlib.contextmanager
def scaled(k):
    """Journal style multiplied by k (fonts, lines, ticks) for poster panels."""
    with mpl.rc_context({key: plt.rcParams[key] * k for key in _SCALED_KEYS}):
        yield


def cm_axes(fig, left, top, width, height):
    """Axes placed in cm from the left / top edges of the figure."""
    fw, fh = fig.get_size_inches() / style.CM
    return fig.add_axes([left / fw, 1 - (top + height) / fh,
                         width / fw, height / fh])


def ax_cm(ax):
    """(width, height) of an axes in cm."""
    fw, fh = ax.figure.get_size_inches() / style.CM
    p = ax.get_position()
    return p.width * fw, p.height * fh


def baseline_shift(s, size):
    """Points to lower a baseline so that capitals are centred on the anchor
    (see fig2_benchmark.py)."""
    return -0.358 * size - (0.45 * K() if "$" in s else 0.0)


def row_text(ax, x, y, s, ha="left", transform=None, **kw):
    """Text on a row, baseline-aligned."""
    size = kw.get("fontsize", plt.rcParams["font.size"])
    tr = mtransforms.offset_copy(transform or ax.transData, fig=ax.figure,
                                 y=baseline_shift(s, size), units="points")
    return ax.text(x, y, s, ha=ha, va="baseline", transform=tr, **kw)


def glyph(ax, x, y, method, ms=3.6, transform=None, **kw):
    """Marker of a method (colour + shape of the team brief)."""
    sty = STYLE[method]
    opts = dict(ls="none", marker=MARKER[sty], ms=ms * K(), color=C[sty],
                mec=C[sty], clip_on=False)
    opts.update(kw)
    return ax.plot([x], [y], transform=transform or ax.transData, **opts)


def host(ax):
    """Turn an axes into an invisible frame for sub-axes and free text."""
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    return ax_cm(ax)


# ---------------------------------------------------------------- panel (a)
def draw_accuracy(ax):
    """(a) Three accuracy measures of five methods on clean38 (ax is the
    frame of the whole strip: method names, three bar blocks, counts)."""
    k = K()
    bench = data()["bench"]
    W, H = host(ax)
    name_w, count_w, gap = 4.85 * k, 1.28 * k, 0.42 * k
    bar_w = (W - name_w - 3 * count_w - 2 * gap) / 3
    n = len(METHODS)
    ylim = (n - 0.42, -0.58)
    ax.set_ylim(ylim)
    yt = ax.get_yaxis_transform()
    for i, (key, *_r) in enumerate(METHODS):
        glyph(ax, 0.10 * k / W, i, key, ms=4.2, transform=yt)
        row_text(ax, 0.32 * k / W, i, MLABEL[key], transform=yt)
    subs = []
    for j, (num, den, label) in enumerate(METRICS):
        x0 = name_w + j * (bar_w + count_w + gap)
        sub = ax.inset_axes([x0 / W, 0, bar_w / W, 1])
        sub.set_xlim(0, 100)
        sub.set_ylim(ylim)
        h = 0.20 * k * abs(ylim[1] - ylim[0]) / H
        for i, (key, *_r) in enumerate(METHODS):
            s = bench[key]["summary"]
            sub.barh(i, 100.0 * s[num] / s[den], height=h, color=C[STYLE[key]],
                     edgecolor=C[STYLE[key]], lw=0.5 * k)
            row_text(sub, 1 + 0.12 * k / bar_w, i, f"{s[num]}/{s[den]}",
                     transform=sub.get_yaxis_transform())
        sub.set_yticks([])
        sub.set_xticks([0, 50, 100])
        sub.set_xticks([25, 75], minor=True)
        sub.set_xlabel(label)
        subs.append(sub)
    return subs


# ---------------------------------------------------------------- panel (b)
def case_rows():
    """Cases in which at least one of the three methods has incorrect atoms."""
    bench = data()["bench"]
    rows = []
    for key, r in bench[THIS]["rows"].items():
        n = r["atoms_labelled"]
        wrong = {m: n - bench[m]["rows"][key]["atoms_ok"] for m in TRIO}
        if any(wrong.values()):
            rows.append({"key": key, "n": n, "wrong": wrong})
    order = lambda r: (-r["wrong"]["flow"] / r["n"], -r["wrong"]["greedy"] / r["n"],
                       r["key"][1].lower())
    clean = sorted((r for r in rows if not r["wrong"][THIS]), key=order)
    rest = sorted((r for r in rows if r["wrong"][THIS]), key=order)
    return clean, rest, len(bench[THIS]["rows"]) - len(rows)


CASE_COLS = (0.40, 0.92, 1.44, 2.02)      # cm right of the axes: three counts, n
CASE_GAP = 0.7                             # rows between the two groups


def draw_cases(ax):
    """(b) Incorrect atoms (% of the labelled atoms) of three methods, case
    by case; counts in the columns on the right."""
    k = K()
    clean, rest, _ = case_rows()
    w_cm, h_cm = ax_cm(ax)
    yt = ax.get_yaxis_transform()
    size = {"greedy": 5.4, "flow": 4.3, THIS: 3.5}     # larger below smaller
    y, y_of = 0.0, []
    for block in (clean, rest):
        for r in block:
            frac = {m: 100.0 * r["wrong"][m] / r["n"] for m in TRIO}
            ax.plot([min(frac.values()), max(frac.values())], [y, y],
                    color=C["muted"], lw=0.8 * k, zorder=1, solid_capstyle="butt")
            for z, m in enumerate(("greedy", "flow", THIS)):
                sty = STYLE[m]
                ax.plot(frac[m], y, ls="none", marker=MARKER[sty], ms=size[m] * k,
                        mfc=C[sty], mec="white", mew=0.3 * k, zorder=2 + z,
                        clip_on=False)
            row_text(ax, -0.16 * k / w_cm, y, pretty(r["key"][1]), ha="right",
                     transform=yt)
            for x, m in zip(CASE_COLS, TRIO):
                row_text(ax, 1 + x * k / w_cm, y, str(r["wrong"][m]), ha="center",
                         transform=yt)
            row_text(ax, 1 + CASE_COLS[3] * k / w_cm, y, str(r["n"]), ha="center",
                     transform=yt)
            y_of.append(y)
            y += 1
        if block is clean and rest:
            ax.axhline(y - 0.5 + CASE_GAP / 2, color=C["muted"], lw=0.5 * k,
                       ls=(0, (4, 2.5)), zorder=0)
            y += CASE_GAP
    ax.set_ylim(y - 1 + 0.75, -0.75)
    ax.set_xlim(-6, 104)
    ax.set_yticks([])
    ax.set_xticks([0, 50, 100])
    ax.set_xticks([25, 75], minor=True)
    ax.set_xlabel("Incorrect atoms (%)")
    # column heads: the three method markers and the number of labelled atoms
    y_head = 1 + 0.22 * k / h_cm
    for x, m in zip(CASE_COLS, TRIO):
        glyph(ax, 1 + x * k / w_cm, y_head, m, ms=4.0, transform=ax.transAxes)
    tr = mtransforms.offset_copy(ax.transAxes, fig=ax.figure,
                                 y=baseline_shift("of", plt.rcParams["font.size"]),
                                 units="points")
    ax.text(1 + CASE_COLS[3] * k / w_cm, y_head, "of", ha="center", va="baseline",
            transform=tr)
    return clean, rest


# ---------------------------------------------------------------- panel (c)
TILE = 0.90                 # tile size in cell units (the rest is the gutter)


def draw_elements(ax):
    """(c) Periodic-table map of clean38: each tile is filled from the bottom
    in proportion to the labelled atoms of the element that are correct; the
    number is the count of labelled atoms."""
    k = K()
    els = data()["analysis"]["elements"]
    fs = plt.rcParams["font.size"]
    col = C[STYLE[THIS]]
    ax.set_xlim(0.5, 18.5)
    ax.set_ylim(F_ROW_Y + 0.5, -0.5)
    ax.set_axis_off()
    half = TILE / 2

    def tile(x, y, frac=None, w=TILE, h=TILE):
        if frac is None:                    # element not in clean38
            ax.add_patch(Rectangle((x - w / 2, y - h / 2), w, h, fc="none",
                                   ec=C["light"], lw=0.5 * k))
            return
        ax.add_patch(Rectangle((x - w / 2, y - h / 2), w, h, fc="white", ec="none"))
        ax.add_patch(Rectangle((x - w / 2, y + h / 2 - h * frac), w, h * frac,
                               fc=col, ec="none", alpha=0.62))
        ax.add_patch(Rectangle((x - w / 2, y - h / 2), w, h, fc="none", ec=col,
                               lw=0.6 * k))

    for sym, (g, row) in TABLE.items():
        if sym not in els:
            tile(g, row)
            continue
        e = els[sym]
        tile(g, row, e["ok"] / e["labelled"])
        ax.text(g, row - 0.08, sym, ha="center", va="baseline", fontsize=fs * 0.95)
        ax.text(g, row + 0.36, str(e["labelled"]), ha="center", va="baseline",
                fontsize=fs * 0.82)
    missing = set(els) - set(TABLE)
    if missing:
        raise RuntimeError(f"elements without a table position: {missing}")

    # key, in the empty region above the d block
    full = [s for s, e in els.items() if e["ok"] == e["labelled"]]
    part = sorted((s for s in els if s not in full),
                  key=lambda s: -els[s]["labelled"])
    x0, th, tw, pitch = 2.95, 0.56, 0.56, 0.72
    lines = [(1.0, f"every labelled atom correct ({len(full)} elements)"),
             (0.6, ", ".join(f"{s} {els[s]['ok']}/{els[s]['labelled']}"
                             for s in part) + " correct"),
             (None, "element not in clean38")]
    for i, (frac, text) in enumerate(lines):
        y = -0.12 + pitch * i
        tile(x0, y, frac, w=tw, h=th)
        ax.text(x0 + 0.52, y + 0.02, text, ha="left", va="center")
    ax.text(x0 - tw / 2, -0.12 + pitch * len(lines) + 0.02,
            "numbers: labelled atoms", ha="left", va="center")
    return full, part


# ---------------------------------------------------------------- panel (d)
CHI_PLACE = {            # label offsets in points (x, y, ha, va)
    "As": (-6, 0, "right", "center"), "I": (-1, 6.5, "center", "bottom"),
    "Se": (4.5, -3.5, "left", "top"), "S": (5.5, 3.5, "left", "bottom"),
    "Cl": (6, 0, "left", "center"), "O": (6, 0, "left", "center"),
    "F": (-6, 0, "right", "center"), "Mo": (6, 0, "left", "center"),
}


def draw_chi(ax):
    """(d) In-situ chi = Z_eff / rho^2 (median over the atoms that carry lone
    pairs at a contested WF; bar: interquartile range over these atoms)
    against the Allen electronegativity."""
    k = K()
    rows = chi_elements()
    col = C[STYLE[THIS]]
    mk = MARKER[STYLE[THIS]]
    for r in rows:
        p = r["block"] == "p"
        ax.plot([r["allen"]] * 2, [r["q1"], r["q3"]], color=col, lw=1.2 * k,
                solid_capstyle="butt", zorder=2)
        ax.plot(r["allen"], r["chi"], ls="none", marker=mk, ms=5.0 * k,
                mfc=col if p else "white", mec=col if not p else "white",
                mew=(0.9 if not p else 0.4) * k, zorder=3)
        dx, dy, ha, va = CHI_PLACE.get(r["el"], (6, 0, "left", "center"))
        ax.annotate(r["el"], (r["allen"], r["chi"]), xytext=(dx * k, dy * k),
                    textcoords="offset points", ha=ha, va=va, zorder=4)
    ax.set_yscale("log")
    lo = min(r["q1"] for r in rows)
    hi = max(r["q3"] for r in rows)
    ax.set_ylim(lo / 1.3, hi * 1.95)           # head room for the legend
    ax.set_xlim(1.2, 4.45)
    ax.set_yticks([t for t in (0.5, 1, 2, 5, 10, 20, 50)
                   if lo / 1.3 <= t <= hi * 1.95])
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _p: f"{v:g}"))
    ax.yaxis.set_minor_locator(FixedLocator(
        [m * 10 ** e for e in (-1, 0, 1) for m in range(1, 10)]))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xticks([1.5, 2, 2.5, 3, 3.5, 4])
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.25))
    ax.set_xlabel("Allen electronegativity")
    ax.set_ylabel(r"In-situ $\chi$ = $Z_{\mathrm{eff}}$/$\rho^{2}$ (e Å$^{-2}$)")
    pb = [r for r in rows if r["block"] == "p"]
    rs_p = spearman([r["allen"] for r in pb], [r["chi"] for r in pb])
    rs_all = spearman([r["allen"] for r in rows], [r["chi"] for r in rows])
    handles = [Line2D([], [], ls="none", marker=mk, ms=5.0 * k, mfc=col,
                      mec="white", mew=0.4 * k),
               Line2D([], [], ls="none", marker=mk, ms=5.0 * k, mfc="white",
                      mec=col, mew=0.9 * k)]
    labels = [rf"p-block elements: $r_{{\mathrm{{s}}}}$ = {rs_p:.2f} "
              f"($n$ = {len(pb)})",
              rf"d-block element; all elements: $r_{{\mathrm{{s}}}}$ = "
              f"{rs_all:.2f} ($n$ = {len(rows)})"]
    ax.legend(handles, labels, loc="upper left", handletextpad=0.3,
              labelspacing=0.45, borderaxespad=0.5)
    return rows, rs_p, rs_all


# ---------------------------------------------------------------- panel (e)
def draw_constants(ax):
    """(e) Atoms correct on clean38 against each of the four constants (ax
    is the frame of the four sub-axes); the other three stay at the default."""
    k = K()
    sc = data()["scan"]["scan"]
    W, H = host(ax)
    gap = 0.32 * k
    unit = (W - gap * (len(SCAN) - 1)) / sum(SCAN_WIDTH)
    col, mk = C[STYLE[THIS]], MARKER[STYLE[THIS]]
    n_all = {p["atoms_labelled"] for s in sc.values() for p in s["points"]}
    assert len(n_all) == 1
    n_all = n_all.pop()
    y_min = min(p["atoms_ok"] for s in sc.values() for p in s["points"])
    y_lo = 10 * math.floor((y_min - 4) / 10)
    x0, subs, info = 0.0, [], {}
    for i, (name, spec) in enumerate(SCAN.items()):
        w = unit * SCAN_WIDTH[i]
        sub = ax.inset_axes([x0 / W, 0, w / W, 1])
        x0 += w + gap
        pts = sc[name]["points"]
        default = sc[name]["default"]
        xs = [p["value"] for p in pts]
        ys = [p["atoms_ok"] for p in pts]
        at0 = ys[xs.index(default)]
        lo, hi = plateau(pts, default)
        if spec.get("log"):
            sub.set_xscale("log")
        sub.axvspan(lo, hi, color=C["light"], alpha=0.55, lw=0, zorder=0)
        sub.axvline(default, color=C["muted"], lw=0.5 * k, ls=(0, (4, 3)), zorder=1)
        sub.plot(xs, ys, color=col, lw=1.0 * k, zorder=3)
        sub.plot(xs, ys, ls="none", marker="o", ms=1.7 * k, mfc=col, mec=col,
                 zorder=3)
        sub.plot(default, at0, ls="none", marker=mk, ms=4.6 * k, mfc=col,
                 mec="white", mew=0.4 * k, zorder=4)
        sub.set_xlim(spec["xlim"])
        sub.set_ylim(y_lo, n_all)
        sub.set_xticks(list(spec["ticks"]))
        if spec.get("log"):
            sub.set_xticklabels([f"10$^{{{int(round(math.log10(t)))}}}$"
                                 for t in spec["ticks"]])
            sub.xaxis.set_minor_locator(FixedLocator(xs))
            sub.xaxis.set_minor_formatter(NullFormatter())
        else:
            sub.set_xticklabels([f"{t:g}" for t in spec["ticks"]])
            sub.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(spec["minor"]))
        sub.set_yticks(range(10 * math.ceil(y_lo / 20) * 2, n_all + 1, 20))
        sub.set_yticks(range(y_lo, n_all + 1, 10), minor=True)
        sub.set_xlabel(spec["label"], linespacing=1.15,
                       fontsize=plt.rcParams["font.size"])
        sub.xaxis.set_label_coords(0.5, -0.50 * k / H)
        if i:
            sub.tick_params(labelleft=False)
        else:
            sub.set_ylabel(f"Atoms correct\n(of {n_all})", linespacing=1.15)
        left = spec.get("default_side") == "left"
        sub.annotate("default\n" + spec["default_text"], (default, 0),
                     xycoords=sub.get_xaxis_transform(),
                     xytext=((-2.2 if left else 2.2) * k, 2.5 * k),
                     textcoords="offset points", ha="right" if left else "left",
                     va="bottom", linespacing=1.1,
                     fontsize=plt.rcParams["xtick.labelsize"])
        subs.append(sub)
        info[name] = dict(default=default, at_default=at0, plateau=(lo, hi),
                          lowest=min(ys), highest=max(ys))
    return subs, info


# ---------------------------------------------------------------- panel (f)
INF_X = 300.0                 # abscissa of the slot "loser fully screened"
INF_LO, INF_HI = 170.0, 530.0
F_COLS = (0.40, 0.95)         # cm right of the axes: counts by force, by chi


def class_label(c):
    if c["single"]:
        return upright(f"{c['winner']}–{c['loser']}")
    h = c["hosts"]
    sub = str(h[0]) if len(h) == 1 else f"{h[0]}–{h[-1]}"
    return rf"$\mathrm{{{c['winner']}}}_{{{sub}}}\mathrm{{–{c['loser']}}}$"


def draw_margins(ax):
    """(f) Winner / loser ratio of the deciding quantity (force Z_eff/d^2 or
    chi = Z_eff/rho^2) of every contested WF with a heteronuclear partner;
    rows are bond classes (winner first)."""
    k = K()
    single, shared = margin_classes()
    w_cm, h_cm = ax_cm(ax)
    yt = ax.get_yaxis_transform()
    col, mk = C[STYLE[THIS]], MARKER[STYLE[THIS]]
    force_kw = dict(ls="none", marker="o", ms=2.7 * k, mfc=col, mec="white",
                    mew=0.25 * k)
    chi_kw = dict(ls="none", marker=mk, ms=4.3 * k, mfc="white", mec=C["ink"],
                  mew=0.7 * k)
    y = 0.0
    blocks = (("One atom takes the pair", single),
              ("Like atoms share the pair", shared))

    def counts(y, n_f, n_c):
        for x, n in zip(F_COLS, (n_f, n_c)):
            row_text(ax, 1 + x * k / w_cm, y, str(n) if n else "–", ha="center",
                     transform=yt, color=C["ink"] if n else C["muted"])

    for title, rows in blocks:
        row_text(ax, 0.10 * k / w_cm, y, title, transform=yt)
        counts(y, sum(len(c["force"]) for c in rows),
               sum(len(c["chi"]) for c in rows))
        # over the rows of this block: the slot of the infinite ratios (loser
        # fully screened) and the line of equal pull (ratio 1)
        ax.add_patch(Rectangle((INF_LO, y + 0.5), INF_HI - INF_LO, len(rows),
                               fc=C["light"], ec="none", alpha=0.55, zorder=0))
        ax.plot([1.0, 1.0], [y + 0.5, y + 0.5 + len(rows)], color=C["ink"],
                lw=0.5 * k, ls=(0, (4, 2.5)), zorder=0)
        for c in rows:
            y += 1
            for vals, kw, z in ((c["force"], force_kw, 2), (c["chi"], chi_kw, 3)):
                if vals:
                    ax.plot([INF_X if math.isinf(v) else v for v in vals],
                            [y] * len(vals), zorder=z, **kw)
            row_text(ax, -0.14 * k / w_cm, y, class_label(c), ha="right",
                     transform=yt)
            counts(y, len(c["force"]), len(c["chi"]))
        y += 1.45
    n_force = sum(len(c["force"]) for c in single + shared)
    n_chi = sum(len(c["chi"]) for c in single + shared)
    y -= 0.25
    row_text(ax, -0.14 * k / w_cm, y, "All", ha="right", transform=yt)
    counts(y, n_force, n_chi)
    ax.set_ylim(y + 0.75, -0.85)
    ax.set_xscale("log")
    ax.set_xlim(0.8, INF_HI)
    ax.set_xticks([1, 10, 100, INF_X])
    ax.set_xticklabels(["1", "10", "100", "∞"])
    ax.xaxis.set_minor_locator(FixedLocator(
        [m * 10 ** e for e in (0, 1) for m in range(2, 10)]))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_yticks([])
    ax.set_xlabel("Ratio winner / loser")
    # legend above the axes; the same markers head the count columns
    big = ({**force_kw, "ms": 3.6 * k}, {**chi_kw, "ms": 4.6 * k})
    labels = [rf"force, $Z_{{\mathrm{{eff}}}}$/$d^{{2}}$",
              rf"$\chi$ = $Z_{{\mathrm{{eff}}}}$/$\rho^{{2}}$"]
    ax.legend([Line2D([], [], **kw) for kw in big], labels, loc="lower left",
              ncol=2, bbox_to_anchor=(-1.55 * k / w_cm, 1), borderaxespad=0.35,
              borderpad=0, columnspacing=1.2, handletextpad=0.25)
    y_head = 1 + 0.22 * k / h_cm
    for x, kw in zip(F_COLS, big):
        ax.plot([1 + x * k / w_cm], [y_head], transform=ax.transAxes,
                clip_on=False, **kw)
    return single, shared, n_force, n_chi


# ------------------------------------------------------------------ layout
LEFT, RIGHT = 0.62, 16.90                 # cm
H_A, H_X, BAND = 2.00, 0.85, 0.55
TOP_A = 0.52
N_CASE_ROWS = 14 + CASE_GAP + 1.5         # rows + gap + margins of panel (b)
H_B = N_CASE_ROWS * 0.30
TOP_B = TOP_A + H_A + H_X + BAND
H_C = (F_ROW_Y + 1.0) * 0.66
H_TOP = TOP_B + H_B + H_X + 0.10          # height of the upper block
H_F = 9.20
H_E = 2.15
H_D = H_F - H_E - 0.92 - BAND
TOP_2 = 0.62                              # axes top inside the lower block
H_BOTTOM = TOP_2 + H_F + 1.22
LABEL_X = (0.45, 10.15)                   # right edge of the panel labels


def panel_letter(ax, letter, x_cm, dy_cm=0.10):
    w, h = ax_cm(ax)
    fw = ax.figure.get_size_inches()[0] / style.CM
    style.panel_label(ax, letter, dx=(x_cm - ax.get_position().x0 * fw) / w,
                      dy=1 + dy_cm / h)


def build_top(fig, y0=0.0):
    ax_a = cm_axes(fig, LEFT, y0 + TOP_A, RIGHT - LEFT, H_A)
    ax_b = cm_axes(fig, 2.40, y0 + TOP_B, 2.60, H_B)
    ax_c = cm_axes(fig, 7.85, y0 + TOP_B + 0.05, RIGHT - 7.85, H_C)
    draw_accuracy(ax_a)
    draw_cases(ax_b)
    draw_elements(ax_c)
    panel_letter(ax_a, "a", LABEL_X[0])
    panel_letter(ax_b, "b", LABEL_X[0], dy_cm=0.10)
    panel_letter(ax_c, "c", 7.72, dy_cm=0.15)
    return ax_a, ax_b, ax_c


def build_bottom(fig, y0=0.0):
    left, right_l = 1.55, 9.25
    top = y0 + TOP_2
    ax_d = cm_axes(fig, left, top, right_l - left, H_D)
    ax_e = cm_axes(fig, left, top + H_F - H_E, right_l - left, H_E)
    ax_f = cm_axes(fig, 12.20, top, 3.60, H_F)
    draw_chi(ax_d)
    subs, _ = draw_constants(ax_e)
    draw_margins(ax_f)
    ax_d.yaxis.set_label_coords(-0.95 / ax_cm(ax_d)[0], 0.5)
    subs[0].yaxis.set_label_coords(-0.95 / ax_cm(subs[0])[0], 0.5)
    panel_letter(ax_d, "d", LABEL_X[0])
    panel_letter(ax_e, "e", LABEL_X[0])
    panel_letter(ax_f, "f", LABEL_X[1])
    return ax_d, ax_e, ax_f


def new_figure(height_cm):
    pad = 0.02                                   # pad_inches of the savefig calls
    fig = plt.figure(figsize=(style.DOUBLE - 2 * pad,
                              height_cm * style.CM - 2 * pad))
    # invisible frame: pins the tight bounding box to the full figure
    fig.add_artist(Rectangle((0, 0), 1, 1, transform=fig.transFigure,
                             fill=False, edgecolor="none", lw=0))
    return fig


# native size (cm) and margins (left, right, top, bottom; cm) of single panels
PANELS = {
    "accuracy": (draw_accuracy, RIGHT - LEFT, H_A, (0.10, 0.20, 0.15, 0.90)),
    "cases": (draw_cases, 2.60, H_B, (1.95, 2.45, 0.50, 0.90)),
    "elements": (draw_elements, RIGHT - 7.85, H_C, (0.10, 0.10, 0.10, 0.10)),
    "chi": (draw_chi, 7.70, H_D, (1.15, 0.20, 0.15, 0.90)),
    "constants": (draw_constants, 7.70, H_E, (1.15, 0.20, 0.15, 1.25)),
    "margins": (draw_margins, 3.60, H_F, (1.75, 1.35, 0.55, 0.90)),
}


def render_panel(name, scale, out):
    """One panel on its own, scaled by `scale` (fonts, lines and size)."""
    fn, w, h, (ml, mr, mt, mb) = PANELS[name]
    with scaled(scale):
        fig = plt.figure(figsize=((w + ml + mr) * scale * style.CM,
                                  (h + mt + mb) * scale * style.CM))
        ax = cm_axes(fig, ml * scale, mt * scale, w * scale, h * scale)
        fn(ax)
        fig.savefig(out, dpi=300, bbox_inches="tight", pad_inches=0.02 * scale)
        plt.close(fig)


# --------------------------------------------------------------- provenance
def print_provenance():
    d = data()
    bench, an, sc = d["bench"], d["analysis"], d["scan"]
    print("\nr3-T-B figure 2 provenance")
    print("=" * 72)
    cases = an["cases"]
    print(f"clean38 = development split of pyleos/sets/clean_set.json: "
          f"{len(cases)} cases, {sum(c['atoms'] for c in cases)} atoms, "
          f"{sum(c['labelled'] for c in cases)} with a per-atom reference, "
          f"{sum(c['wfs'] for c in cases)} WFs, {len(an['elements'])} elements "
          f"(counted from the case directories listed in {rel(ANALYSIS_JSON)})")
    print(f"plugin {an['plugin']}  sha256 {an['plugin_sha256'][:16]}...  "
          f"constants {an['constants']}")

    print("\n(a) totals of the per-case tables (python3 -m pyleos bench clean "
          "--split dev <option> --out <file>)")
    print(f"    {'method':40s} {'atoms correct':>16s} {'cases fully correct':>20s} "
          f"{'element-mean pass':>18s}  file  [option]")
    for key, _sty, label, fname, cmd in METHODS:
        s = bench[key]["summary"]
        print(f"    {label.replace('$', '').replace(chr(92), ''):40s} "
              f"{s['atoms_ok']:>5d}/{s['atoms_labelled']} "
              f"({100 * s['atoms_ok'] / s['atoms_labelled']:5.1f}%) "
              f"{s['full']:>13d}/{s['cases']} ({100 * s['full'] / s['cases']:5.1f}%)"
              f" {s['elem_ok']:>11d}/{s['cases']}  results/{fname}  [{cmd}]")

    clean, rest, n_all_ok = case_rows()
    print(f"\n(b) results/clean38_r3-T-B.tsv, clean38_flow.tsv, clean38_greedy.tsv: "
          f"incorrect = atoms_labelled - atoms_ok; {len(clean) + len(rest)} cases "
          f"plotted, {n_all_ok} cases fully correct with all three methods")
    print(f"    {'case':22s} {'labelled':>8s} {'this work':>10s} {'flow':>12s} "
          f"{'greedy':>12s}")
    for title, block in (("this work fully correct", clean),
                         ("this work not fully correct", rest)):
        print(f"    -- {title}: {len(block)} cases")
        for r in block:
            b, n = r["key"]
            cells = "".join(f" {r['wrong'][m]:4d} ({100 * r['wrong'][m] / r['n']:5.1f}%)"
                            for m in TRIO)
            print(f"    {b + '/' + n:22s} {r['n']:8d}{cells}")
    for m in TRIO:
        tot = sum(r["wrong"][m] for r in clean + rest)
        s = bench[m]["summary"]
        print(f"    total incorrect, {MLABEL[m]}: {tot} = {s['atoms_labelled']} - "
              f"{s['atoms_ok']}: {tot == s['atoms_labelled'] - s['atoms_ok']}")

    els = an["elements"]
    print(f"\n(c) {rel(ANALYSIS_JSON)}: frozen plugin re-run through the "
          f"evaluator's collect_clean / scrub / site_labels; per case equal to "
          f"results/clean38_r3-T-B.tsv (atoms_ok, atoms_labelled, wrong_sites)")
    print(f"    {len(els)} elements; correct/labelled (atoms in the cells, "
          f"cases containing the element):")
    for el, e in sorted(els.items(), key=lambda kv: (-kv[1]["labelled"], kv[0])):
        flag = "" if e["ok"] == e["labelled"] else "   <-- not all correct"
        print(f"    {el:3s} {e['ok']:4d}/{e['labelled']:<4d} ({e['atoms']:3d} atoms, "
              f"{len(e['cases']):2d} cases){flag}")
    n_ok = sum(e["ok"] for e in els.values())
    n_lab = sum(e["labelled"] for e in els.values())
    s = bench[THIS]["summary"]
    print(f"    sum {n_ok}/{n_lab}; equals panel (a): "
          f"{(n_ok, n_lab) == (s['atoms_ok'], s['atoms_labelled'])}; elements with "
          f"every labelled atom correct: "
          f"{sum(e['ok'] == e['labelled'] for e in els.values())}/{len(els)}; "
          f"atoms without a per-atom label: "
          f"{sum(e['atoms'] - e['labelled'] for e in els.values())}")

    rows = chi_elements()
    pb = [r for r in rows if r["block"] == "p"]
    print(f"\n(d) {rel(ANALYSIS_JSON)}: chi = Z_eff/rho^2 at the two ends of every "
          f"contested WF with a heteronuclear partner, for the ends that carry "
          f"lone pairs (plugin's _Shells.seen_from; records equal to "
          f"_analyse()['notes'] and, per element pair, to pair_margins()); "
          f"Allen values: pyleos.capacity.ENEG_ALLEN")
    print(f"    {'el':3s} {'block':5s} {'Allen':>6s} {'chi median':>10s} "
          f"{'quartiles':>15s} {'range (atoms)':>15s} {'range (ends)':>15s} "
          f"{'atoms':>5s} {'ends':>5s}  cases")
    for r in rows:
        print(f"    {r['el']:3s} {r['block']:5s} {r['allen']:6.3f} {r['chi']:10.3f} "
              f"{r['q1']:7.3f}-{r['q3']:<7.3f} {r['lo']:7.3f}-{r['hi']:<7.3f} "
              f"{r['end_lo']:7.3f}-{r['end_hi']:<7.3f} {r['atoms']:5d} {r['ends']:5d}"
              f"  {', '.join(r['cases'])}")
    print(f"    Spearman rank correlation, chi median vs Allen: p-block "
          f"({len(pb)} elements) "
          f"{spearman([r['allen'] for r in pb], [r['chi'] for r in pb]):.4f}; all "
          f"{len(rows)} elements "
          f"{spearman([r['allen'] for r in rows], [r['chi'] for r in rows]):.4f}")
    print(f"    order by chi:   {' < '.join(r['el'] for r in sorted(rows, key=lambda r: r['chi']))}")
    print(f"    order by Allen: {' < '.join(r['el'] for r in rows)}")

    print(f"\n(e) {rel(SCAN_JSON)} (produced {sc['produced']}): {sc['command']}")
    print(f"    {sc['copy']} ({sc.get('copies_dir', 'temporary directory')})")
    for name, s in sc["scan"].items():
        pts = s["points"]
        lo, hi = plateau(pts, s["default"])
        at0 = next(p for p in pts if p["value"] == s["default"])
        print(f"    {name}: default {s['default']:g} -> atoms {at0['atoms_ok']}/"
              f"{at0['atoms_labelled']}, cases fully correct {at0['full']}, "
              f"element-mean pass {at0['elem_ok']}; every case unchanged for "
              f"{lo:g} .. {hi:g}; {len(pts)} values, atoms correct "
              f"{min(p['atoms_ok'] for p in pts)}-{max(p['atoms_ok'] for p in pts)}")
        print("      value      " + " ".join(f"{p['value']:>7g}" for p in pts))
        print("      atoms ok   " + " ".join(f"{p['atoms_ok']:>7d}" for p in pts))
        print("      full cases " + " ".join(f"{p['full']:>7d}" for p in pts))
        print("      elem pass  " + " ".join(f"{p['elem_ok']:>7d}" for p in pts))
        seen = None
        for p in pts:
            sig = [(c["case"], c["atoms_ok"], c["status"]) for c in p["changed"]]
            if sig and sig != seen:
                print(f"      from {p['value']:g}: " + ", ".join(
                    f"{c['case']} {c['atoms_ok']}/{c['atoms_labelled']} {c['status']}"
                    for c in p["changed"]))
            seen = sig

    single, shared = margin_classes()
    wfs = an["wfs"]
    n_force = sum(w["rule"] == "force" for w in wfs)
    n_chi = sum(w["rule"] == "chi" for w in wfs)
    print(f"\n(f) {rel(ANALYSIS_JSON)}: the plugin's _analyse()['notes'] (the "
          f"records behind pair_margins): {len(wfs)} contested WFs with a "
          f"heteronuclear partner out of "
          f"{sum(c['contested'] for c in cases)} contested WFs "
          f"({sum(c['wfs'] for c in cases)} WFs in all); decided by force "
          f"{n_force}, by chi {n_chi}")
    print(f"    ratio = chi_margin (rule chi) or force_margin (rule force); inf = "
          f"Z_eff of the loser <= 0")
    for title, rows_f in (("one atom takes the pair", single),
                          ("like atoms share the pair", shared)):
        n_blk = sum(c["n"] for c in rows_f)
        print(f"    -- {title}: {n_blk} WFs, {len(rows_f)} classes")
        print(f"    {'class':10s} {'hosts':>6s} {'n':>4s} {'force':>5s} {'chi':>4s} "
              f"{'min':>8s} {'median':>8s} {'max':>8s} {'inf':>4s} "
              f"{'Allen w-l':>9s}  cases")
        for c in rows_f:
            allv = c["force"] + c["chi"]
            fin = [v for v in allv if not math.isinf(v)]
            f = lambda v: f"{v:8.3f}" if v is not None and not math.isinf(v) else f"{'inf' if v is not None else '-':>8s}"
            print(f"    {c['winner'] + '-' + c['loser']:10s} "
                  f"{'-'.join(str(h) for h in (c['hosts'][0], c['hosts'][-1])) if len(c['hosts']) > 1 else str(c['hosts'][0]):>6s} "
                  f"{c['n']:4d} {len(c['force']):5d} {len(c['chi']):4d} "
                  f"{f(min(fin) if fin else None)} {f(c['median'])} "
                  f"{f(max(fin) if fin else None)} {len(allv) - len(fin):4d} "
                  f"{c['allen']:+9.3f}  {', '.join(c['cases'])}")
    one = [w for w in wfs if w["hosts"] == 1]
    from pyleos.capacity import ENEG_ALLEN
    agree = sum(ENEG_ALLEN[w["winner"]] > ENEG_ALLEN[w["loser"]] for w in one)
    fin = [w["margin"] for w in wfs if w["margin"] != "inf"]
    for rule in ("force", "chi"):
        v = [w["margin"] for w in wfs if w["rule"] == rule and w["margin"] != "inf"]
        n_inf = sum(w["rule"] == rule and w["margin"] == "inf" for w in wfs)
        print(f"    rule {rule}: {len(v) + n_inf} WFs, finite ratios "
              f"{min(v):.3f} .. {max(v):.3f}, median (all) "
              f"{np.median(v + [math.inf] * n_inf):.3f}, infinite {n_inf}")
    print(f"    one-owner WFs won by the element with the larger Allen "
          f"electronegativity: {agree}/{len(one)}; the others: "
          + ", ".join(f"{k[0]}-{k[1]} x{v} ({k[2]})" for k, v in sorted(Counter(
              (w["winner"], w["loser"], w["case"]) for w in one
              if ENEG_ALLEN[w["winner"]] <= ENEG_ALLEN[w["loser"]]).items())))
    print(f"    ratios below 1.5: {sum(v < 1.5 for v in fin)}, below 2: "
          f"{sum(v < 2 for v in fin)}, smallest {min(fin):.3f} "
          f"({next((w['winner'] + '-' + w['loser'] + ' in ' + w['case']) for w in wfs if w['margin'] == min(fin))})")


class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, s):
        for st in self.streams:
            st.write(s)

    def flush(self):
        for st in self.streams:
            st.flush()


# --------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--recompute", action="store_true",
                    help="rebuild _cache/r3tb/fig2_analysis.json and fig2_scan.json")
    ap.add_argument("--scan-dir", type=Path,
                    help="keep the scan copies of the plugin in this directory "
                         "(outside the repository) instead of a temporary one")
    ap.add_argument("--panel", choices=sorted(PANELS),
                    help="render one panel only (with --out)")
    ap.add_argument("--scale", type=float, default=1.0,
                    help="scale of a single panel (fonts, lines, size)")
    ap.add_argument("--out", type=Path, help="output file of --panel")
    args = ap.parse_args()
    _OPTS.update(recompute=args.recompute, scan_dir=args.scan_dir)

    if args.panel:
        if not args.out:
            ap.error("--panel needs --out")
        render_panel(args.panel, args.scale, args.out)
        print(f"wrote {args.out}")
        return

    CACHE.mkdir(parents=True, exist_ok=True)
    with LOG.open("w") as log:
        stdout = sys.stdout
        sys.stdout = Tee(stdout, log)
        try:
            print(f"{Path(__file__).name}  {datetime.now():%Y-%m-%d %H:%M}  "
                  f"{' '.join(sys.argv[1:])}")
            data()
            height = H_TOP + H_BOTTOM
            fig = new_figure(height)
            build_top(fig, 0.0)
            build_bottom(fig, H_TOP)
            style.save(fig, NAME_OUT)
            for name, h, build in (("fig2_top", H_TOP, build_top),
                                   ("fig2_bottom", H_BOTTOM, build_bottom)):
                half = new_figure(h)
                build(half, 0.0)
                half.savefig(CACHE / f"{name}.png", dpi=300, bbox_inches="tight",
                             pad_inches=0.02)
                plt.close(half)
            print_provenance()
            print(f"\nwrote {rel(HERE / (NAME_OUT + '.pdf'))} and .png "
                  f"({style.DOUBLE / style.CM:.1f} x {height:.2f} cm), "
                  f"{rel(CACHE / 'fig2_top.png')} "
                  f"({style.DOUBLE / style.CM:.1f} x {H_TOP:.2f} cm), "
                  f"{rel(CACHE / 'fig2_bottom.png')} "
                  f"({style.DOUBLE / style.CM:.1f} x {H_BOTTOM:.2f} cm), "
                  f"log {rel(LOG)}")
        finally:
            sys.stdout = stdout


if __name__ == "__main__":
    main()
