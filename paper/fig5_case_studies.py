#!/usr/bin/env python3
"""Figure 5 - Case studies whose result is not fixed by stoichiometry.

(a) CH4 activation on Co4/CeO2 (initial, transition, final state): view of
    the active site; total oxidation state of the Co4 cluster treated as one
    fragment; number of Ce atoms at +3; oxidation state of the transferred
    H atom and of C; C-H and shortest Co-H distance of that H atom.  Two
    Wannierisations of the same three geometries are drawn (filled / open).
(b) The eight MH4 hydrides: oxidation state of the metal, of each H atom of
    one formula unit, and the shortest H-H distance, against the values of
    Comas-Vila et al.
(c) O vacancy in CeO2 (Ce4O8 -> Ce4O7): oxidation state of each Ce atom,
    and the position of the four Ce atoms around the vacancy.

Layout: (a) fills the left column; (b) above (c) in the right column; the
saved page is exactly 17.1 cm wide.

Run from this directory:   python3 fig5_case_studies.py [--recompute]

Every plotted number is computed with pyleos from the run directories
under cases/ (defaults: beta = 10, Allen electronegativity).  The baseline
("per-WF greedy rule", style key "nearest") is the OS(int) column that the
VASP patch writes to each EOSCAR.  The only hand-typed numbers are the
literature values of panel (b), in the table LIT_MH4.

"Occupation-weighted" values are a diagnostic, not a pyleos result: the same
min-cost-flow problem is solved without the cap on the total electron number
(every WF is shipped), and each WF is counted with its occupation from the
#WF block of the EOSCAR instead of as one (or two) electrons.  Set SHOW_OCC
to False to leave it out of the figure.

pyleos results are cached in _cache/fig5_{methane,hydrides,ceo2}.json with
the md5 of every input file and of the pyleos sources; a changed input (or
--recompute) triggers recomputation.  Writes fig5_case_studies.pdf / .png
and prints the tables behind fig5_caption.md.
"""
import argparse
import copy
import hashlib
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

import style
from style import C, MARKER, LABEL, formula

import matplotlib.pyplot as plt                     # noqa: E402
from matplotlib.lines import Line2D                 # noqa: E402
from matplotlib.patches import Circle, Patch, Rectangle  # noqa: E402
from matplotlib.ticker import (FixedLocator, FuncFormatter,  # noqa: E402
                               MultipleLocator, NullLocator)

REPO = style.REPO
DATA = REPO / "cases"                       # <case>/<run>/ directories
CACHE = style.HERE / "_cache"
sys.path.insert(0, str(REPO))
from pyleos.parsers import load_case, parse_poscar  # noqa: E402
from pyleos.transport import assign_case, reliability_margins  # noqa: E402

PYLEOS_SRC = [REPO / "pyleos" / f for f in
              ("parsers.py", "transport.py", "capacity.py")]
CACHE_VERSION = 8           # bump when the layout of the cached data changes
NAME = "fig5_case_studies"
BASE = "nearest"            # style.py key of the per-WF greedy baseline
SHOW_OCC = True             # draw the occupation-weighted diagnostic

# ---------------------------------------------------------------- inputs --
# (a) Methane activation: two Wannierisations of the same three geometries.
#   "gs"     - 2026-07-11 reruns (v31-mag*), number of WFs set per spin
#              channel, in the magnetic ground state listed in
#              results/methane_groundstate_verdict.tsv
#   "legacy" - runs (v3) with the same number of WFs in both spin channels;
#              the values quoted in docs/RESULTS.md section 6 come from these
CH4 = {im: DATA / f"ce18o36co4ch4_{im}" for im in ("is", "ts", "fs")}
CH4_SETS = {
    "gs": {"is": CH4["is"] / "v31-mag6", "ts": CH4["ts"] / "v31-mag6",
           "fs": CH4["fs"] / "v31-mag10"},
    "legacy": {"is": CH4["is"] / "v3", "ts": CH4["ts"] / "v3",
               "fs": CH4["fs"] / "v3"},
}
# replicas of the ground-state runs (second initial moment); reported only
CH4_REPLICAS = {"is": CH4["is"] / "v31-mag10", "ts": CH4["ts"] / "v31-mag10"}
CH4_VERDICT = REPO / "results" / "methane_groundstate_verdict.tsv"
CH4_SET_LABEL = {"gs": r"$N_{\mathrm{WF}}$ set per spin channel",
                 "legacy": r"equal $N_{\mathrm{WF}}$ in both channels"}
IMAGES = [("is", "Initial state"), ("ts", "Transition state"),
          ("fs", "Final state")]
CE3_POP = 8.5               # a Ce with more WF electrons than this is Ce(+3)

# (b) MH4 hydrides, ordered by the oxidation state of the metal; run
# directory cases/<name>/v3.
MH4_RUN = "v3"
MH4 = [("mgh4", "Mg"), ("cah4", "Ca"), ("srh4", "Sr"), ("sch4", "Sc"),
       ("yh4", "Y"), ("zrh4", "Zr"), ("ceh4", "Ce"), ("thh4", "Th")]
H2_CUT = 1.0                # two H atoms closer than this (A) are an H2 unit

# =============================================================================
# LITERATURE VALUES - the only hand-typed numbers in this script.
# G. Comas-Vila et al., J. Chem. Theory Comput. 2025, 21, 7075-7086, as
# transcribed in docs/comas_vila_2025_notes.md (N<k> = line k of that file;
# p. = journal page).  Their structures are I4/mmm with two H sites per
# formula unit: two hydridic H_b and one H_a-H_a unit.  Values per atom:
#   case: (metal, H_b, H_a, notes line, page)
# =============================================================================
LIT_MH4 = {
    "mgh4": (2, -1, 0.0, 352, "7081"),         # M(+2), Hb(-1), neutral Ha-Ha
    "cah4": (2, -1, 0.0, 352, "7081"),
    "srh4": (2, -1, 0.0, 352, "7081"),
    "sch4": (3, -1, -0.5, 353, "7081-7082"),   # M(+3), Hb(-1), Ha(-1/2)
    "yh4": (3, -1, -0.5, 353, "7081-7082"),
    "zrh4": (4, -1, -1.0, 354, "7082"),        # M4+ (H-)4
    "ceh4": (4, -1, -1.0, 354, "7082"),
    "thh4": (4, -1, -1.0, 354, "7082"),
}

# (c) CeO2 with and without one O vacancy; run directory cases/<name>/v3.
CEO2_RUN = "v3"
CEO2 = [("ce4o8", "Ce4O8"), ("ce4o7", "Ce4O7")]

PROV = []       # (panel, quantity, value, source)


def prov(panel, quantity, value, source):
    PROV.append((panel, quantity, str(value), str(source)))


def rel(p):
    return str(Path(p).relative_to(REPO))


# ------------------------------------------------------------ cache layer --
def md5(path):
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def cached(name, inputs, compute, force):
    """Return compute() from _cache/fig5_<name>.json when the md5 of every
    input (and of the pyleos sources) is unchanged; recompute otherwise."""
    path = CACHE / f"fig5_{name}.json"
    sig = {rel(p): md5(p) for p in list(inputs) + PYLEOS_SRC}
    sig["cache_version"] = CACHE_VERSION
    if path.is_file() and not force:
        blob = json.loads(path.read_text())
        if blob.get("inputs") == sig:
            return blob["data"]
    t0 = time.time()
    print(f"[fig5] computing {name} with pyleos ...", flush=True)
    data = compute()
    CACHE.mkdir(exist_ok=True)
    path.write_text(json.dumps(
        {"generated_by": f"{NAME}.py", "inputs": sig, "data": data},
        indent=1))
    print(f"[fig5] {name}: {time.time() - t0:.0f} s -> {rel(path)}")
    return data


def eoscar_header(d, key):
    """Value of a '# key = value' line in the header of an EOSCAR."""
    for ln in (Path(d) / "EOSCAR").read_text().splitlines():
        if ln.startswith("# " + key) and "=" in ln:
            return float(ln.split("=", 1)[1])
    raise KeyError(f"{d}: no {key} in the EOSCAR header")


def case_files(d):
    return [p for p in (Path(d) / f for f in
                        ("EOSCAR", "POSCAR", "wannier90.wout",
                         "wannier90.1.wout", "wannier90.2.wout"))
            if p.is_file()]


# --------------------------------------------------------------- geometry --
SHIFTS = np.array(list(itertools.product((-1, 0, 1), repeat=3)), float)


def min_image_vec(st, frac_from, frac_to):
    """Shortest periodic vector (A) between two fractional positions."""
    lat = np.array(st.lattice)
    d = np.array(frac_to, float) - np.array(frac_from, float)
    d -= np.round(d)
    v = (d + SHIFTS) @ lat
    return v[np.argmin(np.linalg.norm(v, axis=1))]


def min_image_dist(st, i, j):
    return float(np.linalg.norm(min_image_vec(st, st.frac[i], st.frac[j])))


def image_dists(st, i, j, cut):
    """All i-j distances below cut, periodic images included (i = j allowed)."""
    lat = np.array(st.lattice)
    d = np.array(st.frac[j]) - np.array(st.frac[i])
    n = np.linalg.norm((d - np.round(d) + SHIFTS) @ lat, axis=1)
    return [float(x) for x in n if 1e-6 < x < cut]


def idx(elements, el):
    return [i for i, e in enumerate(elements) if e == el]


# ------------------------------------------- occupation-weighted diagnostic --
def uncapped_population(case, fragments=None):
    """Min-cost flow without the electron-number cap (every WF shipped);
    returns per assignment unit the number of electrons counted as integers
    and counted with the WF occupations of the EOSCAR #WF block.  Units are
    the unmerged atoms in POSCAR order followed by the fragments."""
    free = copy.deepcopy(case)
    free.eoscar.nelect_wann = None            # no cap: target = full supply
    res = assign_case(free, fragments=fragments)
    occ = {(w.spin, w.index): w.occ for w in case.eoscar.wf_detail}
    full = 2.0 if case.ispin == 1 else 1.0
    merged = {a for f in (fragments or []) for a in f}
    units = [[a] for a in range(len(case.structure.elements))
             if a not in merged] + [sorted(f) for f in (fragments or [])]
    count = [0] * len(units)
    pop = [0.0] * len(units)
    for w in res.wfs:
        o = occ[(w.channel, w.index)]
        for u, n in w.dests:
            count[u] += n
            pop[u] += o * n / full
    return units, count, pop


# ---------------------------------------------------- panel (a): methane --
def run_methane(d):
    case = load_case(d)
    st, el = case.structure, case.structure.elements
    eo = case.eoscar
    co, hs, cs, ce, ox = (idx(el, x) for x in ("Co", "H", "C", "Ce", "O"))
    res = assign_case(case, fragments=[co])             # Co4 as one unit
    units, count, pop = uncapped_population(case, fragments=[co])
    unit_of = {m[0]: u for u, m in enumerate(units) if len(m) == 1}
    frag = len(units) - 1
    zv = eo.zval
    occ = [sum(w.occ for w in eo.wf_detail if w.spin == s) for s in (1, 2)]
    hole = [c - p for c, p in zip(count, pop)]
    g = eo.greedy_os                                    # EOSCAR OS(int)
    magmom = [ln.split("=", 1)[1].strip()
              for ln in (Path(d) / "INCAR").read_text().splitlines()
              if ln.strip().startswith("MAGMOM")]
    wf_occ = {(w.spin, w.index): w.occ for w in eo.wf_detail}
    near_spin = {i: [sum(w.nearest == i and w.channel == s for w in res.wfs)
                     for s in (1, 2)] for i in ce}
    low = min(eo.wf_detail, key=lambda w: w.occ)
    return {
        "dir": rel(d), "co_atoms": [i + 1 for i in co],
        "h_atoms": [i + 1 for i in hs], "ce_atoms": [i + 1 for i in ce],
        "n_wf": [len(ch.centres) for ch in case.channels],
        "occ_spin": occ, "moment": occ[0] - occ[1],
        "nelect_wann": eo.nelect_wann, "dne": eo.dne,
        "supply": res.supply_total, "target": res.target,
        "ghost": res.ghost, "stray": len(res.stray),
        "magmom": magmom[-1] if magmom else "(not set: VASP default)",
        # capped flow (the pyleos result)
        "frag_ne": res.frag_ne[0], "frag_os": res.frag_os[0],
        "h_flow": [res.os_int[i] for i in hs],
        "c_flow": [res.os_int[i] for i in cs],
        "ce_flow": [res.os_int[i] for i in ce],
        "o_flow": sorted({res.os_int[i] for i in ox}),
        # WFs that the cap leaves out: nearest atom and occupation
        "cap_cut": sorted([f"{el[w.nearest]}{w.nearest + 1}",
                           wf_occ[(w.channel, w.index)]] for w in res.wfs
                          if sum(n for _, n in w.dests) < w.supply),
        "lowest_occ": [f"{el[low.at1 - 1]}{low.at1}", low.spin, low.occ],
        "ce_near_spin": [near_spin[i] for i in ce],
        # greedy column of the EOSCAR
        "co_greedy": [g[i] for i in co], "h_greedy": [g[i] for i in hs],
        "c_greedy": [g[i] for i in cs], "ce_greedy": [g[i] for i in ce],
        "ce_greedy_ne": [eo.greedy_ne[i] for i in ce],
        # uncapped flow, integer count and occupation-weighted
        "co4_count": count[frag], "co4_pop": pop[frag],
        "co4_os_occ": 4 * zv["Co"] - pop[frag],
        "ce_count": [count[unit_of[i]] for i in ce],
        "ce_pop": [pop[unit_of[i]] for i in ce],
        "h_os_occ": [zv["H"] - pop[unit_of[i]] for i in hs],
        "c_os_occ": [zv["C"] - pop[unit_of[i]] for i in cs],
        "holes": {"Co4": hole[frag],
                  "Ce": sum(hole[unit_of[i]] for i in ce),
                  "O": sum(hole[unit_of[i]] for i in ox),
                  "CH4": sum(hole[unit_of[i]] for i in cs + hs)},
        # geometry
        "d_ch": [min_image_dist(st, h, cs[0]) for h in hs],
        "d_coh": [min(min_image_dist(st, h, k) for k in co) for h in hs],
        "d_ce_co": [min(min_image_dist(st, i, k) for k in co) for i in ce],
    }


def compute_methane():
    out = {}
    for s, dirs in CH4_SETS.items():
        for image, d in dirs.items():
            out[f"{s}/{image}"] = run_methane(d)
    for image, d in CH4_REPLICAS.items():
        out[f"replica/{image}"] = run_methane(d)
    return out


def read_tsv(path):
    lines = [ln for ln in Path(path).read_text().splitlines() if ln.strip()]
    head = lines[0].split("\t")
    return [dict(zip(head, ln.split("\t"))) for ln in lines[1:]]


def load_methane(force):
    dirs = [d for s in CH4_SETS.values() for d in s.values()] + \
        list(CH4_REPLICAS.values())
    inputs = [f for d in dirs for f in case_files(d)] + \
        [d / "INCAR" for d in dirs]
    raw = cached("methane", inputs, compute_methane, force)
    # transferred H = the H whose C-H distance is largest in the final state
    fs = raw["gs/fs"]
    k = int(np.argmax(fs["d_ch"]))
    data = {"h_atom": fs["h_atoms"][k], "k": k, "raw": raw, "sets": {}}
    for s in CH4_SETS:
        rows = []
        for image, _ in IMAGES:
            d = raw[f"{s}/{image}"]
            assert d["co_atoms"] == [55, 56, 57, 58], d["co_atoms"]
            rows.append({"co4_flow": d["frag_os"], "co4_occ": d["co4_os_occ"],
                         "h_flow": d["h_flow"][k], "c_flow": d["c_flow"][0],
                         "d_ch": d["d_ch"][k], "d_coh": d["d_coh"][k],
                         "n_ce3_flow": d["ce_flow"].count(3),
                         "n_ce3_occ": sum(p > CE3_POP for p in d["ce_pop"]),
                         "raw": d})
        data["sets"][s] = rows
    a, b = (data["sets"][s] for s in CH4_SETS)
    assert all(abs(p["d_ch"] - q["d_ch"]) < 1e-6 for p, q in zip(a, b)), \
        "the two Wannierisations must share their geometries"
    return data


def report_methane(ch4):
    raw, k, h = ch4["raw"], ch4["k"], ch4["h_atom"]
    print("\n(a) METHANE ON Co4/CeO2 - what differs between the data sets")
    print("    500 valence electrons; Co4 fragment = atoms 55-58; "
          f"transferred H = atom {h}")
    keys = [f"{s}/{im}" for s in CH4_SETS for im, _ in IMAGES] + \
        [f"replica/{im}" for im in CH4_REPLICAS]
    rows = [
        ("directory (cases/...)", lambda d: d["dir"].split("cases/")[1]),
        ("MAGMOM in INCAR", lambda d: d["magmom"][:38]),
        ("N_WF spin 1 / 2", lambda d: f"{d['n_wf'][0]} / {d['n_wf'][1]}"),
        ("occupation spin 1 / 2", lambda d: f"{d['occ_spin'][0]:.2f} / "
                                            f"{d['occ_spin'][1]:.2f}"),
        ("moment (mu_B)", lambda d: f"{d['moment']:.2f}"),
        ("NELECT(Wannier), |dNe|", lambda d: f"{d['nelect_wann']:.3f}, "
                                             f"{d['dne']:.1e}"),
        ("surplus (ghost) e / stray WFs", lambda d: f"{d['ghost']} / "
                                                    f"{d['stray']}"),
        ("-- flow with cap (pyleos) --", lambda d: ""),
        ("Co4 n_e -> OS", lambda d: f"{d['frag_ne']} -> {d['frag_os']:+d}"),
        (f"H{h} / C", lambda d: f"{d['h_flow'][k]:+d} / {d['c_flow'][0]:+d}"),
        ("Ce at +3 (atoms)", lambda d: str([a for a, v in
                                            zip(d["ce_atoms"], d["ce_flow"])
                                            if v == 3])),
        ("WFs left out by the cap (occ.)", lambda d: " ".join(
            f"{a}({o:.2f})" for a, o in d["cap_cut"])[:40]),
        ("  ... continued", lambda d: " ".join(
            f"{a}({o:.2f})" for a, o in d["cap_cut"])[40:]),
        ("least occupied WF: atom, spin, occ.", lambda d: "{} spin {} {:.3f}"
         .format(*d["lowest_occ"])),
        ("WFs nearest Ce6, Ce7, Ce8 (up+down)", lambda d: " ".join(
            f"{d['ce_near_spin'][i][0]}+{d['ce_near_spin'][i][1]}"
            for i in (5, 6, 7))),
        ("-- greedy column of EOSCAR --", lambda d: ""),
        ("sum over the four Co", lambda d: f"{sum(d['co_greedy']):+d} "
                                           f"{d['co_greedy']}"),
        (f"H{h} / C", lambda d: f"{d['h_greedy'][k]:+d} / "
                                f"{d['c_greedy'][0]:+d}"),
        ("Ce not +4 (atom:state)", lambda d: " ".join(
            f"{a}:{v:+d}" for a, v in zip(d["ce_atoms"], d["ce_greedy"])
            if v != 4)),
        ("electrons on Ce beyond 18 x 8", lambda d: "{:.2f}".format(
            sum(d["ce_greedy_ne"]) - 8 * len(d["ce_greedy_ne"]))),
        ("greedy n_e of Ce6, Ce7, Ce8", lambda d: " ".join(
            f"{d['ce_greedy_ne'][i]:.2f}" for i in (5, 6, 7))),
        ("-- flow without cap --", lambda d: ""),
        ("Ce holding 9 WFs (atoms)", lambda d: str([a for a, c in
                                                    zip(d["ce_atoms"],
                                                        d["ce_count"])
                                                    if c == 9])),
        ("their occupation (e)", lambda d: " ".join(
            f"{p:.2f}" for p, c in zip(d["ce_pop"], d["ce_count"])
            if c == 9)),
        ("Ce holding 10 WFs: occupation", lambda d: " ".join(
            f"Ce{a}:{p:.2f}" for a, p, c in zip(d["ce_atoms"], d["ce_pop"],
                                                d["ce_count"]) if c > 9)
            or "-"),
        ("Co4: WFs -> OS", lambda d: f"{d['co4_count']} -> "
                                     f"{60 - d['co4_count']:+d}"),
        ("Co4: occupation -> OS", lambda d: f"{d['co4_pop']:.2f} -> "
                                            f"{d['co4_os_occ']:+.2f}"),
        ("missing e on Co4 / Ce / O / CH4", lambda d: " / ".join(
            f"{d['holes'][x]:.2f}" for x in ("Co4", "Ce", "O", "CH4"))),
        (f"H{h} / C by occupation", lambda d: f"{d['h_os_occ'][k]:+.2f} / "
                                              f"{d['c_os_occ'][0]:+.2f}"),
        ("-- geometry --", lambda d: ""),
        ("nearest Co of Ce6, Ce7, Ce8 (A)", lambda d: " ".join(
            f"{d['d_ce_co'][i]:.2f}" for i in (5, 6, 7))),
        (f"C-H{h} / shortest Co-H{h} (A)", lambda d: f"{d['d_ch'][k]:.3f} / "
                                                     f"{d['d_coh'][k]:.3f}"),
    ]
    w0 = max(len(r[0]) for r in rows) + 2
    for block in (keys[:3], keys[3:6], keys[6:]):
        print()
        print(" " * (w0 + 2) + "".join(f"{key:<42}" for key in block))
        for name, fn in rows:
            print(f"  {name:<{w0}}" + "".join(f"{fn(raw[key]):<42}"
                                              for key in block))
    print("\n  ground-state check (" + rel(CH4_VERDICT) + ")")
    for row in read_tsv(CH4_VERDICT):
        print(f"    {row['image']:<6} {row['init_guess']:<11} E0 = "
              f"{row['E0_eV']} eV  moment {row['mag_converged']}  "
              f"{row['verdict']}")
    for s in CH4_SETS:
        r = ch4["sets"][s]
        prov("a", f"[{s}] Co4 fragment OS, flow (IS/TS/FS)",
             " / ".join(f"{x['co4_flow']:+d}" for x in r),
             "assign_case(fragments=[[54,55,56,57]])")
        prov("a", f"[{s}] Co4 OS, occupation-weighted",
             " / ".join(f"{x['co4_occ']:+.2f}" for x in r),
             "uncapped flow x #WF occupations")
        prov("a", f"[{s}] Ce(+3) count, flow | occupation > 8.5 e",
             " / ".join(str(x["n_ce3_flow"]) for x in r) + " | " +
             " / ".join(str(x["n_ce3_occ"]) for x in r), "same")
        prov("a", f"[{s}] H{h} / C, flow",
             " ; ".join(f"{x['h_flow']:+d} / {x['c_flow']:+d}" for x in r),
             "same")
    r = ch4["sets"]["gs"]
    prov("a", f"C-H{h} (A)", " / ".join(f"{x['d_ch']:.3f}" for x in r),
         "POSCAR, minimum image")
    prov("a", f"shortest Co-H{h} (A)",
         " / ".join(f"{x['d_coh']:.3f}" for x in r), "POSCAR")


# --------------------------------------------------- panel (b): hydrides --
def h_sites(st, hs, ms):
    """Group the H atoms of a cell into symmetry classes by their neighbour
    distances; return per class the atoms, the nearest H-H distance, and the
    pairs of mutually nearest atoms when the nearest H belongs to the class."""
    sig, near = {}, {}
    for i in hs:
        hh = sorted(d for j in hs for d in image_dists(st, i, j, 2.6))
        mh = sorted(d for j in ms for d in image_dists(st, i, j, 3.0))
        sig[i] = (tuple(round(x, 2) for x in hh), tuple(round(x, 2)
                                                         for x in mh))
        cand = [(min(image_dists(st, i, j, 9.0)), j) for j in hs]
        near[i] = min(cand)
    classes = {}
    for i in hs:
        classes.setdefault(sig[i], []).append(i)
    out = []
    for atoms in classes.values():
        pairs = sorted({tuple(sorted((i, near[i][1]))) for i in atoms
                        if near[near[i][1]][1] == i and near[i][1] != i})
        out.append({"atoms": atoms, "d_hh": near[atoms[0]][0],
                    "pairs": pairs})
    return sorted(out, key=lambda c: c["d_hh"])


def compute_hydrides():
    out = {}
    for name, metal in MH4:
        d = DATA / name / MH4_RUN
        case = load_case(d)
        st, el = case.structure, case.structure.elements
        eo = case.eoscar
        res = reliability_margins(case, res=assign_case(case))
        hs, ms = idx(el, "H"), idx(el, metal)
        sites = h_sites(st, hs, ms)
        out[name] = {
            "dir": rel(d), "metal": metal, "z": len(ms),
            "m_flow": [res.os_int[i] for i in ms],
            "m_greedy": [eo.greedy_os[i] for i in ms],
            "h_flow": {str(i + 1): res.os_int[i] for i in hs},
            "h_greedy": {str(i + 1): eo.greedy_os[i] for i in hs},
            "sites": [{"atoms": [i + 1 for i in c["atoms"]],
                       "d_hh": c["d_hh"],
                       "pairs": [[i + 1, j + 1] for i, j in c["pairs"]]}
                      for c in sites],
            "d_mh": min(min(image_dists(st, i, j, 9.0)) for i in hs
                        for j in ms),
            "ghost": res.ghost, "stray": len(res.stray),
            "n_wf": len(res.wfs), "dne": eo.dne,
            "dne_warn": eoscar_header(d, "EOS_NELECT_WARN"),
            "nelect_wann": eo.nelect_wann,
            "r_site": res.case_conf, "r_elem": res.case_conf_elem,
            "r_h": {str(i + 1): res.atom_conf[i] for i in hs},
        }
    return out


def load_hydrides(force):
    inputs = [f for n, _ in MH4 for f in case_files(DATA / n / MH4_RUN)]
    raw = cached("hydrides", inputs, compute_hydrides, force)
    rows = []
    for name, metal in MH4:
        d = raw[name]
        z, hf = d["z"], d["h_flow"]
        two = len(d["sites"]) == 2
        d_min = min(site["d_hh"] for site in d["sites"])
        row = {"name": name, "metal": metal, "two_sites": two, "raw": d,
               "m_flow": d["m_flow"], "m_greedy": d["m_greedy"],
               "stray": d["stray"], "lit": LIT_MH4[name],
               "lost": d["dne"] > d["dne_warn"],
               "d_min": d_min, "h2": d_min < H2_CUT}
        if two:
            # site with the short H-H contact = H_a (pairs), other = H_b
            ha, hb = d["sites"]
            pairs = sorted(tuple(sorted(hf[str(i)] for i in p))
                           for p in ha["pairs"])
            assert len(ha["pairs"]) == z and len(set(pairs)) == 1, pairs
            hb_vals = sorted(hf[str(i)] for i in hb["atoms"])
            assert len(set(hb_vals)) == 1 and len(hb_vals) == 2 * z
            assert abs(ha["d_hh"] - d_min) < 1e-9, "shortest contact is Ha-Ha"
            row.update(hb=hb_vals[:2], ha=list(pairs[0]),
                       d_ha=ha["d_hh"], d_hb=hb["d_hh"])
        else:
            # one H site: the four atoms form two equivalent pairs of mutual
            # nearest neighbours.  The pair with the lower total is listed
            # (and drawn) first.  This order is a drawing convention only:
            # the pairs are equivalent, neither is the H_b or the H_a site of
            # the literature structure, and none is an H2 unit (the caption
            # states the rule).
            (site,) = d["sites"]
            assert z == 1 and len(site["pairs"]) == 2, site
            pairs = sorted((sorted(hf[str(i)] for i in p)
                            for p in site["pairs"]), key=lambda v: (sum(v), v))
            row.update(h_pairs=pairs, h_all=pairs[0] + pairs[1],
                       d_h=site["d_hh"])
        rows.append(row)
    return rows


def report_hydrides(rows):
    print("\n(b) MH4 HYDRIDES - per formula unit; literature = Comas-Vila "
          "et al. (LIT_MH4)")
    print("    case   H sites  nearest H-H (A)        M flow/lit/greedy   "
          "H atoms, flow              lit: Hb, Ha     greedy H         "
          "stray/N_WF ghost |dNe|    R_elem R_site")
    for r in rows:
        d, lit = r["raw"], r["lit"]
        if r["two_sites"]:
            dist = f"Ha-Ha {r['d_ha']:.3f}, Hb {r['d_hb']:.3f}"
            hf = (f"Hb {r['hb'][0]:+d},{r['hb'][1]:+d}; Ha-Ha "
                  f"{r['ha'][0]:+d},{r['ha'][1]:+d}")
        else:
            dist = f"all H {r['d_h']:.3f}"
            hf = "2 pairs: " + "; ".join(",".join(f"{v:+d}" for v in pr)
                                         for pr in r["h_pairs"])
        hg = ",".join(f"{v:+d}" for v in d["h_greedy"].values())
        mg = "/".join(f"{v:+d}" for v in sorted(set(r["m_greedy"])))
        print(f"    {r['name']:<6} {len(d['sites']):^7}  {dist:<22} "
              f"{r['m_flow'][0]:+d} / {lit[0]:+d} / {mg:<9} {hf:<26} "
              f"{lit[1]:+d}, {lit[2]:+.1f}   {hg:<22} "
              f"{d['stray']:>3}/{d['n_wf']:<3} {d['ghost']:>4}  "
              f"{d['dne']:.1e}"
              f"{'*' if r['lost'] else ' '} {d['r_elem']:5.1f} "
              f"{d['r_site']:5.1f}")
        prov("b", f"{r['name']}: M flow / lit / greedy; nearest H-H",
             f"{r['m_flow'][0]:+d} / {lit[0]:+d} / {mg}; {dist}",
             f"{d['dir']}; notes N{lit[3]}, p. {lit[4]}")
    one = [r["name"] for r in rows if not r["two_sites"]]
    print(f"    cells with ONE H site (not the two-site I4/mmm structure of "
          f"the literature): {', '.join(one)}")
    print("      their four H atoms form two equivalent pairs; the pair with "
          "the lower total is drawn on the left (drawing convention, no site "
          "assignment); H-H in a pair: " + ", ".join(
              f"{r['name']} {r['d_h']:.2f} A" for r in rows
              if not r["two_sites"]))
    stray = [r for r in rows if r["stray"]]
    print("    dagger = stray WFs (spread > 3 A^2 or centre > 1.6 A from "
          "every atom), of N_WF: " + ", ".join(
              "{} {}/{} ({:.0f} %)".format(
                  r["name"], r["stray"], r["raw"]["n_wf"],
                  100 * r["stray"] / r["raw"]["n_wf"]) for r in stray))
    prov("b", "dagger: stray WFs / N_WF",
         ", ".join(f"{r['name']} {r['stray']}/{r['raw']['n_wf']}"
                   for r in stray),
         "len(res.stray) / len(res.wfs), assign_case defaults")
    thh = next(r for r in rows if len(set(r["m_greedy"])) > 1)
    prov("b", f"{thh['name']}: greedy, one square per metal atom of the cell",
         " and ".join(f"{v:+d}" for v in thh["m_greedy"]),
         "OS(int) of the two metal atoms; drawn unjoined")
    print("    * |dNe| above EOS_NELECT_WARN of the EOSCAR header; H2 unit = "
          f"H-H below {H2_CUT} A: "
          + ", ".join(r["name"] for r in rows if r["h2"]))
    print("    mean H oxidation state, flow / literature: " + "  ".join(
        "{} {:+.2f}/{:+.2f}".format(
            r["name"], float(np.mean(list(r["raw"]["h_flow"].values()))),
            (2 * r["lit"][1] + 2 * r["lit"][2]) / 4) for r in rows))


# ------------------------------------------------------ panel (c): CeO2 ---
def compute_ceo2():
    out = {}
    for name, _ in CEO2:
        d = DATA / name / CEO2_RUN
        case = load_case(d)
        el, eo = case.structure.elements, case.eoscar
        res = reliability_margins(case, res=assign_case(case))
        units, count, pop = uncapped_population(case)
        ce = idx(el, "Ce")
        occ = [sum(w.occ for w in eo.wf_detail if w.spin == s)
               for s in (1, 2)]
        out[name] = {
            "dir": rel(d), "ce_atoms": [i + 1 for i in ce],
            "ce_flow": [res.os_int[i] for i in ce],
            "ce_greedy": [eo.greedy_os[i] for i in ce],
            "ce_greedy_ne": [eo.greedy_ne[i] for i in ce],
            "ce_count": [count[i] for i in ce],
            "ce_os_occ": [eo.zval["Ce"] - pop[i] for i in ce],
            "o_flow": sorted({res.os_int[i] for i in idx(el, "O")}),
            "o_greedy": sorted({eo.greedy_os[i] for i in idx(el, "O")}),
            "n_wf": [len(ch.centres) for ch in case.channels],
            "occ_spin": occ, "moment": occ[0] - occ[1],
            "nelect_wann": eo.nelect_wann, "dne": eo.dne,
            "supply": res.supply_total, "ghost": res.ghost,
            "stray": len(res.stray),
            "r_site": res.case_conf, "r_elem": res.case_conf_elem,
            "sum_os": sum(res.os_int),
        }
    # vacancy = the O of ce4o8 that has no counterpart in ce4o7
    full = parse_poscar(DATA / "ce4o8" / CEO2_RUN / "POSCAR")
    vac = parse_poscar(DATA / "ce4o7" / CEO2_RUN / "POSCAR")
    missing = [full.frac[i] for i in idx(full.elements, "O")
               if all(np.linalg.norm(min_image_vec(full, full.frac[i],
                                                   vac.frac[j])) > 0.2
                      for j in idx(vac.elements, "O"))]
    assert len(missing) == 1, missing
    ce = idx(vac.elements, "Ce")
    out["vacancy"] = {
        "frac": list(missing[0]),
        "d_ce": [float(np.linalg.norm(min_image_vec(vac, missing[0],
                                                    vac.frac[i])))
                 for i in ce],
        "vec_ce": [list(min_image_vec(vac, missing[0], vac.frac[i]))
                   for i in ce],
        "max_shift": max(float(np.linalg.norm(min_image_vec(
            full, full.frac[i], vac.frac[i]))) for i in ce),
        "ce_o_mean": [float(np.mean(sorted(
            d for j in idx(vac.elements, "O")
            for d in image_dists(vac, i, j, 3.0)))) for i in ce],
        "ce_o_n": [len([d for j in idx(vac.elements, "O")
                        for d in image_dists(vac, i, j, 3.0)]) for i in ce],
    }
    return out


def load_ceo2(force):
    inputs = [f for n, _ in CEO2 for f in case_files(DATA / n / CEO2_RUN)]
    return cached("ceo2", inputs, compute_ceo2, force)


def report_ceo2(ce):
    print("\n(c) O VACANCY IN CeO2")
    for name, _ in CEO2:
        d = ce[name]
        print(f"    {name}: N_WF {d['n_wf'][0]}/{d['n_wf'][1]}, occupation "
              f"{d['occ_spin'][0]:.2f}/{d['occ_spin'][1]:.2f} (moment "
              f"{d['moment']:.2f}), NELECT(Wannier) {d['nelect_wann']:.3f}, "
              f"|dNe| {d['dne']:.1e}, surplus e {d['ghost']}, stray "
              f"{d['stray']}, sum of OS {d['sum_os']}, R_elem "
              f"{d['r_elem']:.1f}, R_site {d['r_site']:.1f}")
        print(f"        Ce1-Ce4  flow {d['ce_flow']}  greedy "
              f"{d['ce_greedy']}  occupation-weighted "
              f"{[round(x, 2) for x in d['ce_os_occ']]}  (WFs per Ce "
              f"{d['ce_count']}; greedy n_e "
              f"{[round(x, 2) for x in d['ce_greedy_ne']]});  O flow "
              f"{d['o_flow']} greedy {d['o_greedy']}")
        prov("c", f"{name}: Ce flow / greedy / occupation-weighted",
             f"{d['ce_flow']} / {d['ce_greedy']} / "
             f"{[round(x, 2) for x in d['ce_os_occ']]}", d["dir"])
    v = ce["vacancy"]
    print(f"    vacancy at fractional {np.round(v['frac'], 3)}; Ce-vacancy "
          f"distances (A) {[round(x, 3) for x in v['d_ce']]}; largest Ce "
          f"displacement between the two POSCARs {v['max_shift']:.3f} A; "
          f"Ce-O in Ce4O7: {v['ce_o_n']} neighbours, mean "
          f"{[round(x, 3) for x in v['ce_o_mean']]} A")
    prov("c", "Ce-vacancy distances (A)",
         [round(x, 3) for x in v["d_ce"]], "POSCAR of ce4o8 and ce4o7")
    # the sketch is a projection along [001]: the drawn dashed lines are
    # shorter than the distances that the label "4 x 2.34 A" gives
    drawn = sorted({round(float(np.hypot(x[0], x[1])), 2)
                    for x in v["vec_ce"]})
    print(f"    sketch: Ce-vacancy vectors projected along [001], length on "
          f"the page {drawn} A (label of the sketch: 3-D distance "
          f"{sorted({round(x, 2) for x in v['d_ce']})} A)")
    prov("c", "sketch: projected Ce-vacancy length (A), not labelled",
         drawn, "x, y components of the same vectors")


# ================================================================ drawing ==
# All positions are in centimetres on the page, measured from its top-left
# corner.  The canvas of the figure is the page minus the pad that style.save
# adds on every side, so that the saved file is exactly PAGE_W wide.
PAGE_W = style.DOUBLE / style.CM       # 17.1 cm
PAD = 0.02 / style.CM                  # pad_inches of style.save, in cm
FIG_W = PAGE_W - 2 * PAD

TOP = 0.95                             # top edge of the first row of (a), (b)
GAP = 0.12                             # gap between stacked axes
LABEL_DX = -0.10                       # panel label: right end from axes edge

# (a) left column: three views of the active site above four stacked axes
A_L, A_W = 1.25, 6.35
# view direction of the site: azimuth about the surface normal and elevation
# above the surface plane (degrees), chosen so that the drawn atoms overlap
# least in the three images
VIEW_AZ, VIEW_TILT = 110.0, 15.0
SITE_U = 7.3                           # width of one view, in A
SITE_V = (-2.75, 5.35)                 # vertical range of a view, in A
A_SITE_H = (SITE_V[1] - SITE_V[0]) * (A_W / 3) / SITE_U
A_NAMES = 0.42                         # row of state names under the views
A_H = (2.25, 1.35, 2.00, 1.60)         # Co4, Ce(+3), H/C and distance axes
A_TOPS = [TOP + A_SITE_H + A_NAMES]
for _h in A_H[:-1]:
    A_TOPS.append(A_TOPS[-1] + _h + GAP)
A_BOTTOM = A_TOPS[-1] + A_H[-1]

# (b), (c) right column
B_L, B_W = 9.55, 7.35
B_H = (2.55, 2.30, 1.40)               # metal, H atoms, H-H distance
B_TOPS = [TOP, TOP + B_H[0] + GAP, TOP + B_H[0] + B_H[1] + 2 * GAP]
C_H = 2.50                             # (c): height of axes and sketch
C_W = 5.30                             # (c): width of the axes
C_TOP = A_BOTTOM - C_H

PAGE_H = A_BOTTOM + 0.54               # room for one line of tick labels
FIG_H = PAGE_H - 2 * PAD

FLOW, GREEDY, REFC, OCCC = C["flow"], C[BASE], C["ref"], C["muted"]
# Marker of the diagnostic.  Deviation from the one-colour-one-marker mapping
# of style.py, on purpose: the marker is MARKER["aux2"] but the colour is
# C["muted"], not C["aux2"].  C["aux2"] (#E69F00) is an orange close to the
# orange of the greedy rule (#D55E00), with which the plus signs share panel
# (c); grey also sets the diagnostic apart from the two methods.
OCC_MK, OCC_MS = MARKER["aux2"], 4.4
OCC_LABEL = "WF-occupation-weighted"
LIT_LABEL = "Comas-Vilà et al."
ONE_SITE_LABEL = "cell with one H site"
HT = r"H$_{\mathrm{t}}$"               # the transferred H atom
CE3_ATOMS = r"Ce$^{3{+}}\!$ atoms"


def ax_cm(fig, left, top, width, height, **kw):
    """Axes placed in centimetres measured from the top-left corner."""
    return fig.add_axes([(left - PAD) / FIG_W,
                         1 - (top - PAD + height) / FIG_H,
                         width / FIG_W, height / FIG_H], **kw)


def label_panel(ax, letter, x_cm=LABEL_DX, y_cm=0.10):
    """style.panel_label with its right end x_cm from the left axes edge
    (negative = to the left) and its baseline y_cm above the top edge."""
    w, h = (ax.get_position().width * FIG_W, ax.get_position().height * FIG_H)
    style.panel_label(ax, letter, dx=x_cm / w, dy=1 + y_cm / h)


def ylabel_cm(fig, axes_left, top, height, text, offset=0.80):
    """Rotated y label centred on an axes (or a block of stacked axes)."""
    fig.text((axes_left - offset - PAD) / FIG_W,
             1 - (top - PAD + height / 2) / FIG_H,
             text, rotation=90, ha="center", va="center", fontsize=8)


def signed(v, _pos=None):
    v = int(round(v))
    return "0" if v == 0 else (f"+{v}" if v > 0 else f"−{-v}")


# --- view of the active site -----------------------------------------------
# element: (drawn radius in A, face colour); the transferred H is black
ATOM = {
    "Co": (0.62, C["muted"]),
    "O": (0.46, "white"),
    "C": (0.44, C["ref"]),
    "H": (0.27, "white"),
}
BOND_CUT = {("C", "H"): 1.30, ("Co", "H"): 2.05, ("C", "Co"): 2.30,
            ("Co", "O"): 2.05, ("Co", "Co"): 2.90}
FORMING = 2.60      # a C-H_t contact between 1.30 A and this is drawn dotted


def site_atoms():
    """The cluster, the molecule and the O atoms bound to the cluster
    (Co-O below 2.05 A) in all three images."""
    keep = None
    for d in CH4_SETS["gs"].values():
        st = parse_poscar(d / "POSCAR")
        el = st.elements
        co = idx(el, "Co")
        near = {i for i in idx(el, "O")
                if min(min_image_dist(st, i, k) for k in co)
                < BOND_CUT[("Co", "O")]}
        keep = near if keep is None else keep & near
    el = parse_poscar(CH4_SETS["gs"]["is"] / "POSCAR").elements
    return sorted(keep) + [i for i, e in enumerate(el)
                           if e in ("Co", "C", "H")]


def site_projection(poscar, atoms):
    """Positions of the atoms relative to the centroid of the Co4 cluster,
    seen from VIEW_TILT degrees above the surface plane at azimuth VIEW_AZ:
    {atom: (u, v, depth)}, depth increasing away from the viewer."""
    st = parse_poscar(poscar)
    co = idx(st.elements, "Co")
    lat = np.array(st.lattice)
    cen = np.array(st.frac[co[0]]) @ lat + np.mean(
        [min_image_vec(st, st.frac[co[0]], st.frac[i]) for i in co], axis=0)
    cen_frac = cen @ np.linalg.inv(lat)
    a, t = np.radians(VIEW_AZ), np.radians(VIEW_TILT)
    proj = {}
    for i in atoms:
        x, y, z = min_image_vec(st, cen_frac, st.frac[i])
        xr, yr = x * np.cos(a) + y * np.sin(a), -x * np.sin(a) + y * np.cos(a)
        proj[i] = (xr, yr * np.sin(t) + z * np.cos(t),
                   yr * np.cos(t) - z * np.sin(t))
    return st, proj


def draw_site(ax, poscar, atoms, h_t):
    st, proj = site_projection(poscar, atoms)
    el = st.elements
    radii = {i: ATOM[el[i]][0] for i in atoms}
    order = sorted(atoms, key=lambda i: -proj[i][2])    # far atoms first
    rank = {i: n for n, i in enumerate(order)}
    c_atom = idx(el, "C")[0]
    for i, j in itertools.combinations(atoms, 2):
        d = min_image_dist(st, i, j)
        cut = BOND_CUT.get(tuple(sorted((el[i], el[j]))))
        forming = {i, j} == {c_atom, h_t} and BOND_CUT[("C", "H")] <= d \
            < FORMING
        if not ((cut and d < cut) or forming):
            continue
        # stick between the two discs (nothing when they overlap)
        p, q = np.array(proj[i][:2]), np.array(proj[j][:2])
        length = float(np.linalg.norm(q - p))
        if length <= radii[i] + radii[j]:
            continue
        u = (q - p) / length
        p, q = p + u * radii[i], q - u * radii[j]
        ax.plot([p[0], q[0]], [p[1], q[1]], color=C["ink"], lw=0.6,
                ls=(0, (1, 1.5)) if forming else "-",
                solid_capstyle="butt", zorder=2 + 0.5 * (rank[i] + rank[j]))
    for i in order:
        r, face = ATOM[el[i]]
        ax.add_patch(Circle(proj[i][:2], r, zorder=2 + rank[i],
                            facecolor=C["ink"] if i == h_t else face,
                            edgecolor=C["ink"], lw=0.5))
    ax.set_xlim(-SITE_U / 2, SITE_U / 2)
    ax.set_ylim(*SITE_V)
    ax.set_aspect("equal")
    ax.axis("off")


def draw_atom_key(ax):
    """Element key of the site views, drawn to the scale of the views in
    the free strip above the second and third view (right-aligned)."""
    pt = 72 * (A_W / 3) / SITE_U * style.CM       # points per A
    entries = [(el, ATOM[el][0], ATOM[el][1], el) for el in
               ("Co", "O", "C", "H")] + [("Ht", ATOM["H"][0], C["ink"], HT)]
    # approximate widths of the labels at 7 pt, in points (spacing only)
    text_w = {"Co": 9.0, "O": 5.5, "C": 5.1, "H": 5.1, "Ht": 7.4}
    pad, gap = 0.16, 0.55                         # A
    total = sum(2 * r + pad + text_w[k] / pt for k, r, _, _ in entries) \
        + gap * (len(entries) - 1)
    x = 2 * SITE_U - total - 0.25
    y = SITE_V[1] - 0.78
    for k, r, face, text in entries:
        ax.add_patch(Circle((x + r, y), r, facecolor=face,
                            edgecolor=C["ink"], lw=0.5))
        ax.text(x + 2 * r + pad, y, text, ha="left", va="center_baseline",
                fontsize=7)
        x += 2 * r + pad + text_w[k] / pt + gap
    ax.set_xlim(0, 2 * SITE_U)
    ax.set_ylim(*SITE_V)
    ax.set_aspect("equal")
    ax.axis("off")


# --- panel (a) ---------------------------------------------------------------
def draw_a(fig, ch4):
    atoms = site_atoms()
    h_t = ch4["h_atom"] - 1
    w3 = A_W / 3
    for n, (image, _) in enumerate(IMAGES):
        ax = ax_cm(fig, A_L + n * w3, TOP, w3, A_SITE_H)
        draw_site(ax, CH4_SETS["gs"][image] / "POSCAR", atoms, h_t)
        if n == 0:
            label_panel(ax, "a")
    draw_atom_key(ax_cm(fig, A_L + w3, TOP, 2 * w3, A_SITE_H))
    ax_co = ax_cm(fig, A_L, A_TOPS[0], A_W, A_H[0])
    ax_ce = ax_cm(fig, A_L, A_TOPS[1], A_W, A_H[1], sharex=ax_co)
    ax_h = ax_cm(fig, A_L, A_TOPS[2], A_W, A_H[2], sharex=ax_co)
    ax_d = ax_cm(fig, A_L, A_TOPS[3], A_W, A_H[3], sharex=ax_co)
    x = np.arange(len(IMAGES))

    def flow_kw(filled):
        return dict(color=FLOW, marker=MARKER["flow"],
                    lw=1.0 if filled else 0.8,
                    ls="-" if filled else (0, (1, 1.5)),
                    ms=3.6 if filled else 5.4,
                    mfc=FLOW if filled else "white",
                    mew=0 if filled else 0.8, zorder=4 if filled else 3)

    def occ_kw(filled):
        return dict(color=OCCC, marker=OCC_MK, lw=0, ms=OCC_MS,
                    mfc=OCCC if filled else "white", mew=0.7, mec=OCCC,
                    zorder=5)

    for s, filled in (("legacy", False), ("gs", True)):
        rows = ch4["sets"][s]
        dx = 0.075 if filled else -0.075
        ax_co.plot(x, [r["co4_flow"] for r in rows], **flow_kw(filled))
        ax_ce.plot(x, [r["n_ce3_flow"] for r in rows], **flow_kw(filled))
        if SHOW_OCC:
            ax_co.plot(x + dx, [r["co4_occ"] for r in rows], **occ_kw(filled))
            ax_ce.plot(x + dx, [r["n_ce3_occ"] for r in rows],
                       **occ_kw(filled))
        ax_h.plot(x, [r["h_flow"] for r in rows], **flow_kw(filled))
        ax_h.plot(x, [r["c_flow"] for r in rows], **flow_kw(filled))

    ax_co.set_ylim(-1.85, 5.95)
    ax_co.yaxis.set_major_locator(FixedLocator([0, 2, 4]))
    ax_co.yaxis.set_minor_locator(FixedLocator([-1, 1, 3, 5]))
    ax_co.yaxis.set_major_formatter(FuncFormatter(signed))
    ax_co.text(0.03, 0.93, formula("Co4") + " cluster",
               transform=ax_co.transAxes, ha="left", va="top")
    handles = [Line2D([], [], color=C["ink"], lw=0, marker="o", ms=3.6,
                      mew=0, label=CH4_SET_LABEL["gs"]),
               Line2D([], [], color=C["ink"], lw=0, marker="o", ms=5.4,
                      mfc="none", mew=0.8, label=CH4_SET_LABEL["legacy"])]
    ax_co.legend(handles=handles, loc="lower right",
                 bbox_to_anchor=(1.0, 0.0), borderaxespad=0.45,
                 handletextpad=0.3, labelspacing=0.3)

    ax_ce.set_ylim(-0.6, 3.6)
    ax_ce.yaxis.set_major_locator(MultipleLocator(1))
    ax_ce.yaxis.set_minor_locator(NullLocator())
    ax_ce.text(0.03, 1.0, CE3_ATOMS,
               transform=ax_ce.get_yaxis_transform(), ha="left", va="center")

    ax_h.set_ylim(-4.95, 1.95)
    ax_h.yaxis.set_major_locator(FixedLocator([-4, -1, 1]))
    ax_h.yaxis.set_minor_locator(MultipleLocator(1))
    ax_h.yaxis.set_major_formatter(FuncFormatter(signed))
    ax_h.annotate(HT, (1.5, -1), xytext=(0, 4), textcoords="offset points",
                  ha="center", va="bottom")
    ax_h.annotate("C", (1.5, -4), xytext=(0, 4), textcoords="offset points",
                  ha="center", va="bottom")

    rows = ch4["sets"]["gs"]
    for key, ls, text, dy, va in (("d_ch", "-", "C–" + HT, 4, "bottom"),
                                  ("d_coh", (0, (4, 2)), "Co–" + HT, -4,
                                   "top")):
        y = [r[key] for r in rows]
        ax_d.plot(x, y, color=C["ink"], marker="o", ms=3.0, mew=0, lw=0.9,
                  ls=ls, zorder=3)
        ax_d.annotate(text, (1.5, 0.5 * (y[1] + y[2])), xytext=(0, dy),
                      textcoords="offset points", ha="center", va=va)
    ax_d.set_ylim(0.3, 3.9)
    ax_d.yaxis.set_major_locator(MultipleLocator(1))
    ax_d.yaxis.set_minor_locator(MultipleLocator(0.5))

    ax_co.set_xlim(-0.5, 2.5)
    ax_co.set_xticks(x)
    ax_co.set_xticklabels([lab for _, lab in IMAGES])
    for ax in (ax_co, ax_ce, ax_h, ax_d):
        ax.xaxis.set_minor_locator(NullLocator())
        ax.tick_params(labelbottom=False)
    # state names under the views (top of the stack) and under the last axes
    ax_co.tick_params(labeltop=True, pad=3)
    ax_d.tick_params(labelbottom=True)
    for top, h, text in zip(A_TOPS, A_H, (
            "Oxidation state", "Number", "Oxidation state",
            "Distance (Å)")):
        ylabel_cm(fig, A_L, top, h, text)


# --- panel (b) ---------------------------------------------------------------
def draw_b(fig, rows):
    ax_m = ax_cm(fig, B_L, B_TOPS[0], B_W, B_H[0])
    ax_h = ax_cm(fig, B_L, B_TOPS[1], B_W, B_H[1], sharex=ax_m)
    ax_d = ax_cm(fig, B_L, B_TOPS[2], B_W, B_H[2], sharex=ax_m)
    lit_kw = dict(color=REFC, marker=MARKER["ref"], ms=6.0, mfc="none",
                  mew=0.8, lw=0, zorder=2)
    flow_kw = dict(color=FLOW, marker=MARKER["flow"], ms=3.4, mew=0, lw=0,
                   zorder=4)
    xb, xa, dx = -0.22, 0.22, 0.075       # H_b and H_a slots, atom offset
    for n, r in enumerate(rows):
        if not r["two_sites"]:         # cell is not the literature structure
            for ax in (ax_m, ax_h, ax_d):
                ax.axvspan(n - 0.5, n + 0.5, color=C["light"], alpha=0.45,
                           lw=0, zorder=0)
        # metal
        ax_m.plot([n], [r["lit"][0]], **lit_kw)
        ax_m.plot([n] * len(set(r["m_flow"])), sorted(set(r["m_flow"])),
                  **flow_kw)
        if set(r["m_greedy"]) != set(r["m_flow"]):
            g = sorted(set(r["m_greedy"]))
            # one square per distinct value among the metal atoms of the cell
            # (ThH4: +8 and +12 for its two Th atoms); not joined, because a
            # joining line reads as a range or an error bar
            ax_m.plot([n] * len(g), g, color=GREEDY, marker=MARKER[BASE],
                      ms=3.8, mew=0, lw=0, zorder=3)
        # hydrogen: literature per site, flow per atom of one formula unit
        ax_h.plot([n + xb, n + xa], [r["lit"][1], r["lit"][2]], **lit_kw)
        if r["two_sites"]:
            ax_h.plot([n + xb - dx, n + xb + dx], r["hb"], **flow_kw)
            ax_h.plot([n + xa - dx, n + xa + dx], r["ha"], color=FLOW,
                      lw=0.8 if r["h2"] else 0, marker=MARKER["flow"],
                      ms=3.4, mew=0, zorder=4)
        else:
            ax_h.plot(n + np.array([xb - dx, xb + dx, xa - dx, xa + dx]),
                      r["h_all"], **flow_kw)
        ax_d.plot([n], [r["d_min"]], color=C["ink"], marker="o", ms=3.4,
                  mew=0, lw=0, zorder=3)
    ax_m.set_ylim(0.6, 13.4)
    ax_m.yaxis.set_major_locator(MultipleLocator(2))
    ax_m.yaxis.set_minor_locator(MultipleLocator(1))
    ax_m.yaxis.set_major_formatter(FuncFormatter(signed))
    ax_m.text(0.015, 0.93, "Metal", transform=ax_m.transAxes, ha="left",
              va="top")
    ax_h.set_ylim(-1.6, 1.95)
    ax_h.yaxis.set_major_locator(FixedLocator([-1, 0, 1]))
    ax_h.yaxis.set_minor_locator(NullLocator())
    ax_h.yaxis.set_major_formatter(FuncFormatter(signed))
    ax_h.text(0.015, 0.93, "H atoms", transform=ax_h.transAxes, ha="left",
              va="top")
    # name the two H sites once, above a compound in which both are at -1
    k = next(k for k, r in enumerate(rows) if r["two_sites"]
             and r["lit"][2] == -1 and set(r["ha"]) == {-1})
    for off, name in ((xb, "H$_b$"), (xa, "H$_a$")):
        ax_h.annotate(name, (k + off, -1), xytext=(0, 6),
                      textcoords="offset points", ha="center", va="bottom",
                      fontsize=7)
    ax_d.set_ylim(0.35, 2.25)
    ax_d.yaxis.set_major_locator(MultipleLocator(1))
    ax_d.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax_m.set_xlim(-0.55, len(rows) - 0.45)
    ax_m.set_xticks(range(len(rows)))
    ax_d.set_xticklabels([formula(r["metal"] + "H4")
                          + ("$^{\\dagger}$" if r["stray"] else "")
                          + ("$^{\\ddagger}$" if r["lost"] else "")
                          for r in rows])
    for ax in (ax_m, ax_h, ax_d):
        ax.xaxis.set_minor_locator(NullLocator())
        ax.tick_params(axis="x", which="both", top=False)
    ax_m.tick_params(labelbottom=False)
    ax_h.tick_params(labelbottom=False)
    ax_d.tick_params(axis="x", pad=10.5)
    for t in ax_d.get_xticklabels():
        t.set_verticalalignment("baseline")
    ylabel_cm(fig, B_L, B_TOPS[0], B_H[0] + GAP + B_H[1], "Oxidation state")
    ylabel_cm(fig, B_L, B_TOPS[2], B_H[2], r"$d_{\mathrm{H–H}}$ (Å)")
    label_panel(ax_m, "b")


# --- panel (c) ---------------------------------------------------------------
def draw_c(fig, ce):
    ax = ax_cm(fig, B_L, C_TOP, C_W, C_H)
    xs = {"ce4o8": np.arange(4), "ce4o7": np.arange(4) + 4.6}
    off = 0.29
    for name, lab in CEO2:
        d = ce[name]
        x = xs[name]
        if SHOW_OCC:
            ax.plot(x, d["ce_os_occ"], color=OCCC, marker=OCC_MK, ms=OCC_MS,
                    mew=0, lw=0, zorder=5)
        ax.plot(x + off, d["ce_greedy"], color=GREEDY, marker=MARKER[BASE],
                ms=3.8, mew=0, lw=0, zorder=3)
        ax.plot(x - off, d["ce_flow"], color=FLOW, marker=MARKER["flow"],
                ms=3.6, mew=0, lw=0, zorder=4)
        ax.text(float(np.mean(x)), 4.42, formula(lab), ha="center",
                va="center")
    ax.axvline(3.8, color=C["ink"], lw=0.6)
    ax.set_xlim(-0.6, 8.2)
    ax.set_xticks(list(xs["ce4o8"]) + list(xs["ce4o7"]))
    ax.set_xticklabels([f"Ce{n}" for n in ce["ce4o8"]["ce_atoms"]]
                       + [f"Ce{n}" for n in ce["ce4o7"]["ce_atoms"]])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.tick_params(axis="x", which="both", top=False)
    ax.set_ylim(2.72, 4.68)
    ax.yaxis.set_major_locator(FixedLocator([3, 3.5, 4]))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"+{v:g}"))
    ylabel_cm(fig, B_L, C_TOP, C_H, "Oxidation state")
    label_panel(ax, "c")

    # sketch: the vacancy and its four Ce neighbours, seen along [001]
    left = B_L + C_W + 0.20
    width = B_L + B_W - left
    scale = 0.40                                   # cm per A
    ax_s = ax_cm(fig, left, C_TOP, width, C_H)
    v = ce["vacancy"]
    for n, vec in enumerate(v["vec_ce"]):
        ax_s.plot([0, vec[0]], [0, vec[1]], color=C["ink"], lw=0.6,
                  ls=(0, (2, 1.5)), zorder=1)
        ax_s.add_patch(Circle(vec[:2], 0.70, facecolor=C["light"],
                              edgecolor=C["ink"], lw=0.5, zorder=2))
        ax_s.text(vec[0], vec[1], f"Ce{n + 1}", ha="center",
                  va="center_baseline", fontsize=7, zorder=3)
    ax_s.add_patch(Circle((0, 0), 0.60, facecolor="white",
                          edgecolor=C["ink"], lw=0.6, ls=(0, (2, 1.5)),
                          zorder=2))
    ax_s.text(0, 0, r"V$_{\mathrm{O}}$", ha="center", va="center_baseline",
              fontsize=7, zorder=3)
    d = sorted(set(round(x, 2) for x in v["d_ce"]))
    assert len(d) == 1, d
    ymax = 0.5 * C_H / scale
    ax_s.text(0, -ymax + 0.10, f"4 × {d[0]:.2f} Å", ha="center",
              va="bottom", fontsize=7)
    ax_s.set_xlim(-0.5 * width / scale, 0.5 * width / scale)
    ax_s.set_ylim(-ymax - 0.25, ymax - 0.25)
    ax_s.set_aspect("equal")
    ax_s.axis("off")


def method_legend(fig):
    handles = [
        Line2D([], [], color=FLOW, marker=MARKER["flow"], ms=3.6, mew=0,
               lw=0, label=LABEL["flow"]),
        Line2D([], [], color=GREEDY, marker=MARKER[BASE], ms=3.8, mew=0,
               lw=0, label=LABEL[BASE]),
    ]
    if SHOW_OCC:
        handles.append(Line2D([], [], color=OCCC, marker=OCC_MK, ms=OCC_MS,
                              mew=0, lw=0, label=OCC_LABEL))
    handles += [
        Line2D([], [], color=REFC, marker=MARKER["ref"], ms=6.0, mfc="none",
               mew=0.8, lw=0, label=LIT_LABEL),
        Patch(facecolor=C["light"], alpha=0.45, lw=0, label=ONE_SITE_LABEL),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=len(handles),
               bbox_to_anchor=(0.5, 1.0), columnspacing=1.6,
               handletextpad=0.4, borderaxespad=0.15, handlelength=1.3,
               handleheight=0.8)


def print_provenance():
    w = [max(len(r[i]) for r in PROV) for i in range(3)]
    w[1], w[2] = min(w[1], 62), min(w[2], 66)
    print("\nPROVENANCE OF THE PLOTTED VALUES (panel | quantity | value | "
          "source)")
    print("-" * 150)
    for p, q, v, s in PROV:
        print(f"({p}) {q:<{w[1]}}  {v:<{w[2]}}  {s}")
    print("-" * 150)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--recompute", action="store_true",
                    help="ignore _cache/fig5_*.json and rerun pyleos")
    args = ap.parse_args()

    ch4 = load_methane(args.recompute)
    mh4 = load_hydrides(args.recompute)
    ce = load_ceo2(args.recompute)
    report_methane(ch4)
    report_hydrides(mh4)
    report_ceo2(ce)

    fig = plt.figure(figsize=(FIG_W * style.CM, FIG_H * style.CM))
    # invisible frame: pins the tight bounding box of style.save to the canvas
    fig.add_artist(Rectangle((0, 0), 1, 1, transform=fig.transFigure,
                             fill=False, edgecolor="none", lw=0))
    method_legend(fig)
    draw_a(fig, ch4)
    draw_b(fig, mh4)
    draw_c(fig, ce)
    style.save(fig, NAME)
    print_provenance()
    print(f"wrote {NAME}.pdf / {NAME}.png (page {PAGE_W:.1f} x "
          f"{PAGE_H:.2f} cm)")


if __name__ == "__main__":
    main()
