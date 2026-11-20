#!/usr/bin/env python3
"""Fig. 3 - comparison with the atomic-overlap-matrix method on 40 solids.

Run from this directory:

    python3 figS1_comparison.py              # uses _cache/fig3_leos.json if valid
    python3 figS1_comparison.py --recompute  # re-run pyleos on every case

Panels
  (a) parity plot of element-mean oxidation states, Comas-Vila et al. (x)
      versus this work (y), one point per element per compound;
  (b) Wannier-function spread versus distance of the WF centre from the
      nearest nucleus for two electrides and two ionic references;
  (c) number of solids per chemical class that agree / differ / cannot be
      compared.

Every LEOS number is computed here from the run directories under cases/
with the pyleos defaults (beta = 10, Allen electronegativity, plain
geometric + electronegativity cost).  Literature values are typed once, in
the table LIT below.  Output: figS1_comparison.pdf / .png (style.save) and a
provenance table on stdout.

Needs numpy + matplotlib; pyleos itself is standard-library only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import defaultdict
from fractions import Fraction

import numpy as np

import style
from style import C, CM, MARKER, formula, panel_label

import matplotlib.pyplot as plt                              # noqa: E402
from matplotlib.lines import Line2D                          # noqa: E402
from matplotlib.patches import Patch, Rectangle              # noqa: E402
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter  # noqa: E402

sys.path.insert(0, str(style.REPO))
from pyleos.cases import load_index                         # noqa: E402
from pyleos.parsers import load_case, parse_poscar           # noqa: E402
from pyleos.transport import STRAY_SPREAD, assign_case       # noqa: E402

DATA = style.REPO / "cases"                       # <case>/<run>/ directories
GOLDEN_TSV = style.REPO / "results" / "bench74_golden.tsv"
CACHE = style.HERE / "_cache" / "fig3_leos.json"
# role / bench_set / ref_key of every run directory, keyed "<case>/<run>"
INDEX = {f"{r['case']}/{r['run']}": r for r in load_index(style.REPO)}
NAME = "figS1_comparison"

# =============================================================================
# LITERATURE VALUES - the only hand-typed numbers in this script.
#
# Source: G. Comas-Vila et al., J. Chem. Theory Comput. 2025, 21, 7075-7086,
# as transcribed in docs/comas_vila_2025_notes.md, table "Comparison by class".
# "N<k>" is line k of that file; "p." is the journal page given there (each
# value was also re-read in docs/refs/Comas-Vilà ... 2025 ... .pdf).
#
# One row per solid:
#   id        number of the solid in Table 1 of the paper
#   dir       LEOS run directory under cases/
#   formula, tag   label (tag distinguishes polymorphs)
#   cls       chemical class used in panel (c)
#   struct    how the LEOS cell relates to the structure used in the paper
#               same      same phase and cell content
#               lattice   same phase, lattice parameters approximate
#               different different cell content or phase
#               unmatched correspondence with the paper's cell not established
#   lit       {species: element-mean oxidation state}; "a..b" = the paper gives
#             site-resolved values only and the element mean lies in (a, b);
#             {} = the paper reports no assignment
#   basis     how each element mean follows from what the paper states
# Species "OH" is the hydroxide fragment (the paper assigns OH as one unit).
# =============================================================================
LIT = [
    # id  dir                      formula      tag        cls         struct       lit                                                N    p.           basis
    (1,  "cro3/v3",                "CrO3",      "",        "salt",     "same",      {"Cr": "6", "O": "-2"},                            341, "7078-7079", "Cr stated; O by charge neutrality"),
    (2,  "kmgf3/v3",               "KMgF3",     "",        "salt",     "same",      {"K": "1", "Mg": "2", "F": "-1"},                  342, "7079",      "K+ Mg2+ (F-)3 as printed"),
    (3,  "nabf4/v3",               "NaBF4",     "",        "salt",     "same",      {"Na": "1", "B": "3", "F": "-1"},                  342, "7079",      "Na+ B3+ (F-)4 as printed"),
    (4,  "srbe_oh_4/v3",           "SrBe(OH)4", "",        "salt",     "same",      {"Sr": "2", "Be": "2", "OH": "-1"},                342, "7079",      "Sr2+ Be2+ (OH-)4 as printed; OH treated as a fragment"),
    (5,  "cs3auo/v3",              "Cs3AuO",    "",        "auride",   "same",      {"Cs": "1", "Au": "-1", "O": "-2"},                343, "7079",      "stated consistent with the IUPAC report (Cs+ Au- O2-)"),
    (6,  "rb3auo/v3",              "Rb3AuO",    "",        "auride",   "same",      {"Rb": "1", "Au": "-1", "O": "-2"},                343, "7079",      "stated consistent with the IUPAC report (Rb+ Au- O2-)"),
    (7,  "ca3aun/v3",              "Ca3AuN",    "",        "electride", "same",     {"Ca": "2", "Au": "-3", "N": "-3"},                344, "7079",      "Au(-3) stated; Ca valence eff-AO unoccupied, i.e. Ca(+2); N(-3) by charge neutrality"),
    (8,  "k_hp/v3",                "K",         "",        "electride", "different", {"K": "1"},                                       345, "7079",      "K(+1) stated; 2 e- on each non-nuclear attractor (no LEOS counterpart)"),
    (9,  "sega/v3",                "GaSe",      "",        "salt",     "same",      {"Ga": "2", "Se": "-2"},                           346, "7079-7080", "Ga(+2) Se(-2) stated (Ga-Ga unit +4)"),
    (10, "ba3si4/v3",              "Ba3Si4",    "",        "zintl",    "same",      {"Ba": "2", "Si": "-3/2"},                         347, "7080",      "Ba(+2) stated; Si4 fragment -6, i.e. -6/4 per Si"),
    (11, "cs2pt/v3",               "Cs2Pt",     "",        "platinide", "lattice",  {"Cs": "1", "Pt": "-2"},                           348, "7080",      "Pt(-2) stated; Cs by charge neutrality"),
    (12, "wcl4/v3",                "WCl4",      "",        "salt",     "same",      {"W": "4", "Cl": "-1"},                            349, "7080",      "W(4+) Cl(1-) stated"),
    (13, "rb2cucl4/v3",            "Rb2CuCl4",  "",        "salt",     "same",      {"Rb": "1", "Cu": "2", "Cl": "-1"},                350, "7080-7081", "Cu(+2) stated; Cl(-1), Rb(+1) as in the text"),
    (14, "ybafe2o5_a/v3",          "YBaFe2O5",  "ordered", "mixed",    "same",      {"Y": "3", "Ba": "2", "Fe": "5/2", "O": "-2"},     351, "7081",      "Fe1(+3), Fe2(+2) stated, 1:1; Y, Ba, O as in the IUPAC report"),
    (15, "ybafe2o5_b/v3",          "YBaFe2O5",  "mixed",   "mixed",    "same",      {"Y": "3", "Ba": "2", "Fe": "5/2", "O": "-2"},     351, "7081",      "Fe(+2.5) stated; Y, Ba, O as in the IUPAC report"),
    (16, "mgh4/v3",                "MgH4",      "",        "mh4",      "same",      {"Mg": "2", "H": "-1/2"},                          352, "7081",      "M(+2), Hb(-1), neutral Ha-Ha; cell = 2 M + 4 Hb + 4 Ha"),
    (17, "cah4/v3",                "CaH4",      "",        "mh4",      "same",      {"Ca": "2", "H": "-1/2"},                          352, "7081",      "as MgH4"),
    (18, "srh4/v3",                "SrH4",      "",        "mh4",      "same",      {"Sr": "2", "H": "-1/2"},                          352, "7081",      "as MgH4"),
    (19, "sch4/v3",                "ScH4",      "",        "mh4",      "same",      {"Sc": "3", "H": "-3/4"},                          353, "7081-7082", "M(+3), Hb(-1), Ha(-1/2); cell = 2 M + 4 Hb + 4 Ha"),
    (20, "yh4/v3",                 "YH4",       "",        "mh4",      "same",      {"Y": "3", "H": "-3/4"},                           353, "7081-7082", "as ScH4"),
    (21, "zrh4/v3",                "ZrH4",      "",        "mh4",      "same",      {"Zr": "4", "H": "-1"},                            354, "7082",      "M4+ (H-)4 stated"),
    (22, "ceh4/v3",                "CeH4",      "",        "mh4",      "same",      {"Ce": "4", "H": "-1"},                            354, "7082",      "M4+ (H-)4 stated"),
    (23, "thh4/v3",                "ThH4",      "",        "mh4",      "same",      {"Th": "4", "H": "-1"},                            354, "7082",      "M4+ (H-)4 stated"),
    (24, "xeo_a/v3",               "XeO",       "100 GPa", "xeo",      "same",      {"Xe": "2", "O": "-2"},                            355, "7082",      "Xe(2+) O(2-) stated"),
    (25, "xeo_b/v3",               "XeO",       "200 GPa", "xeo",      "same",      {"Xe": "2", "O": "-2"},                            355, "7082",      "Xe1(+4), Xe2(0), O(-2); 1:1 by charge neutrality"),
    (26, "xeo2_a/v3",              "XeO2",      "150 GPa", "xeo",      "same",      {"Xe": "4", "O": "-2"},                            356, "7082",      "Xe(+4) O(-2)2 stated"),
    (27, "xeo2_b/v3",              "XeO2",      "200 GPa", "xeo",      "same",      {"Xe": "4", "O": "-2"},                            356, "7082",      "Xe(+4) O(-2)2 stated"),
    (28, "xeo3_a/v3",              "XeO3",      "130 GPa", "xeo",      "same",      {"Xe": "4", "O": "-4/3"},                          357, "7082-7083", "Xe(+4), O1(-2) x2, O2(0) x1 per formula unit (atomic partition)"),
    (29, "xeo3_b/v3",              "XeO3",      "150 GPa", "xeo",      "same",      {"Xe": "4..6", "O": "-2..-4/3"},                   358, "7083",      "Xe1(+6), Xe2(+4) stated without multiplicities; O interval by charge neutrality"),
    (30, "xef2/v3",                "XeF2",      "",        "xef",      "same",      {"Xe": "2", "F": "-1"},                            359, "7083",      "Xe(+2) stated; F by charge neutrality"),
    (31, "xef3/v3",                "XeF3",      "",        "xef",      "same",      {"Xe": "3", "F": "-1"},                            359, "7083",      "Xe1(+2), Xe2(+4) in XeF2.XeF4, 1:1; F by charge neutrality"),
    (32, "xef4/v3",                "XeF4",      "",        "xef",      "same",      {"Xe": "4", "F": "-1"},                            359, "7083",      "Xe(4+) stated; F by charge neutrality"),
    (33, "mgxe/v3",                "MgXe",      "",        "ngan",     "lattice",   {"Mg": "2", "Xe": "-2"},                           360, "7083",      "Xe(2-) Mg(2+) stated"),
    (34, "hgf4_a/v3",              "HgF4",      "cell a",  "hg4",      "unmatched", {"Hg": "4", "F": "-1"},                            361, "7083",      "Hg(+4) stated for all three structures; F by charge neutrality"),
    (35, "hgf4_b/v3",              "HgF4",      "cell b",  "hg4",      "unmatched", {"Hg": "4", "F": "-1"},                            361, "7083",      "as cell a"),
    (36, "hgf4_c/v3",              "HgF4",      "cell c",  "hg4",      "unmatched", {"Hg": "4", "F": "-1"},                            361, "7083",      "as cell a"),
    (37, "csh3/v3",                "CsH3",      "",        "poly",     "different", {"Cs": "1", "H": "-1/3"},                          362, "7083",      "H3 fragment -1 (atomic partition: -1, +1, -1)"),
    (38, "csf3/v3",                "CsF3",      "",        "poly",     "same",      {"Cs": "1", "F": "-1/3"},                          363, "7083",      "F3 fragment -1"),
    (39, "csf5/v3",                "CsF5",      "",        "poly",     "same",      {"Cs": "1", "F": "-1/5"},                          364, "7084",      "Cs(+1) F5(-1) stated"),
    (40, "lin5/v3",                "LiN5",      "",        "poly",     "different", {},                                                365, "7084",      "no assignment: R = 50 %, declared completely ambiguous"),
]

# Why a structure is not flagged "same" (shown in the provenance table).
STRUCT_NOTE = {
    "k_hp/v3": "LEOS cell is oP8 K (8 atoms); the paper's cell has 4 K "
               "atoms and 2 non-nuclear attractors (p. 7079)",
    "cs2pt/v3": "Ni2In-type cell, but c/a is 16.6 % below experiment "
                "(docs/testset_registry.tsv, SUSPECT)",
    "mgxe/v3": "Pm-3m cell with lattice parameter inferred from an "
               "equation of state, +-1.5 % (docs/RESULTS.md, sec. 4)",
    "hgf4_a/v3": "paper: three I4/m structures (50, 100, 200 GPa); "
                 "this cell: square-planar HgF4 molecules, not tetragonal",
    "hgf4_b/v3": "paper: I4/m; this cell: C-centred monoclinic with "
                 "F2 units (F-F 1.43 A)",
    "hgf4_c/v3": "paper: I4/m; this cell: cubic, 5 atoms",
    "csh3/v3": "paper: linear H3 units; this cell: two isolated "
               "triangular H3 (docs/testset_findings.md)",
    "lin5/v3": "paper: pentazolate N5 ring; this cell: N3 + N2 "
               "(docs/testset_findings.md)",
}

CLASSES = [                     # (key, label) in the order of the paper's IDs
    ("salt", "Oxides, halides, GaSe"),
    ("auride", "Aurides"),
    ("electride", "Electrides"),
    ("zintl", "Zintl phase"),
    ("platinide", "Platinide"),
    ("mixed", "Mixed valence"),
    ("mh4", "MH$_4$ hydrides"),
    ("xeo", "High-pressure Xe oxides"),
    ("xef", "Xe fluorides"),
    ("ngan", "Noble-gas anion"),
    ("hg4", "Hg(IV) fluoride"),
    ("poly", "Polyatomic anions"),
]

# Panel (b): two electrides, each under its conventional ionic counterpart.
SPREAD_CASES = [
    ("cs3auo/v3", "Cs3AuO"),
    ("ca3aun/v3", "Ca3AuN"),
    ("kmgf3/v3", "KMgF3"),
    ("k_hp/v3", "K"),
]

# Spot values quoted in docs/RESULTS.md, section 4 (lines 257-295); used only
# to cross-check the numbers computed here, never plotted.
RESULTS_MD_MEANS = {
    ("ba3si4", "Si"): -1.5, ("cs2pt", "Pt"): -2, ("k_hp", "K"): 0,
    ("ybafe2o5_a", "Fe"): 2.5, ("ybafe2o5_b", "Fe"): 2.5,
    ("mgh4", "H"): -0.5, ("sch4", "H"): -0.75, ("zrh4", "H"): -1,
    ("xef3", "Xe"): 3, ("xeo_b", "Xe"): 2, ("xeo2_a", "Xe"): 4,
    ("xeo3_a", "Xe"): 5, ("xeo3_b", "Xe"): 6,
    ("hgf4_a", "Hg"): 4, ("hgf4_b", "Hg"): 3, ("hgf4_c", "Hg"): 4,
    ("csf3", "F"): -1 / 3, ("csf5", "F"): -0.2, ("csh3", "H"): -1 / 3,
    ("lin5", "N"): -0.2, ("mgxe", "Xe"): -2,
}
RESULTS_MD_STRAY = {"ba3si4": 36, "xeo3_b": 9, "ceh4": 8, "thh4": 8,
                    "k_hp": 6, "mgxe": 5, "ybafe2o5_b": 4, "csh3": 3,
                    "mgh4": 3, "sch4": 3, "lin5": 1, "xeo3_a": 1}

TOL = 0.005        # two element means count as equal below this difference


# ----------------------------------------------------------------------------
# LEOS values: pyleos on the case directories, cached as JSON
# ----------------------------------------------------------------------------
def _digest(case_dir):
    h = hashlib.sha1()
    for fn in ("POSCAR", "EOSCAR", "wannier90.wout", "wannier90.1.wout",
               "wannier90.2.wout"):
        p = case_dir / fn
        if p.is_file():
            h.update(fn.encode())
            h.update(p.read_bytes())
    return h.hexdigest()


def _run_case(rel):
    """Default pyleos assignment for one case directory."""
    case = load_case(DATA / rel)
    res = assign_case(case)          # beta = 10, Allen, tier 1 (defaults)
    return {
        "elements": list(case.structure.elements),
        "os": [int(v) for v in res.os_int],
        "n_wf": len(res.wfs),
        "n_stray": len(res.stray),
        "ghost": int(res.ghost),
        "spread": [round(w.spread, 4) for w in res.wfs],
        "dist": [round(w.d1, 4) for w in res.wfs],
    }


def leos_values(recompute=False):
    cache = {}
    if CACHE.is_file() and not recompute:
        cache = json.loads(CACHE.read_text())
    out, dirty = {}, False
    for rel in [row[1] for row in LIT]:
        dig = _digest(DATA / rel)
        hit = cache.get(rel)
        if hit is None or hit.get("sha1") != dig:
            print(f"  pyleos: {rel}", flush=True)
            hit = _run_case(rel)
            hit["sha1"] = dig
            dirty = True
        out[rel] = hit
    if dirty:
        CACHE.parent.mkdir(exist_ok=True)
        CACHE.write_text(json.dumps(out, indent=0, sort_keys=True))
    return out


def element_stats(rec):
    by = defaultdict(list)
    for el, v in zip(rec["elements"], rec["os"]):
        by[el].append(v)
    return {el: (Fraction(sum(v), len(v)), min(v), max(v), len(v))
            for el, v in by.items()}


def leos_species(rec, sp):
    """(mean, min, max) of one species; 'OH' is the sum over O and H."""
    st = element_stats(rec)
    if sp == "OH":
        assert st["O"][3] == st["H"][3], "OH fragment needs as many O as H"
        return (st["O"][0] + st["H"][0], st["O"][1] + st["H"][1],
                st["O"][2] + st["H"][2])
    return st[sp][:3]


def parse_lit(s):
    """'-3/2' -> Fraction; 'a..b' -> (Fraction, Fraction) open interval."""
    if ".." in s:
        a, b = s.split("..")
        return (Fraction(a), Fraction(b))
    return Fraction(s)


# ----------------------------------------------------------------------------
# Comparison table
# ----------------------------------------------------------------------------
def build_points(leos):
    """One record per species per solid, plus one status per solid."""
    points, solids = [], []
    for (sid, rel, form, tag, cls, struct, lit, nline, page, basis) in LIT:
        rec = leos[rel]
        comparable = struct in ("same", "lattice")
        states = []
        for sp, raw in lit.items():
            x = parse_lit(raw)
            y, lo, hi = leos_species(rec, sp)
            if isinstance(x, tuple):
                # the literature mean lies strictly inside (x0, x1)
                state = "differ" if (y <= x[0] or y >= x[1]) else "open"
            else:
                state = "agree" if abs(float(y - x)) < TOL else "differ"
            states.append(state)
            points.append(dict(id=sid, rel=rel, formula=form, tag=tag,
                               species=sp, x=x, y=y, ymin=lo, ymax=hi,
                               state=state, comparable=comparable,
                               struct=struct, nline=nline, page=page,
                               basis=basis))
        if not lit or not comparable:
            status = "nc"
        elif "differ" in states:
            status = "differ"
        else:
            status = "agree"
        solids.append(dict(id=sid, rel=rel, formula=form, tag=tag, cls=cls,
                           struct=struct, status=status, nline=nline,
                           page=page, basis=basis, n_stray=rec["n_stray"],
                           ghost=rec["ghost"], n_wf=rec["n_wf"],
                           has_lit=bool(lit)))
    return points, solids


def fmt_frac(v):
    if isinstance(v, tuple):
        return f"({fmt_frac(v[0])},{fmt_frac(v[1])})"
    if v.denominator == 1:
        return f"{int(v):+d}"
    return f"{float(v):+.3f}"


def print_provenance(points, solids, leos):
    print("\nProvenance - literature value (Comas-Vila et al. 2025) versus "
          "LEOS element mean")
    print("  struct: same / lattice (approximate lattice) / different / "
          "unmatched;   notes = line of docs/comas_vila_2025_notes.md")
    head = (f"{'id':>2} {'solid':<18} {'LEOS dir (cases/)':<22} "
            f"{'struct':<9} {'sp':<3} {'lit':>10} {'LEOS':>7} "
            f"{'site range':>10} {'state':<7} {'stray':>5} {'ghost':>5} "
            f"{'notes':>5}  basis")
    print(head)
    print("-" * len(head))
    by_id = defaultdict(list)
    for p in points:
        by_id[p["id"]].append(p)
    for s in solids:
        name = s["formula"] + (f" [{s['tag']}]" if s["tag"] else "")
        rows = by_id.get(s["id"], [])
        if not rows:
            st = element_stats(leos[s["rel"]])
            own = ", ".join(f"{el} {fmt_frac(v[0])}" for el, v in st.items())
            print(f"{s['id']:>2} {name:<18} {s['rel']:<22} {s['struct']:<9} "
                  f"{'-':<3} {'none':>10} {'':>7} {'':>10} {'omitted':<7} "
                  f"{s['n_stray']:>5} {s['ghost']:>5} {'N' + str(s['nline']):>5}"
                  f"  {s['basis']}  [LEOS: {own}]")
        for k, p in enumerate(rows):
            lead = (f"{s['id']:>2} {name:<18} {s['rel']:<22} {s['struct']:<9}"
                    if k == 0 else f"{'':>2} {'':<18} {'':<22} {'':<9}")
            tail = (f"{s['n_stray']:>5} {s['ghost']:>5} "
                    f"{'N' + str(s['nline']):>5}  {s['basis']}" if k == 0 else "")
            rng = (f"{p['ymin']:+d}..{p['ymax']:+d}"
                   if p["ymin"] != p["ymax"] else "")
            print(f"{lead} {p['species']:<3} {fmt_frac(p['x']):>10} "
                  f"{fmt_frac(p['y']):>7} {rng:>10} {p['state']:<7} {tail}")
    notes = [(s["rel"], STRUCT_NOTE[s["rel"]]) for s in solids
             if s["rel"] in STRUCT_NOTE]
    print("\nStructure caveats")
    for rel, txt in notes:
        print(f"  {rel:<20} {txt}")

    n_pt = len(points)
    same = [p for p in points if p["comparable"]]
    other = [p for p in points if not p["comparable"]]
    print("\nPanel (a) points (species per solid)")
    print(f"  plotted                        {n_pt}  "
          f"({len({p['id'] for p in points})} solids)")
    print(f"  same structure, agree          "
          f"{sum(p['state'] == 'agree' for p in same)}")
    print(f"  same structure, differ         "
          f"{sum(p['state'] == 'differ' for p in same)}  "
          f"(of which intervals: "
          f"{sum(isinstance(p['x'], tuple) for p in same)})")
    print(f"  other structure, on diagonal   "
          f"{sum(p['state'] == 'agree' for p in other)}")
    print(f"  other structure, off diagonal  "
          f"{sum(p['state'] == 'differ' for p in other)}")
    print(f"  omitted (no literature value)  "
          f"{', '.join(s['formula'] for s in solids if not s['has_lit'])}")
    print("\nPanel (c) solids per class  (agree / differ / not comparable)")
    tot = defaultdict(int)
    for key, label in CLASSES:
        cnt = [sum(s["cls"] == key and s["status"] == st for s in solids)
               for st in ("agree", "differ", "nc")]
        for st, c in zip(("agree", "differ", "nc"), cnt):
            tot[st] += c
        members = ", ".join(
            s["formula"] + (f" [{s['tag']}]" if s["tag"] else "")
            + {"agree": "", "differ": " (differ)", "nc": " (n.c.)"}[s["status"]]
            for s in solids if s["cls"] == key)
        print(f"  {label.replace('$', '').replace('_', ''):<24} "
              f"{cnt[0]:>2} / {cnt[1]} / {cnt[2]}   {members}")
    print(f"  {'total':<24} {tot['agree']:>2} / {tot['differ']} / {tot['nc']}")

    stray = {s["rel"].split("/")[0]: s["n_stray"] for s in solids
             if s["n_stray"]}
    print("\nStray-WF warnings (spread > 3 A^2 or centre > 1.6 A from every "
          "atom): " + ", ".join(f"{k} {v}" for k, v in
                                sorted(stray.items(), key=lambda kv: -kv[1])))


def structure_facts():
    """Geometry behind the 'different' / 'unmatched' flags that is not already
    documented in docs/: read straight from the POSCAR files."""
    shifts = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1)
                       for k in (-1, 0, 1)], float)

    def dists(st, a, b):
        """Sorted a-b distances (A) over all pairs and periodic images."""
        lat = np.array(st.lattice)
        fr = np.array(st.frac)
        ia = [i for i, e in enumerate(st.elements) if e == a]
        ib = [i for i, e in enumerate(st.elements) if e == b]
        out = []
        for i in ia:
            for j in ib:
                d = np.linalg.norm((fr[j] - fr[i] + shifts) @ lat, axis=1)
                out += [x for x in d if x > 1e-6]
        return np.sort(out)

    print("\nStructure facts read from the POSCAR files")
    st = parse_poscar(DATA / "k_hp/v3/POSCAR")
    print(f"  k_hp     {len(st.elements)} K atoms in the cell "
          f"(the paper's cell: 4 K + 2 non-nuclear attractors, p. 7079)")
    for tag in "abc":
        st = parse_poscar(DATA / f"hgf4_{tag}/v3/POSCAR")
        n_hg = st.elements.count("Hg")
        hf = dists(st, "Hg", "F")
        near = hf[hf < 2.6]
        ff = dists(st, "F", "F")
        print(f"  hgf4_{tag}   {len(st.elements)} atoms; Hg-F below 2.6 A: "
              f"{len(near) // n_hg} per Hg, {near.min():.2f}-{near.max():.2f} A;"
              f" shortest F-F {ff.min():.2f} A")


def spread_summary(leos):
    """Numbers behind panel (b)."""
    print("\nPanel (b) - WF spreads (A^2) and centre-to-nearest-nucleus "
          "distances (A)")
    for rel, form in SPREAD_CASES:
        rec = leos[rel]
        sp, dist = np.array(rec["spread"]), np.array(rec["dist"])
        big = sp > STRAY_SPREAD
        line = (f"  {form:<7} {rel:<17} n = {len(sp):>2}; spread "
                f"{sp.min():.2f}-{sp.max():.2f}; above {STRAY_SPREAD:g} A^2: "
                f"{int(big.sum())}")
        if big.any():
            pairs = ", ".join(f"{a:.2f}/{b:.2f}" for a, b in
                              sorted(zip(sp[big], dist[big])))
            line += f"  (spread/distance: {pairs})"
        print(line)
        print(f"  {'':<7} {'':<17} other WFs: spread <= {sp[~big].max():.2f},"
              f" distance <= {dist[~big].max():.2f}")


def cross_check(leos):
    """Compare with docs/RESULTS.md sec. 4 and results/bench74_golden.tsv."""
    print("\nCross-checks")
    bad = 0
    for (case, el), ref in RESULTS_MD_MEANS.items():
        mean = float(element_stats(leos[f"{case}/v3"])[el][0])
        if abs(mean - ref) > TOL:
            bad += 1
            print(f"  MISMATCH {case} {el}: computed {mean:+.3f}, "
                  f"RESULTS.md {ref:+.3f}")
    print(f"  docs/RESULTS.md sec. 4 element means: "
          f"{len(RESULTS_MD_MEANS) - bad}/{len(RESULTS_MD_MEANS)} reproduced")
    got = {r.split("/")[0]: v["n_stray"] for r, v in leos.items()
           if INDEX[r]["role"] == "comas" and v["n_stray"]}
    print("  docs/RESULTS.md sec. 4 stray counts (12 of 30 cases): "
          + ("reproduced" if got == RESULTS_MD_STRAY else f"MISMATCH {got}"))
    tsv = {}
    for line in GOLDEN_TSV.read_text().splitlines()[1:]:
        q = line.split("\t")
        tsv[(q[0], q[1])] = q
    for rel, rec in leos.items():
        row = INDEX[rel]
        if row["role"] != "bench74":
            continue
        case = row["ref_key"]
        q = tsv[(row["bench_set"], case)]
        ok = int(q[6]) == rec["n_stray"] and int(q[7]) == rec["ghost"]
        print(f"  bench74_golden.tsv {case:<10} status {q[2]:<4} dev {q[3]}  "
              f"n_stray/n_ghost {q[6]}/{q[7]} "
              f"{'= computed' if ok else '!= computed'}")


# ----------------------------------------------------------------------------
# Drawing
# ----------------------------------------------------------------------------
FLOW = C["flow"]              # this work: colour ...
MK = MARKER["flow"]           # ... and marker, as in the other figures
D1 = 2.5                      # marker diameter (pt) of a single point
EDGE = 0.7                    # edge width of open markers (pt)


def signed(v, _pos=None):
    v = int(round(v))
    return "0" if v == 0 else (f"+{v}" if v > 0 else f"−{-v}")


def _sector(frac_filled):
    """Vertices of a disc sector covering frac_filled of the circle; the gap
    is centred on 315 deg (perpendicular to the diagonal)."""
    half_gap = 180.0 * (1.0 - frac_filled)
    ang = np.radians(np.linspace(315 + half_gap, 675 - half_gap, 60))
    verts = np.vstack([[0.0, 0.0], np.column_stack([np.cos(ang), np.sin(ang)])])
    return verts, np.abs(verts).max()


def draw_node(ax, x, y, n_same, n_other):
    """n_same filled + n_other open coincident points; area ~ n.  Smaller
    nodes are drawn above larger ones."""
    n = n_same + n_other
    d = D1 * np.sqrt(n)
    z = 3.0 + 1.0 / (n + 1)
    if n_other == 0:
        ax.scatter([x], [y], s=d ** 2, marker=MK, facecolor=FLOW,
                   edgecolor="white", linewidths=0.3, zorder=z)
        return
    ax.scatter([x], [y], s=(d - EDGE) ** 2, marker=MK, facecolor="white",
               edgecolor=FLOW, linewidths=EDGE, zorder=z)
    if n_same:
        verts, ext = _sector(n_same / n)
        ax.scatter([x], [y], s=(d * ext) ** 2, marker=verts, facecolor=FLOW,
                   edgecolor="none", zorder=z + 0.001)


def pt_label(p, star=False):
    name = formula(p["formula"]) + ("*" if star else "")
    if p["tag"] and not star:
        name += f", {p['tag']}"
    return name if p["formula"] == p["species"] else f"{name} ({p['species']})"


def panel_a(ax, points):
    lim = (-3.7, 6.7)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal", adjustable="box")
    ticks = list(range(-3, 7))
    for t in ticks:                                  # hairline integer grid
        ax.axhline(t, color=C["light"], lw=0.3, zorder=0)
        ax.axvline(t, color=C["light"], lw=0.3, zorder=0)
    ax.plot(lim, lim, color=C["muted"], lw=0.6, zorder=1)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.xaxis.set_major_formatter(FuncFormatter(signed))
    ax.yaxis.set_major_formatter(FuncFormatter(signed))
    ax.set_xlabel("Element-mean oxidation state (Comas-Vilà et al.)")
    ax.set_ylabel("Element-mean oxidation state (this work)")

    nodes = defaultdict(lambda: [0, 0])
    bars = []
    for p in points:
        if isinstance(p["x"], tuple):
            bars.append(p)
            continue
        key = (float(p["x"]), float(p["y"]))
        nodes[key][0 if p["comparable"] else 1] += 1
    for (x, y), (ns, no) in sorted(nodes.items(),
                                   key=lambda kv: -sum(kv[1])):
        draw_node(ax, x, y, ns, no)
    ebar = None
    for p in bars:
        x0, x1 = float(p["x"][0]), float(p["x"][1])
        ebar = ax.errorbar([(x0 + x1) / 2], [float(p["y"])],
                           xerr=[(x1 - x0) / 2], fmt="none", ecolor=FLOW,
                           elinewidth=1.0, capsize=2.2, capthick=0.8,
                           zorder=5)

    # Labels of every off-diagonal point: (id, species) -> offset in points.
    place = {
        (7, "Au"): dict(xytext=(-11, 6), ha="left", va="bottom"),
        (7, "Ca"): dict(xytext=(6, 0), ha="left", va="center"),
        (8, "K"): dict(xytext=(6, 0), ha="left", va="center"),
        (28, "Xe"): dict(xytext=(-6, 0), ha="right", va="center"),
        (28, "O"): dict(xytext=(6, -0.5), ha="left", va="center"),
        (29, "Xe"): dict(xytext=(-6, 0), ha="right", va="center"),
        (29, "O"): dict(xytext=(4, -8), ha="left", va="top"),
        (35, "Hg"): dict(xytext=(6, 0), ha="left", va="center"),
        (35, "F"): dict(xytext=(-3, 21), ha="right", va="bottom",
                        arrowprops=dict(arrowstyle="-", lw=0.4,
                                        color=C["muted"], shrinkA=0.5,
                                        shrinkB=2.5)),
    }
    off = [p for p in points if p["state"] == "differ"]
    missing = [(p["id"], p["species"]) for p in off
               if (p["id"], p["species"]) not in place]
    assert not missing, f"off-diagonal points without a label slot: {missing}"
    for p in off:
        kw = place[(p["id"], p["species"])]
        # an interval is labelled at its lower end
        x = float(p["x"][0] if isinstance(p["x"], tuple) else p["x"])
        # HgF4*: the one HgF4 cell off the diagonal (cell b, with F2 units)
        ax.annotate(pt_label(p, star=(p["id"] == 35)), xy=(x, float(p["y"])),
                    textcoords="offset points", fontsize=7, zorder=6, **kw)

    # Legends: kind of comparison, then marker-area key.
    h_same = Line2D([], [], ls="none", marker=MK, ms=D1 * np.sqrt(3),
                    mfc=FLOW, mec="white", mew=0.3)
    h_other = Line2D([], [], ls="none", marker=MK,
                     ms=D1 * np.sqrt(3) - EDGE, mfc="white", mec=FLOW, mew=EDGE)
    leg1 = ax.legend(
        [h_same, h_other, ebar],
        ["same crystal structure",
         "different or unmatched\ncrystal structure",
         "interval implied by\nsite-resolved values"],
        loc="lower right", bbox_to_anchor=(1.0, 0.115), handletextpad=0.6,
        labelspacing=0.55, borderaxespad=0.5)
    ax.add_artist(leg1)
    sizes = (1, 5, 15)
    hs = [Line2D([], [], ls="none", marker=MK, ms=D1 * np.sqrt(n), mfc=FLOW,
                 mec="white", mew=0.3) for n in sizes]
    ax.legend(hs, [str(n) for n in sizes], title="Coincident points",
              title_fontsize=7, loc="lower right", ncol=3, columnspacing=1.3,
              handletextpad=0.35, borderaxespad=0.5,
              bbox_to_anchor=(1.0, 0.0))


def panel_b(axes, leos):
    xlim = (0.3, 14.0)
    for ax, (rel, form) in zip(axes, SPREAD_CASES):
        rec = leos[rel]
        sp, dist = np.array(rec["spread"]), np.array(rec["dist"])
        ax.axvline(STRAY_SPREAD, color=C["muted"], lw=0.4, ls=(0, (4, 3)),
                   zorder=1)
        ax.scatter(sp, dist, s=2.8 ** 2, marker=MK, facecolor=FLOW,
                   edgecolor="none", alpha=0.55, zorder=3)
        ax.set_xscale("log")
        ax.set_xlim(xlim)
        ax.set_ylim(-0.3, 2.5)
        ax.set_yticks([0, 1, 2])
        ax.xaxis.set_major_locator(FixedLocator([0.5, 1, 2, 5, 10]))
        ax.xaxis.set_minor_locator(FixedLocator(
            [0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 3, 4, 6, 7, 8, 9]))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{v:g}"))
        name = formula(form) + (" ($oP$8)" if form == "K" else "")
        ax.text(0.022, 0.80, f"{name}, $n$ = {len(sp)}",
                transform=ax.transAxes, ha="left", va="top", fontsize=7)
        if ax is not axes[-1]:
            ax.tick_params(labelbottom=False)
    axes[-1].set_xlabel("Wannier-function spread (Å²)")


def panel_c(ax, solids):
    fills = {
        "agree": dict(facecolor=FLOW, edgecolor=FLOW),
        "differ": dict(facecolor="white", edgecolor=FLOW, hatch="//////"),
        "nc": dict(facecolor="white", edgecolor=FLOW),
    }
    names = {"agree": "agree", "differ": "differ", "nc": "not comparable"}
    for row, (key, _label) in enumerate(CLASSES):
        left = 0
        for st in ("agree", "differ", "nc"):
            n = sum(s["cls"] == key and s["status"] == st for s in solids)
            if n:
                ax.barh(row, n, left=left, height=0.62, lw=0.5, zorder=2,
                        **fills[st])
                left += n
    ax.set_ylim(len(CLASSES) - 0.45, -0.55)
    ax.set_yticks(range(len(CLASSES)))
    ax.set_yticklabels([lab for _k, lab in CLASSES],
                       va="center_baseline")
    ax.tick_params(axis="y", which="both", length=0, pad=4)
    ax.set_xlim(0, 9)
    ax.set_xticks([0, 2, 4, 6, 8])
    ax.xaxis.set_minor_locator(FixedLocator([1, 3, 5, 7, 9]))
    ax.set_xlabel("Number of solids")
    ax.legend([Patch(lw=0.5, **fills[st]) for st in ("agree", "differ", "nc")],
              [names[st] for st in ("agree", "differ", "nc")],
              loc="center right", bbox_to_anchor=(1.0, 0.70),
              handlelength=1.5, handleheight=0.75, labelspacing=0.45,
              borderaxespad=0.6)


def make_figure(points, solids, leos):
    pad = 0.02                               # inches, as used by style.save
    w_in = style.DOUBLE - 2 * pad
    h_cm = 9.05
    h_in = h_cm * CM - 2 * pad
    fig = plt.figure(figsize=(w_in, h_in))
    # Invisible frame: keeps the tight bounding box equal to the canvas, so
    # the saved figure is exactly one double column (17.1 cm) wide.
    fig.add_artist(Rectangle((0, 0), 1, 1, transform=fig.transFigure,
                             fill=False, edgecolor="none", lw=0))

    def axes_cm(left, bottom, width, height):
        return fig.add_axes([left * CM / w_in, bottom * CM / h_in,
                             width * CM / w_in, height * CM / h_in])

    bottom, top = 0.92, 8.52                 # common baseline / top (cm)
    side = top - bottom
    ax_a = axes_cm(1.08, bottom, side, side)
    panel_a(ax_a, points)
    panel_label(ax_a, "a", dx=-0.15 / side, dy=1 + 0.10 / side)

    right = 17.1 - 2 * pad / CM - 0.12       # right spine of (b) and (c)
    col = 9.95                               # right end of the (b)/(c) labels
    # (b): four rows sharing the spread axis
    b_left, row_h, gap = 10.75, 0.70, 0.10
    axes_b = []
    for i in range(4):
        y0 = top - (i + 1) * row_h - i * gap
        axes_b.append(axes_cm(b_left, y0, right - b_left, row_h))
    panel_b(axes_b, leos)
    b_bottom = top - 4 * row_h - 3 * gap
    fig.text((b_left - 0.78) * CM / w_in,
             (b_bottom + top) / 2 * CM / h_in,
             "Distance to nearest\nnucleus (Å)",
             rotation=90, ha="center", va="center", fontsize=8,
             linespacing=1.15)
    panel_label(axes_b[0], "b", dx=(col - b_left) / (right - b_left),
                dy=1 + 0.10 / row_h)

    # (c): solids per class
    c_left, c_top = 12.45, 4.15
    ax_c = axes_cm(c_left, bottom, right - c_left, c_top - bottom)
    panel_c(ax_c, solids)
    panel_label(ax_c, "c", dx=(col - c_left) / (right - c_left),
                dy=1 + 0.10 / (c_top - bottom))
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--recompute", action="store_true",
                    help="ignore the cache and re-run pyleos on all cases")
    args = ap.parse_args()

    assert len(LIT) == 40 and len({r[1] for r in LIT}) == 40
    leos = leos_values(args.recompute)
    points, solids = build_points(leos)
    print_provenance(points, solids, leos)
    spread_summary(leos)
    structure_facts()
    cross_check(leos)

    fig = make_figure(points, solids, leos)
    w, h = fig.get_size_inches()
    style.save(fig, NAME)
    print(f"\nwrote {NAME}.pdf / {NAME}.png  "
          f"(canvas {(w + 0.04) / CM:.2f} x {(h + 0.04) / CM:.2f} cm)")


if __name__ == "__main__":
    main()
