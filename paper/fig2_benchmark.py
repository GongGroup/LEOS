#!/usr/bin/env python3
"""Figure 2 - benchmark results of the min-cost-flow assignment.

(a) Verdicts of the two assignment methods on the scored cases of the
    74-case benchmark, per batch and in total.
(b) Largest element-mean deviation (dev) of both methods for every scored
    case that at least one of them does not pass.
(c) Per-atom accuracy on the clean development subset: three methods (left)
    and, for the min-cost flow, the cases that are not fully correct (right).
(d) Reliability index R_elem of the scored cases, grouped by outcome, with
    every case that does not pass named.
(e) Number of passing cases against the electronegativity weight beta.
(f) Element-mean oxidation states of this work against those of
    Comas-Vila et al. (J. Chem. Theory Comput. 2025) for their 40 solids.

Run from this directory:
    python3 fig2_benchmark.py               # uses the caches under _cache/
    python3 fig2_benchmark.py --recompute   # rebuild every cache first
    python3 fig2_benchmark.py --check-beta  # also repeat the beta scan on the
                                            # current data (printed, not plotted)

Needs numpy + matplotlib + the standard library, and the repository itself
(pyleos, including `pyleos bench clean`).  Nothing plotted is typed in by
hand, with one marked exception: the literature values of panel (f), table
LIT.  Every other number is read from results/ or computed by project code on
cases/.  The script prints a provenance table and the cross-checks against
docs/.
"""
import argparse
import csv
import hashlib
import inspect
import itertools
import json
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime
from fractions import Fraction

sys.dont_write_bytecode = True      # leave no __pycache__ in the repository

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import FuncFormatter, MultipleLocator

import style
from style import C, MARKER, LABEL, REPO, HERE, formula

sys.path.insert(0, str(REPO))
from pyleos.cases import load_index                         # noqa: E402
from pyleos.parsers import load_case, parse_poscar           # noqa: E402
from pyleos.transport import assign_case                     # noqa: E402

# ------------------------------------------------------------------ inputs
DATA = REPO / "cases"                             # <case>/<run>/ directories
RESULTS = REPO / "results"
TSV = {"flow": RESULTS / "bench74_golden.tsv",    # min-cost flow
       "nearest": RESULTS / "bench74_greedy.tsv"}  # per-WF greedy baseline
BETA_TSV = RESULTS / "scan_beta.tsv"
BETA_SCRIPT = REPO / "tools/beta_scan.py"
CACHE = HERE / "_cache"
BETA_CHECK = CACHE / "beta_scan_check.tsv"        # printed comparison only
LEOS_CACHE = CACHE / "fig3_leos.json"             # shared with figS1
EVALUATOR = REPO / "pyleos/bench.py"              # run as `pyleos bench clean`
BENCH_ARGS = {"flow": [],
              "geo": ["--method", str(REPO / "pyleos/methods/baseline_geo.py")],
              "nearest": ["--greedy"]}
FLOW_PLUGIN = REPO / "pyleos/methods/baseline.py"
RESULTS_MD = REPO / "docs/RESULTS.md"
# role / bench_set / ref_key of every run directory, keyed "<case>/<run>"
INDEX = {f"{r['case']}/{r['run']}": r for r in load_index(REPO)}

# Panel (b): also list the cases that neither method passes (set False to
# keep only the cases that exactly one method passes).
SHOW_NEITHER = True

# ------------------------------------------------------------- definitions
PASS_MAX, MID_MAX = 0.25, 1.0          # verdict limits of docs/RESULTS.md
VERDICTS = ("Pass", "Intermediate", "Fail")
STATUS = {"OK": 0, "~": 1, "FAIL": 2}  # status column -> verdict index
GROUPS = (("Legacy set", ("v1_mol", "v1_oxide")),
          ("Extension 1", ("wave1",)),
          ("Extension 2", ("wave2",)),
          ("All", ("v1_mol", "v1_oxide", "wave1", "wave2")))
TAG = {"v1_mol": "L1", "v1_oxide": "L2", "wave1": "E1", "wave2": "E2"}
# A Wannier centre counts as shared by two atoms when the two atoms nearest
# to it are within this distance ratio of each other (nearer / farther).
SHARED_RATIO = 0.5
# Method names in running text, taken from the shared style module.
SHORT = {m: LABEL[m].replace(" (this work)", "") for m in ("flow", "nearest")}

# Display formulas for the case directory names (labels only, no data);
# L = PPh3 (triphenylphosphine).
NAME = {
    "as4s4": "As4S4", "reo3": "ReO3", "cse2": "CSe2", "li2o2": "Li2O2",
    "ref6": "ReF6", "cro5": "CrO5", "ca3aun": "Ca3AuN", "h3po3": "H3PO3",
    "hi3": "HI3", "kmgf3": "KMgF3", "mocl2": "MoCl2", "n4s4": "S4N4",
    "rb2cucl4": "Rb2CuCl4", "rbpb": "RbPb", "rhh": "RhH",
    "rhcl_co__p_c6h5_3_2": "RhCl(CO)L2", "srbe_oh_4": "SrBe(OH)4",
    "wau12": "WAu12", "wcl3": "WCl3", "wcl4": "WCl4", "wcl5": "WCl5",
    "wcl6": "WCl6", "re2": "Re2", "nan3": "NaN3", "h2s2": "H2S2",
    "h2s3": "H2S3", "h2s5": "H2S5", "k2s2o3": "K2S2O3", "sega": "GaSe",
    "co": "CO", "tio2": "TiO2", "tio": "TiO", "fe3o4": "Fe3O4",
    "tmao": "(CH3)3NO", "ch2f2": "CH2F2", "cr_co_6": "Cr(CO)6",
}
# Panel (d): cases whose reference follows another convention than the
# electronegativity rule (docs/RESULTS.md, sections 3 and 8); asterisk.
CONVENTION = {("wave1", "ca3aun"), ("wave2", "wau12")}
# Panel (d): named cases whose R_elem differ by less than this (about a
# quarter of the marker diameter) are bracketed and share one leader line.
R_TIE = 1.0

# Values quoted in docs/RESULTS.md (sections 2, 7 and 8) - used for the
# printed cross-check only, never plotted.
QUOTED = {"n_ok": 64, "n_not": 8, "med_ok": 83.0, "med_not": 25.55,
          "beta10_scan": 62}

# =============================================================================
# LITERATURE VALUES of panel (f) - the only hand-typed numbers in this script.
# Copied unchanged from figS1_comparison.py.
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
#   cls       chemical class (provenance table only)
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
    ("salt", "Oxides, halides, GaSe"), ("auride", "Aurides"),
    ("electride", "Electrides"), ("zintl", "Zintl phase"),
    ("platinide", "Platinide"), ("mixed", "Mixed valence"),
    ("mh4", "MH4 hydrides"), ("xeo", "High-pressure Xe oxides"),
    ("xef", "Xe fluorides"), ("ngan", "Noble-gas anion"),
    ("hg4", "Hg(IV) fluoride"), ("poly", "Polyatomic anions"),
]
# Spot values quoted in docs/RESULTS.md, section 4; used only to cross-check
# the numbers computed here, never plotted.
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


def rel(path):
    try:
        return str(path.relative_to(REPO))
    except ValueError:
        return str(path)


def lc(s):
    """Lower-case first letter, for a method name inside a phrase."""
    return s[0].lower() + s[1:]


def upright(s):
    """Plain text as mathtext, so that it is set like the formulas."""
    return rf"$\mathrm{{{s}}}$"


def pretty(case):
    s = formula(NAME.get(case.lower(), case))
    return s if "$" in s else upright(s)


# ------------------------------------------------------- data: (a, b, d, e)
def read_rows(path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def verdict_of(dev):
    return "OK" if dev <= PASS_MAX else "~" if dev <= MID_MAX else "FAIL"


def read_tsv(path):
    """{(batch, case): (status, dev)} for the scored cases of one table."""
    out = {}
    for row in read_rows(path):
        if row["status"] == "run":          # no reference answer
            continue
        dev = float(row["dev"])
        # the status column is the verdict; dev is printed to 0.01, so
        # only check that the two agree within that rounding
        if row["status"] not in {verdict_of(dev + e) for e in (-.005, 0, .005)}:
            raise ValueError(f"{path.name}: {row['case']} status "
                             f"{row['status']} disagrees with dev {dev}")
        out[(row["batch"], row["case"])] = (row["status"], dev)
    return out


def verdict_counts(table, batches):
    c = Counter(STATUS[s] for (b, _), (s, _) in table.items() if b in batches)
    return [c[0], c[1], c[2]]


def docs_numbers():
    """Pass counts quoted in docs/RESULTS.md (for the cross-check only)."""
    text = RESULTS_MD.read_text(encoding="utf-8")
    keys = {"Legacy regression set": "Legacy set",
            "Extension wave 1": "Extension 1",
            "Extension wave 2": "Extension 2", "Total": "All"}
    out = {}
    for row, en in keys.items():
        m = re.search(rf"\|\s*{row} \((\d+)\)\s*\|\s*(\d+)[^|]*\|\s*(\d+)",
                      text)
        if m:
            out[en] = tuple(int(x) for x in m.groups())
    m = re.search(r"(\d+)\s+cases\s+pass\s+only\s+with\s+the\s+flow,\s+(\d+)\s+"
                  r"pass\s+only\s+with\s+greedy", text)
    out["only"] = tuple(int(x) for x in m.groups()) if m else None
    m = re.search(r"Within\s+wave\s+2\s+these\s+two\s+numbers\s+are\s+(\d+)\s+and\s+"
                  r"(\d+)", text)
    out["only_wave2"] = tuple(int(x) for x in m.groups()) if m else None
    return out


def beta_check(run=False):
    """Beta scan of the data tree as it is now (tools/beta_scan.py), kept in
    _cache/.  It is compared with the stored scan in the printed provenance
    and is never plotted.  Returns its rows, or None if it was never run."""
    if run:
        CACHE.mkdir(exist_ok=True)
        cmd = [sys.executable, str(BETA_SCRIPT), "--out", str(BETA_CHECK)]
        print("running:", " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                       env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
    return read_rows(BETA_CHECK) if BETA_CHECK.is_file() else None


# --------------------------------------------------------------- data: (c)
def clean_bench(method, recompute=False):
    """Dev-split result of the frozen evaluator (cached JSON, else run it)."""
    path = CACHE / f"{method}.json"
    if recompute or not path.is_file():
        CACHE.mkdir(exist_ok=True)
        cmd = [sys.executable, "-m", "pyleos", "bench", "clean",
               "--split", "dev", *BENCH_ARGS[method], "--json", str(path)]
        print("running:", " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, cwd=REPO,
                       env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
    res = json.loads(path.read_text())
    if res["split"] != "dev":
        raise ValueError(f"{path} was not produced with --split dev")
    return res


_IMAGES = np.array(list(itertools.product(range(-2, 3), repeat=3)), float)


def _distances(lattice, frac, point):
    """Distances (n_atom, n_image) from a Cartesian point to all atoms."""
    df = frac - point @ np.linalg.inv(lattice)
    df -= np.round(df)
    return np.linalg.norm((df[:, None, :] + _IMAGES[None]) @ lattice, axis=2)


def like_atom_contacts(flow, recompute=False):
    """Structure check of the per-atom errors of the min-cost flow.

    For each development case that is not fully correct, re-run the frozen
    flow plugin to locate the wrong sites and report (i) the shortest
    distance from a wrong site to another atom of the same element and
    (ii) whether every wrong site is one of the two atoms sharing a Wannier
    centre with an atom of the same element.  Cached in _cache/.
    """
    path = CACHE / "flow_contacts.json"
    todo = {(r["batch"], r["case"]): r for r in flow["cases"]
            if r["atoms_ok"] != r["atoms_labelled"]}
    if path.is_file() and not recompute:
        cached = json.loads(path.read_text())
        if (isinstance(cached, dict)
                and cached.get("shared_ratio") == SHARED_RATIO
                and {(r["batch"], r["case"]) for r in cached["rows"]} == set(todo)):
            return cached["rows"]
    import pyleos.bench as cb                      # frozen evaluator
    plugin = cb.load_plugin(FLOW_PLUGIN)
    spec = {(b, n): s for b, cases in cb.CLEAN["cases"].items()
            for n, s in cases.items()}
    rows = []
    for batch, name, d, ref in cb.collect_clean():   # same dirs as the evaluator
        if (batch, name) not in todo:
            continue
        case = load_case(d)
        st = case.structure
        lat, frac, el = np.array(st.lattice), np.array(st.frac), st.elements
        centres = [np.array(c) for ch in case.channels for c in ch.centres]
        labels = cb.site_labels(name, st, ref, spec[(batch, name)])
        os_vals = plugin.assign(cb.scrub(case))
        wrong = [i for i, (lab, o) in enumerate(zip(labels, os_vals))
                 if lab is not None and o != lab]
        r = todo[(batch, name)]
        if len(wrong) != r["atoms_labelled"] - r["atoms_ok"]:
            raise RuntimeError(f"{name}: wrong sites do not match the evaluator")
        cart = frac @ lat
        shortest = []
        for i in wrong:
            dist = _distances(lat, frac, cart[i])
            same = [dist[j][dist[j] > 1e-6].min() for j in range(len(el))
                    if el[j] == el[i]]
            shortest.append(float(min(same)))
        on_shared, n_shared = set(), 0
        for c in centres:
            dist = _distances(lat, frac, c)
            k1, k2 = np.argsort(dist, axis=None)[:2]
            j1, j2 = k1 // dist.shape[1], k2 // dist.shape[1]
            if el[j1] == el[j2] and dist.flat[k1] >= SHARED_RATIO * dist.flat[k2]:
                n_shared += 1
                on_shared.update((int(j1), int(j2)))
        rows.append({
            "batch": batch, "case": name, "dir": str(d.relative_to(REPO)),
            "elements": sorted({el[i] for i in wrong}),
            "n_wrong": len(wrong), "n_labelled": r["atoms_labelled"],
            "d_min": round(min(shortest), 4), "d_max": round(max(shortest), 4),
            "n_wf": len(centres), "n_shared_wf": n_shared,
            "wrong_on_shared_wf": len(set(wrong) & on_shared),
            "evaluator_wrong": r["wrong"]})
    CACHE.mkdir(exist_ok=True)
    path.write_text(json.dumps({"shared_ratio": SHARED_RATIO, "rows": rows},
                               indent=1) + "\n")
    return rows


# --------------------------------------------------------------- data: (f)
def _digest(case_dir):
    h = hashlib.sha1()
    for fn in ("POSCAR", "EOSCAR", "wannier90.wout", "wannier90.1.wout",
               "wannier90.2.wout"):
        p = case_dir / fn
        if p.is_file():
            h.update(fn.encode())
            h.update(p.read_bytes())
    return h.hexdigest()


def _run_case(rel_dir):
    """Default pyleos assignment for one case directory (same record as
    figS1_comparison.py writes, so that the two scripts share the cache)."""
    case = load_case(DATA / rel_dir)
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
    """LEOS result of the 40 solids; cache keyed by the SHA-1 of the inputs."""
    cache = {}
    if LEOS_CACHE.is_file() and not recompute:
        cache = json.loads(LEOS_CACHE.read_text())
    out, dirty = {}, False
    for rel_dir in [row[1] for row in LIT]:
        dig = _digest(DATA / rel_dir)
        hit = cache.get(rel_dir)
        if hit is None or hit.get("sha1") != dig:
            print(f"  pyleos: {rel_dir}", flush=True)
            hit = _run_case(rel_dir)
            hit["sha1"] = dig
            dirty = True
        out[rel_dir] = hit
    if dirty:
        LEOS_CACHE.parent.mkdir(exist_ok=True)
        LEOS_CACHE.write_text(json.dumps(out, indent=0, sort_keys=True))
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


def build_points(leos):
    """One record per species per solid, plus one status per solid."""
    points, solids = [], []
    for (sid, rel_dir, form, tag, cls, struct, lit, nline, page, basis) in LIT:
        rec = leos[rel_dir]
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
            points.append(dict(id=sid, rel=rel_dir, formula=form, tag=tag,
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
        solids.append(dict(id=sid, rel=rel_dir, formula=form, tag=tag, cls=cls,
                           struct=struct, status=status, nline=nline,
                           page=page, basis=basis, n_stray=rec["n_stray"],
                           ghost=rec["ghost"], n_wf=rec["n_wf"],
                           has_lit=bool(lit)))
    return points, solids


# ----------------------------------------------------------------- drawing
PT = 2.54 / 72          # cm per point


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


def bar_height(ax, thickness_cm=0.20):
    """Bar height in data units that gives the same thickness everywhere."""
    y0, y1 = ax.get_ylim()
    return thickness_cm * abs(y1 - y0) / ax_cm(ax)[1]


def method_glyph(ax, y, method, offset_cm=0.17):
    """Method marker left of a bar: second cue besides the bar colour."""
    ax.plot([-offset_cm / ax_cm(ax)[0]], [y], ls="none", marker=MARKER[method],
            ms=3.6, color=C[method], mec=C[method], clip_on=False,
            transform=ax.get_yaxis_transform())


def categorical(ax):
    ax.tick_params(axis="y", which="both", left=False, right=False)


def baseline_shift(s, size):
    """Points to lower a baseline so that capitals are centred on the anchor.

    Vector backends draw mathtext about 0.9 pt above the requested baseline
    (the box height is rounded up to a whole point), the raster backend does
    not; strings containing mathtext are lowered by half of that, which keeps
    them within 0.45 pt of plain text in both the PDF and the PNG."""
    return -0.358 * size - (0.45 if "$" in s else 0.0)


def row_text(ax, x, y, s, ha="left", transform=None):
    """Text on a row, baseline-aligned (centre alignment would move every
    label by its own extent: subscripts, parentheses, the ring of the
    angstrom sign)."""
    size = plt.rcParams["font.size"]
    tr = mtransforms.offset_copy(transform or ax.transData, fig=ax.figure,
                                 y=baseline_shift(s, size), units="points")
    return ax.text(x, y, s, ha=ha, va="baseline", transform=tr)


def method_legend(fig, x_cm, y_cm, methods, gap_cm=0.65):
    """One row of marker + label per method; labels share one baseline."""
    fw, fh = fig.get_size_inches() / style.CM
    size = plt.rcParams["legend.fontsize"]
    renderer = fig.canvas.get_renderer()
    y = 1 - y_cm / fh
    for m in methods:
        fig.add_artist(Line2D([(x_cm + 0.09) / fw], [y], ls="none",
                              marker=MARKER[m], ms=4.2, color=C[m], mec=C[m],
                              transform=fig.transFigure))
        tr = mtransforms.offset_copy(fig.transFigure, fig=fig,
                                     y=baseline_shift(LABEL[m], size),
                                     units="points")
        t = fig.text((x_cm + 0.30) / fw, y, LABEL[m], ha="left", va="baseline",
                     fontsize=size, transform=tr)
        x_cm = t.get_window_extent(renderer).x1 / fig.dpi / style.CM + gap_cm


def align_ylabel(ax, ax_ref):
    """Shift the y label of ax so that its left edge (the top of its first
    line) lies on the left edge of the y label of ax_ref."""
    renderer = ax.figure.canvas.get_renderer()
    lab, ref = ax.yaxis.label, ax_ref.yaxis.label
    dx = (ref.get_window_extent(renderer).x0
          - lab.get_window_extent(renderer).x0)              # pixels
    x, y = lab.get_position()
    ax.yaxis.set_label_coords(x + dx / ax.bbox.width, y)


def top_legend(ax, handles, labels, **kw):
    """Legend in the band above an axes, starting at its left edge."""
    opts = dict(loc="lower left", ncol=len(handles), bbox_to_anchor=(0, 1),
                borderaxespad=0.35, borderpad=0, columnspacing=1.4)
    opts.update(kw)
    return ax.legend(handles, labels, **opts)


# ---------------------------------------------------------------- panel (a)
def panel_a(ax, tables):
    step, gap = 1.0, 0.9
    ticks, labels, prov = [], [], []
    col_cm = (0.50, 1.15, 1.80)             # column centres right of the axes
    w_cm = ax_cm(ax)[0]
    n_groups = len(GROUPS)
    ax.set_ylim(n_groups * (2 * step + gap) - gap - step + 0.75, -0.75)
    ax.set_xlim(0, 100)
    h = bar_height(ax)
    for g, (name, batches) in enumerate(GROUPS):
        y0 = g * (2 * step + gap)
        n_tot = None
        for k, method in enumerate(("flow", "nearest")):
            y = y0 + k * step
            counts = verdict_counts(tables[method], batches)
            n = sum(counts)
            n_tot = n if n_tot is None else n_tot
            assert n == n_tot
            left = 0.0
            for v, cnt in enumerate(counts):
                width = 100.0 * cnt / n
                ax.barh(y, width, left=left, height=h, lw=0.5,
                        color=C[method] if v == 0 else "white",
                        edgecolor=C[method], hatch="//////" if v == 1 else None)
                left += width
                row_text(ax, 1 + col_cm[v] / w_cm, y, str(cnt), ha="center",
                         transform=ax.get_yaxis_transform())
            method_glyph(ax, y, method)
            prov.append((name, method, n, counts))
        ticks.append(y0 + step / 2)
        labels.append(f"{name}\n($n$ = {n_tot})")
    ax.set_yticks(ticks, labels)
    ax.tick_params(axis="y", pad=9)
    categorical(ax)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Scored cases (%)")
    # fill styles: legend above the bars, the same swatches head the columns
    fills = [dict(facecolor=C["muted"] if v == 0 else "white",
                  edgecolor=C["muted"], lw=0.5,
                  hatch="//////" if v == 1 else None) for v in range(3)]
    top_legend(ax, [Patch(**f) for f in fills], VERDICTS, handleheight=0.75)
    h_cm = ax_cm(ax)[1]
    gw, gh = 0.34 / w_cm, 0.18 / h_cm
    for v, f in enumerate(fills):
        ax.add_patch(Rectangle((1 + col_cm[v] / w_cm - gw / 2,
                                1 + 0.21 / h_cm - gh / 2), gw, gh,
                               clip_on=False, transform=ax.transAxes, **f))
    return prov


# ---------------------------------------------------------------- panel (b)
def panel_b(ax, tables):
    flow, near = tables["flow"], tables["nearest"]
    blocks = [(f"Pass: {lc(SHORT['flow'])} only", [], lambda f, n: (-n, -f)),
              (f"Pass: {lc(SHORT['nearest'])} only", [], lambda f, n: (-f, -n)),
              ("Pass: neither", [], lambda f, n: (-f, -n))]
    n_both = 0
    for key in flow:
        (sf, df), (sn, dn) = flow[key], near[key]
        if sf == "OK" and sn == "OK":
            n_both += 1
            continue
        blocks[0 if sf == "OK" else 1 if sn == "OK" else 2][1].append(
            (key, df, dn))
    counts = [len(rows) for _, rows, _ in blocks]
    if not SHOW_NEITHER:
        blocks = blocks[:2]
    y, prov = 0.0, []
    label_x = -0.50 / ax_cm(ax)[0]
    tag_x = -0.42 / ax_cm(ax)[0]
    yt = ax.get_yaxis_transform()
    for title, rows, order in blocks:
        rows.sort(key=lambda r: (order(r[1], r[2]), r[0][1].lower(), r[0][0]))
        # the white box interrupts the threshold line behind the heading
        row_text(ax, 0.10 / ax_cm(ax)[0], y, f"{title} ($n$ = {len(rows)})",
                 transform=yt).set_bbox(dict(fc="white", ec="none", pad=1.0))
        for key, df, dn in rows:
            y += 1
            ax.plot([df, dn], [y, y], color=C["muted"], lw=0.8, zorder=1,
                    solid_capstyle="butt")
            ax.plot(dn, y, ls="none", marker=MARKER["nearest"], ms=4.4,
                    color=C["nearest"], mec=C["nearest"], zorder=2)
            ax.plot(df, y, ls="none", marker=MARKER["flow"], ms=3.6,
                    color=C["flow"], mec=C["flow"], zorder=3)
            row_text(ax, label_x, y, pretty(key[1]), ha="right", transform=yt)
            row_text(ax, tag_x, y, upright(TAG[key[0]]), transform=yt)
            prov.append((title, key, df, dn))
        y += 1.4
    ax.set_ylim(y - 1.4 + 0.65, -0.95)
    ax.set_xlim(-0.42, 8.5)
    ax.axvline(PASS_MAX, color=C["ink"], lw=0.5, ls=(0, (4, 2.5)),
               zorder=0)                  # pass threshold
    ax.text(PASS_MAX, 1 + 0.10 / ax_cm(ax)[1], f"dev = {PASS_MAX:g}",
            ha="left", va="bottom", transform=ax.get_xaxis_transform())
    ax.set_yticks([])
    ax.set_xticks(range(0, 9, 2))
    ax.set_xticks(range(0, 9, 1), minor=True)
    ax.set_xlabel("Deviation from reference, dev")
    return prov, n_both, counts


# ---------------------------------------------------------------- panel (c)
def panel_c_summary(ax, bench):
    methods = ("flow", "nearest", "geo")
    rows = (("Atoms\ncorrect", "atoms_ok", "atoms_labelled"),
            ("Cases fully\ncorrect", "all_atoms_correct_cases", "cases"))
    step, gap = 1.0, 1.0
    ax.set_ylim(len(rows) * (3 * step + gap) - gap - step + 0.8, -0.8)
    ax.set_xlim(0, 100)
    h = bar_height(ax)
    ticks = []
    for g, (_, num, den) in enumerate(rows):
        y0 = g * (3 * step + gap)
        for k, m in enumerate(methods):
            s = bench[m]["summary"]
            y = y0 + k * step
            ax.barh(y, 100.0 * s[num] / s[den], height=h, color=C[m],
                    edgecolor=C[m], lw=0.5)
            method_glyph(ax, y, m)
            row_text(ax, 1 + 0.12 / ax_cm(ax)[0], y, f"{s[num]}/{s[den]}",
                     transform=ax.get_yaxis_transform())
        ticks.append(y0 + step)
    ax.set_yticks(ticks, [r[0] for r in rows])
    ax.tick_params(axis="y", pad=9)
    categorical(ax)
    ax.set_xticks([0, 50, 100])
    ax.set_xticks([25, 75], minor=True)
    ax.set_xlabel("Fraction (%)")


def panel_c_cases(ax, contacts):
    rows = sorted(contacts, key=lambda r: (-r["n_wrong"], -r["n_labelled"]))
    ax.set_ylim(len(rows) - 1 + 0.75, -0.75)
    xmax = max(r["n_wrong"] for r in rows)
    ax.set_xlim(0, xmax * 1.08)
    h = bar_height(ax)
    w_cm = ax_cm(ax)[0]
    yt = ax.get_yaxis_transform()
    for y, r in enumerate(rows):
        ax.barh(y, r["n_wrong"], height=h, color=C["flow"],
                edgecolor=C["flow"], lw=0.5)
        method_glyph(ax, y, "flow")
        row_text(ax, -0.34 / w_cm, y, pretty(r["case"]), ha="right",
                 transform=yt)
        row_text(ax, 1 + 0.75 / w_cm, y, f"{r['n_wrong']}/{r['n_labelled']}",
                 ha="right", transform=yt)
        if len(r["elements"]) != 1:
            raise RuntimeError(f"{r['case']}: wrong sites of several elements")
        el = r["elements"][0]
        bonded = r["wrong_on_shared_wf"] == r["n_wrong"]
        link = "–" if bonded else "···"
        row_text(ax, 1 + 0.90 / w_cm, y,
                 f"{el}{link}{el} {r['d_min']:.2f} Å", transform=yt)
    ax.set_yticks([])
    ax.set_xticks(range(0, xmax + 1, 8))
    ax.set_xticks(range(0, xmax + 1, 4), minor=True)
    ax.set_xlabel("Incorrect atoms")
    return rows


# ---------------------------------------------------------------- panel (d)
def swarm(y_pt, dist, step=0.5):
    """Deterministic bee swarm: horizontal offsets (points) such that no two
    marker centres are closer than `dist` points."""
    xs = np.zeros(len(y_pt))
    placed = []
    for i in np.argsort(y_pt, kind="stable"):
        k = 0
        while True:
            hit = False
            for cand in ((0.0,) if k == 0 else (k * step, -k * step)):
                if all((cand - px) ** 2 + (y_pt[i] - py) ** 2 >= dist ** 2
                       for px, py in placed):
                    xs[i] = cand
                    hit = True
                    break
            if hit:
                break
            k += 1
        placed.append((xs[i], y_pt[i]))
    return xs


def spread_labels(targets, sep, lo, hi):
    """Label positions close to `targets`, at least `sep` apart, in [lo, hi]."""
    order = np.argsort(targets, kind="stable")
    blocks = []                                    # [centre, members]
    for i in order:
        blocks.append([float(targets[i]), [i]])
        while len(blocks) > 1:
            (c1, m1), (c2, m2) = blocks[-2], blocks[-1]
            if c1 + (len(m1) - 1) / 2 * sep + sep <= c2 - (len(m2) - 1) / 2 * sep:
                break
            n1, n2 = len(m1), len(m2)
            blocks[-2:] = [[(c1 * n1 + c2 * n2) / (n1 + n2), m1 + m2]]
    out = np.zeros(len(targets))
    for centre, members in blocks:
        half = (len(members) - 1) / 2 * sep
        centre = min(max(centre, lo + half), hi - half)
        for k, i in enumerate(members):
            out[i] = centre - half + k * sep
    return out


def panel_d(ax, rows):
    """R_elem of the scored cases by outcome (x axis in cm from the left)."""
    scored = [dict(r) for r in rows if r["status"] != "run"]
    groups = [[r for r in scored if r["status"] == "OK"],
              [r for r in scored if r["status"] != "OK"]]
    repeated = Counter(r["case"].lower() for r in scored)

    def name(r):
        key = (r["batch"], r["case"])
        tag = f" ({TAG[r['batch']]})" if repeated[r["case"].lower()] > 1 else ""
        return pretty(r["case"]) + tag + ("*" if key in CONVENTION else "")

    w_cm, h_cm = ax_cm(ax)
    ax.set_xlim(0, w_cm)
    ax.set_ylim(-4, 104)
    pt_per_y = h_cm / PT / 108.0
    centre = (1.60, 4.35)                   # cm: the two groups
    ms, med = 3.3, []
    for gi, grp in enumerate(groups):
        r_el = np.array([float(r["R_elem"]) for r in grp])
        stray = np.array([int(r["n_stray"]) > 0 for r in grp])
        q1, q2, q3 = np.percentile(r_el, [25, 50, 75])
        med.append(float(np.median(r_el)))
        ax.add_patch(Rectangle((centre[gi] - 0.60, q1), 1.20, q3 - q1,
                               fc=C["light"], ec="none", alpha=0.75, zorder=1))
        ax.plot([centre[gi] - 0.70, centre[gi] + 0.70], [q2, q2],
                color=C["ink"], lw=0.9, solid_capstyle="butt", zorder=5)
        x = centre[gi] + swarm(r_el * pt_per_y, dist=ms + 0.5) * PT
        ax.plot(x[~stray], r_el[~stray], ls="none", marker=MARKER["flow"],
                ms=ms, mfc=C["flow"], mec="white", mew=0.3, zorder=4)
        ax.plot(x[stray], r_el[stray], ls="none", marker=MARKER["aux3"],
                ms=ms + 1.0, mfc=C["flow"], mec="white", mew=0.3, zorder=4)
        for r, xi in zip(grp, x):
            r["_x"] = xi

    def leader(x0, y0, x1, y1):
        dx, dy = (x1 - x0) / PT, (y1 - y0) * pt_per_y            # points
        gap = 2.6 / np.hypot(dx, dy)
        ax.plot([x0 + gap * (x1 - x0), x1], [y0 + gap * (y1 - y0), y1],
                color=C["muted"], lw=0.4, zorder=3, solid_capstyle="round")

    def labelled(cases, x_lab):
        """Names right of the markers.  Cases whose R_elem lie within R_TIE
        of each other cannot be told apart at this scale: they share one
        leader, from the right-most marker to the middle of a bracket around
        their names (largest R_elem at the top).  Returns the tied groups."""
        cases = sorted(cases, key=lambda r: float(r["R_elem"]))
        y_pts = np.array([float(r["R_elem"]) for r in cases])
        y_lab = spread_labels(y_pts, 9.2 / pt_per_y, lo=1.5, hi=98.5)
        runs = [[0]]
        for i in range(1, len(cases)):
            if y_pts[i] - y_pts[i - 1] < R_TIE:
                runs[-1].append(i)
            else:
                runs.append([i])
        for run in runs:
            for i in run:
                ax.text(x_lab, y_lab[i], name(cases[i]), ha="left",
                        va="center")
            if len(run) == 1:
                i = run[0]
                leader(cases[i]["_x"], y_pts[i], x_lab - 0.08, y_lab[i])
                continue
            xb, tick, half = x_lab - 0.07, 0.05, 4.3 / pt_per_y
            lo_b, hi_b = y_lab[run[0]] - half, y_lab[run[-1]] + half
            ax.plot([xb + tick, xb, xb, xb + tick], [lo_b, lo_b, hi_b, hi_b],
                    color=C["muted"], lw=0.5, zorder=3,
                    solid_capstyle="butt", solid_joinstyle="miter")
            i = max(run, key=lambda i: cases[i]["_x"])   # right-most marker
            leader(cases[i]["_x"], y_pts[i], xb, (lo_b + hi_b) / 2)
        return [[cases[i] for i in run] for run in runs if len(run) > 1]

    # every case that does not pass, and the passing cases below the median
    # of that group
    low_ok = [r for r in groups[0] if float(r["R_elem"]) < med[1]]
    tied = labelled(groups[1], centre[1] + 0.88)
    labelled(low_ok, centre[0] + 0.60)

    ax.set_xticks(list(centre))
    ax.set_xticklabels([f"{t}\n($n$ = {len(g)})" for t, g in
                        zip(("Pass", "Intermediate or fail"), groups)])
    ax.tick_params(axis="x", which="both", length=0, pad=4)
    ax.yaxis.set_major_locator(MultipleLocator(20))
    ax.yaxis.set_minor_locator(MultipleLocator(10))
    ax.set_ylabel(r"Reliability index $R_{\mathrm{elem}}$")
    ax.yaxis.set_label_coords(-0.95 / w_cm, 0.5)
    handles = [
        Line2D([], [], ls="none", marker=MARKER["flow"], ms=ms, mfc=C["flow"],
               mec="white", mew=0.3),
        Line2D([], [], ls="none", marker=MARKER["aux3"], ms=ms + 1.0,
               mfc=C["flow"], mec="white", mew=0.3)]
    top_legend(ax, handles, ["no stray WF", "≥ 1 stray WF"],
               handletextpad=0.2)
    return dict(groups=groups, med=med, low_ok=low_ok, tied=tied,
                n_stray=[sum(int(r["n_stray"]) > 0 for r in g) for g in groups])


# ---------------------------------------------------------------- panel (e)
E_STEP = 8               # spacing of the labelled y ticks (cases)
E_LINESPACING = 1.15     # of the two-line y label


def panel_e(ax, rows, n_pass_now):
    beta = np.array([float(r["beta"]) for r in rows])
    n_ok = np.array([int(r["n_OK"]) for r in rows])
    n_all = {int(r["n_OK"]) + int(r["n_mid"]) + int(r["n_FAIL"]) for r in rows}
    assert len(n_all) == 1, "beta scan: number of scored cases varies"
    n_all = n_all.pop()
    # the reference (default) beta is encoded in the column name
    col = [k for k in rows[0] if k.startswith("changed_vs_beta")][0]
    beta0 = float(re.search(r"beta([\d.]+)$", col).group(1))
    at0 = int(n_ok[np.argmin(np.abs(beta - beta0))])

    ax.axvline(beta0, color=C["muted"], lw=0.5, ls=(0, (4, 3)), zorder=1)
    ax.plot(beta, n_ok, color=C["flow"], lw=1.0, marker=MARKER["flow"], ms=3.8,
            mfc=C["flow"], mec="white", mew=0.4, zorder=3)
    ax.set_xlim(1.5, 21.5)
    # The y axis ends at the number of scored cases (the top of the frame is
    # "all cases pass") and starts on the tick grid, at least two cases below
    # the lowest count: the changes along the scan are seen against the full
    # count and are not magnified by a tight range.
    y_lo = n_all - E_STEP * int(np.ceil((n_all - (int(n_ok.min()) - 2))
                                        / E_STEP))
    ax.set_ylim(y_lo, n_all)
    ax.set_yticks(range(y_lo, n_all + 1, E_STEP))
    ax.set_yticks(range(y_lo, n_all + 1, E_STEP // 4), minor=True)
    ax.annotate("default", (beta0, 0), xycoords=ax.get_xaxis_transform(),
                xytext=(3.5, 3), textcoords="offset points", ha="left",
                va="bottom")
    if n_pass_now != at0:
        # the stored scan predates the benchmark table: show its count too
        ax.plot(beta0, n_pass_now, ls="none", marker=MARKER["flow"], ms=3.8,
                mfc="white", mec=C["flow"], mew=0.8, zorder=4)
        ax.annotate("present benchmark", (beta0, n_pass_now),
                    xytext=(4.5, 2.0), textcoords="offset points", ha="left",
                    va="bottom")
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.set_xlabel(r"Electronegativity weight $\beta$ (" + "Å²)")
    # two lines; main() puts the first one on the line of the y label of (d)
    ax.set_ylabel(f"Cases passed\n(of {n_all})", linespacing=E_LINESPACING)
    ax.yaxis.set_label_coords(-0.95 / ax_cm(ax)[0], 0.5)
    return dict(beta=beta, n_ok=n_ok, n_all=n_all, beta0=beta0, at0=at0,
                col=col, ylim=(y_lo, n_all))


# ---------------------------------------------------------------- panel (f)
FLOW = C["flow"]              # this work: colour ...
MK = MARKER["flow"]           # ... and marker, as in the other panels
D1 = 2.5                      # marker diameter (pt) of a single point
EDGE = 0.7                    # edge width of open markers (pt)
HALO = 0.5                    # white rim around all-open markers (pt)
LEG_N = 5                     # node size shown in the legend of marker kinds
LEG_SHARE = 2 / 3             # filled share of the two-kind legend marker
# Hairline grid at the integer oxidation states.  Off: the style of the
# figure set allows gridlines only where they are essential, and the integer
# ticks on all four sides locate the points.
GRID_F = False


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
    if not n_same:
        # white rim: keeps an open marker apart from a larger node it touches
        ax.scatter([x], [y], s=(d + 2 * HALO) ** 2, marker=MK,
                   facecolor="white", edgecolor="none", zorder=z - 0.001)
    ax.scatter([x], [y], s=(d - EDGE) ** 2, marker=MK, facecolor="white",
               edgecolor=FLOW, linewidths=EDGE, zorder=z)
    if n_same:
        verts, ext = _sector(n_same / n)
        ax.scatter([x], [y], s=(d * ext) ** 2, marker=verts, facecolor=FLOW,
                   edgecolor="none", zorder=z + 0.001)


def pt_label(p, star=False, wrap=False):
    name = formula(p["formula"]) + ("*" if star else "")
    if p["tag"] and not star:
        name += f", {p['tag']}"
    if p["formula"] == p["species"]:
        return name
    return f"{name}{chr(10) if wrap else ' '}({p['species']})"


# Panel (f): label of every off-diagonal point, (id, species) -> placement
# (offsets in points from the point; an interval is labelled at its lower end).
LEADER = dict(arrowstyle="-", lw=0.4, color=C["muted"], shrinkA=0.5,
              shrinkB=2.5)
PLACE = {
    (7, "Au"): dict(xytext=(-9, 5), ha="left", va="bottom", wrap=True),
    (7, "Ca"): dict(xytext=(6, 0), ha="left", va="center"),
    (8, "K"): dict(xytext=(6, 0), ha="left", va="center"),
    (28, "Xe"): dict(xytext=(-6, 0), ha="right", va="center"),
    (28, "O"): dict(xytext=(6, -0.5), ha="left", va="center"),
    (29, "Xe"): dict(xytext=(-6, 0), ha="right", va="center"),
    (29, "O"): dict(xytext=(4, -8), ha="left", va="top"),
    (35, "Hg"): dict(xytext=(5, -5), ha="left", va="top"),
    (35, "F"): dict(xytext=(-4, 30), ha="right", va="bottom",
                    arrowprops=LEADER),
}


def panel_f(ax, points):
    lim = (-3.7, 6.7)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal", adjustable="box")
    ticks = list(range(-3, 7))
    for t in ticks if GRID_F else ():                # hairline integer grid
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
    for (x, y), (ns, no) in sorted(nodes.items(), key=lambda kv: -sum(kv[1])):
        draw_node(ax, x, y, ns, no)
    ebar = None
    for p in bars:
        x0, x1 = float(p["x"][0]), float(p["x"][1])
        ebar = ax.errorbar([(x0 + x1) / 2], [float(p["y"])],
                           xerr=[(x1 - x0) / 2], fmt="none", ecolor=FLOW,
                           elinewidth=1.0, capsize=2.2, capthick=0.8, zorder=5)

    off = [p for p in points if p["state"] == "differ"]
    missing = [(p["id"], p["species"]) for p in off
               if (p["id"], p["species"]) not in PLACE]
    assert not missing, f"off-diagonal points without a label slot: {missing}"
    for p in off:
        kw = dict(PLACE[(p["id"], p["species"])])
        wrap = kw.pop("wrap", False)
        x = float(p["x"][0] if isinstance(p["x"], tuple) else p["x"])
        # HgF4*: the one HgF4 cell off the diagonal (cell b, with F2 units)
        ax.annotate(pt_label(p, star=(p["id"] == 35), wrap=wrap),
                    xy=(x, float(p["y"])), textcoords="offset points",
                    zorder=6, linespacing=1.05, **kw)

    # Legends: kind of comparison (upper left), marker-area key (lower right)
    d_leg = D1 * np.sqrt(LEG_N)

    def ring():
        return Line2D([], [], ls="none", marker=MK, ms=d_leg - EDGE,
                      mfc="white", mec=FLOW, mew=EDGE)

    h_same = Line2D([], [], ls="none", marker=MK, ms=d_leg, mfc=FLOW,
                    mec="white", mew=0.3)
    handles = [h_same, ring()]
    # the first two entries have short lines: the labels of the two XeO3 (Xe)
    # points start to the right of this legend, on the same rows
    labels = ["same\nstructure", "different or\nunmatched\nstructure"]
    if any(ns and no for ns, no in nodes.values()):
        # a node holding both kinds, drawn as draw_node draws it
        verts, ext = _sector(LEG_SHARE)
        handles.append((ring(), Line2D([], [], ls="none", marker=verts,
                                       ms=d_leg * ext, mfc=FLOW, mec="none")))
        labels.append("both kinds coincide,\nsectors in proportion")
    if ebar is not None:
        handles.append(ebar)
        labels.append("interval implied by\nsite-resolved values")
    leg1 = ax.legend(handles, labels, loc="upper left", handletextpad=0.5,
                     labelspacing=0.5, borderaxespad=0.5,
                     handler_map={tuple: HandlerTuple(ndivide=1, pad=0)})
    ax.add_artist(leg1)
    sizes = (1, 5, 15)
    hs = [Line2D([], [], ls="none", marker=MK, ms=D1 * np.sqrt(n), mfc=FLOW,
                 mec="white", mew=0.3) for n in sizes]
    ax.legend(hs, [str(n) for n in sizes], title="Coincident points",
              title_fontsize=plt.rcParams["legend.fontsize"],
              loc="lower right", ncol=3, columnspacing=1.3,
              handletextpad=0.35, borderaxespad=0.5)


# --------------------------------------------------------------- provenance
def fmt_frac(v):
    if isinstance(v, tuple):
        return f"({fmt_frac(v[0])},{fmt_frac(v[1])})"
    if v.denominator == 1:
        return f"{int(v):+d}"
    return f"{float(v):+.3f}"


def print_abc(prov_a, prov_b, n_both, counts_b, bench, rows_c):
    names = {m: SHORT[m] for m in ("flow", "nearest")}
    print(f"(a) {rel(TSV['flow'])} ({names['flow']}), "
          f"{rel(TSV['nearest'])} ({names['nearest']})")
    print(f"    verdict: pass dev <= {PASS_MAX}, intermediate <= {MID_MAX}, "
          "fail above; status 'run' rows excluded")
    print(f"    {'group':12s} {'method':20s} {'n':>3s} {'pass':>5s} "
          f"{'int.':>5s} {'fail':>5s}")
    for name, method, n, c in prov_a:
        print(f"    {name:12s} {names[method]:20s} {n:3d} {c[0]:5d} "
              f"{c[1]:5d} {c[2]:5d}")
    docs = docs_numbers()
    print(f"    cross-check against {rel(RESULTS_MD)}:")
    for name, _ in GROUPS:
        mine = (next(n for g, m, n, c in prov_a if g == name),
                next(c[0] for g, m, n, c in prov_a if g == name and m == "flow"),
                next(c[0] for g, m, n, c in prov_a
                     if g == name and m == "nearest"))
        quoted = docs.get(name)
        verdict = ("not found in docs" if quoted is None
                   else "match" if quoted == mine else "MISMATCH")
        print(f"      {name:12s} n / flow passes / greedy passes: data {mine},"
              f" docs {quoted} -> {verdict}")

    print(f"\n(b) same two tables; {len(prov_b)} cases plotted, "
          f"{n_both} cases pass with both methods (not plotted); "
          f"passed by one method only / the other only / neither: {counts_b}")
    only = tuple(counts_b[:2])
    print(f"    only flow passes / only greedy passes: data {only}, docs "
          f"{docs.get('only')} -> "
          f"{'match' if docs.get('only') == only else 'MISMATCH'}")
    titles = []
    for t, *_ in prov_b:
        if t not in titles:
            titles.append(t)
    only2 = tuple(sum(1 for t, (b, _), _, _ in prov_b
                      if t == titles[k] and b == "wave2") for k in (0, 1))
    print(f"    of which Extension 2: data {only2}, docs "
          f"{docs.get('only_wave2')} -> "
          f"{'match' if docs.get('only_wave2') == only2 else 'MISMATCH'}")
    for title, (batch, case), df, dn in prov_b:
        print(f"    {title:34s} {TAG[batch]:3s} {batch:9s}{case:22s} "
              f"flow {df:5.2f}  greedy {dn:5.2f}")

    print(f"\n(c) frozen evaluator {rel(EVALUATOR)} (pyleos bench clean) "
          f"--split dev (cached in {rel(CACHE)}/)")
    for m in ("flow", "nearest", "geo"):
        s = bench[m]["summary"]
        print(f"    {LABEL[m].replace('$', ''):28s} atoms {s['atoms_ok']}/"
              f"{s['atoms_labelled']} "
              f"({100 * s['atoms_ok'] / s['atoms_labelled']:.1f}%)  cases with "
              f"all labelled atoms correct {s['all_atoms_correct_cases']}/"
              f"{s['cases']} ({100 * s['all_atoms_correct_cases'] / s['cases']:.1f}%)"
              f"  element-mean pass {s['elem_ok']}/{s['cases']}"
              f"  [--split dev {' '.join(BENCH_ARGS[m]) or '(default plugin)'}]")
    print("    min-cost flow, cases not fully correct (wrong sites located by "
          "re-running the frozen plugin;")
    print(f"    shared WF = Wannier centre whose two nearest atoms are of one "
          f"element, distance ratio >= {SHARED_RATIO}):")
    print(f"    {'case':10s} {'tag':3s} {'wrong':>7s}  {'el':2s} "
          f"{'d(X-X) min-max / A':>19s} {'shared WFs':>10s} "
          f"{'wrong sites on one':>18s}  evaluator detail")
    for r in rows_c:
        print(f"    {r['case']:10s} {TAG[r['batch']]:3s} "
              f"{r['n_wrong']:3d}/{r['n_labelled']:<3d}  "
              f"{'+'.join(r['elements']):2s} "
              f"{r['d_min']:9.3f}-{r['d_max']:<9.3f} {r['n_shared_wf']:10d} "
              f"{r['wrong_on_shared_wf']:8d}/{r['n_wrong']:<9d}  "
              f"{r['evaluator_wrong']}")
    n_wrong = sum(r["n_wrong"] for r in rows_c)
    s = bench["flow"]["summary"]
    bonded = [r for r in rows_c if r["wrong_on_shared_wf"] == r["n_wrong"]]
    print(f"    total wrong atoms {n_wrong} = {s['atoms_labelled']} - "
          f"{s['atoms_ok']}: {n_wrong == s['atoms_labelled'] - s['atoms_ok']}; "
          f"{len(bonded)}/{len(rows_c)} cases have every wrong site on a "
          f"like-atom shared WF ({sum(r['n_wrong'] for r in bonded)} atoms)")


def print_de(d, e, n_pass_now, beta_chk=None):
    groups, med = d["groups"], d["med"]
    print(f"\n(d) {rel(TSV['flow'])}: columns status, R_elem, n_stray; "
          "status 'run' rows excluded")
    print(f"    pass (OK) / intermediate or fail (~, FAIL): "
          f"{len(groups[0])} / {len(groups[1])}")
    print(f"    median R_elem: {med[0]:.2f} / {med[1]:.2f}")
    for gi, lab in enumerate(("pass", "intermediate or fail")):
        r_el = [float(r["R_elem"]) for r in groups[gi]]
        q1, q3 = np.percentile(r_el, [25, 75])
        print(f"    {lab}: quartiles {q1:.2f} - {q3:.2f}, range "
              f"{min(r_el):.1f} - {max(r_el):.1f}, cases with >= 1 stray WF "
              f"{d['n_stray'][gi]}")
    for r in sorted(groups[1], key=lambda r: float(r["R_elem"])):
        star = " (*)" if (r["batch"], r["case"]) in CONVENTION else ""
        print(f"    not passing: {TAG[r['batch']]} {r['batch']}/{r['case']:<22s}"
              f" R_elem {r['R_elem']:>5s}  dev {r['dev']}  status {r['status']:4s}"
              f" n_stray {r['n_stray']}{star}")
    for r in d["low_ok"]:
        print(f"    passing, R_elem < {med[1]:.2f}: {r['batch']}/{r['case']}"
              f"  R_elem {r['R_elem']}")
    for grp in d["tied"]:
        print(f"    bracketed (R_elem within {R_TIE:g} of the next), top to "
              "bottom: " + ", ".join(
                  f"{TAG[r['batch']]} {r['case']} {r['R_elem']}"
                  for r in reversed(grp)))
    worst_clean = max((float(r["R_elem"]) for r in groups[1]
                       if int(r["n_stray"]) == 0), default=float("nan"))
    print(f"    every non-passing case has >= 1 stray WF or R_elem <= "
          f"{worst_clean:g}: "
          f"{all(int(r['n_stray']) > 0 or float(r['R_elem']) <= worst_clean for r in groups[1])}")

    print(f"\n(e) {rel(BETA_TSV)}: columns beta, n_OK")
    print("    beta  " + " ".join(f"{b:>4g}" for b in e["beta"]))
    print("    n_OK  " + " ".join(f"{n:>4d}" for n in e["n_ok"]))
    print(f"    scored cases per beta: {e['n_all']}; default beta (column "
          f"name {e['col']}): {e['beta0']:g}; n_OK there: {e['at0']}")
    code_default = inspect.signature(assign_case).parameters["beta"].default
    print(f"    default beta of pyleos.transport.assign_case: {code_default}"
          f" -> {'match' if code_default == e['beta0'] else 'MISMATCH'}")
    for r in read_rows(BETA_TSV):
        if r[e["col"]] != "-":
            print(f"    changes at beta = {r['beta']}: {r[e['col']]}")
    print(f"    open symbol: {n_pass_now} passes in {rel(TSV['flow'])} at the "
          f"default beta; the stored scan gives {e['at0']} "
          f"({'same data' if n_pass_now == e['at0'] else 'scan predates the current benchmark data'})")
    print(f"    y axis: {e['ylim'][0]} to {e['ylim'][1]} (all {e['n_all']} "
          f"scored cases); the scan spans {e['n_ok'].min()}-{e['n_ok'].max()}"
          f" = {100 * e['n_ok'].min() / e['n_all']:.1f}-"
          f"{100 * e['n_ok'].max() / e['n_all']:.1f}% of the cases")
    cmd = f"python3 {rel(BETA_SCRIPT)} --out {rel(BETA_CHECK)}"
    if beta_chk is None:
        print(f"    check scan of the data tree: not present "
              f"(--check-beta runs: {cmd})")
    else:
        chk = np.array([int(r["n_OK"]) for r in beta_chk])
        same_beta = [float(r["beta"]) for r in beta_chk] == list(e["beta"])
        when = datetime.fromtimestamp(BETA_CHECK.stat().st_mtime)
        print(f"    check scan of the data tree, run {when:%Y-%m-%d %H:%M} "
              f"(NOT plotted; {cmd}):")
        print("      n_OK  " + " ".join(f"{n:>4d}" for n in chk))
        if same_beta:
            print("      minus the stored scan  "
                  + " ".join(f"{a - b:>+4d}" for a, b in zip(chk, e["n_ok"])))
            same_changes = ([r[e["col"]] for r in beta_chk]
                            == [r[e["col"]] for r in read_rows(BETA_TSV)])
            print(f"      status changes relative to the default identical to "
                  f"the stored scan: {same_changes}")
            at0 = int(chk[np.argmin(np.abs(e["beta"] - e["beta0"]))])
            print(f"      n_OK at the default: {at0}; equals the passes of "
                  f"{rel(TSV['flow'])}: {at0 == n_pass_now}")
        else:
            print("      beta values differ from the stored scan")

    print(f"\n    cross-check of (d), (e) against {rel(RESULTS_MD)} "
          "(sections 2, 7, 8):")
    got = {"n_ok": len(groups[0]), "n_not": len(groups[1]), "med_ok": med[0],
           "med_not": med[1], "beta10_scan": e["at0"]}
    for key, label in (("n_ok", "passing cases"), ("n_not", "non-passing cases"),
                       ("med_ok", "median R_elem, passing"),
                       ("med_not", "median R_elem, non-passing"),
                       ("beta10_scan", "n_OK at default beta (scan)")):
        ok = abs(float(got[key]) - float(QUOTED[key])) < 1e-6
        print(f"      {label:<30s} data {got[key]:g}   quoted {QUOTED[key]:g}"
              f"   {'match' if ok else 'MISMATCH'}")


def print_f(points, solids, leos):
    print("\n(f) literature value (Comas-Vila et al. 2025, table LIT) versus "
          f"LEOS element mean (pyleos defaults on cases/; cache "
          f"{rel(LEOS_CACHE)})")
    print("    struct: same / lattice (approximate lattice) / different / "
          "unmatched;   notes = line of docs/comas_vila_2025_notes.md")
    head = (f"    {'id':>2} {'solid':<18} {'LEOS dir (cases/)':<22} "
            f"{'struct':<9} {'sp':<3} {'lit':>10} {'LEOS':>7} "
            f"{'site range':>10} {'state':<7} {'stray':>5} {'ghost':>5} "
            f"{'notes':>5}  basis")
    print(head)
    print("    " + "-" * (len(head) - 4))
    by_id = defaultdict(list)
    for p in points:
        by_id[p["id"]].append(p)
    for s in solids:
        name = s["formula"] + (f" [{s['tag']}]" if s["tag"] else "")
        rows = by_id.get(s["id"], [])
        if not rows:
            st = element_stats(leos[s["rel"]])
            own = ", ".join(f"{el} {fmt_frac(v[0])}" for el, v in st.items())
            print(f"    {s['id']:>2} {name:<18} {s['rel']:<22} {s['struct']:<9} "
                  f"{'-':<3} {'none':>10} {'':>7} {'':>10} {'omitted':<7} "
                  f"{s['n_stray']:>5} {s['ghost']:>5} {'N' + str(s['nline']):>5}"
                  f"  {s['basis']}  [LEOS: {own}]")
        for k, p in enumerate(rows):
            lead = (f"    {s['id']:>2} {name:<18} {s['rel']:<22} {s['struct']:<9}"
                    if k == 0 else f"    {'':>2} {'':<18} {'':<22} {'':<9}")
            tail = (f"{s['n_stray']:>5} {s['ghost']:>5} "
                    f"{'N' + str(s['nline']):>5}  {s['basis']}" if k == 0 else "")
            rng = (f"{p['ymin']:+d}..{p['ymax']:+d}"
                   if p["ymin"] != p["ymax"] else "")
            print(f"{lead} {p['species']:<3} {fmt_frac(p['x']):>10} "
                  f"{fmt_frac(p['y']):>7} {rng:>10} {p['state']:<7} {tail}")
    print("\n    structure caveats")
    for s in solids:
        if s["rel"] in STRUCT_NOTE:
            print(f"      {s['rel']:<20} {STRUCT_NOTE[s['rel']]}")

    same = [p for p in points if p["comparable"]]
    other = [p for p in points if not p["comparable"]]
    print("\n    points (species per solid)")
    print(f"      plotted                        {len(points)}  "
          f"({len({p['id'] for p in points})} solids)")
    print(f"      same structure, agree          "
          f"{sum(p['state'] == 'agree' for p in same)}")
    print(f"      same structure, differ         "
          f"{sum(p['state'] == 'differ' for p in same)}  (of which intervals: "
          f"{sum(isinstance(p['x'], tuple) for p in same)})")
    print(f"      other structure, on diagonal   "
          f"{sum(p['state'] == 'agree' for p in other)}")
    print(f"      other structure, off diagonal  "
          f"{sum(p['state'] == 'differ' for p in other)}")
    print(f"      omitted (no literature value)  "
          f"{', '.join(s['formula'] for s in solids if not s['has_lit'])}")
    print("    solids per class  (agree / differ / not comparable)")
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
        print(f"      {label:<24} {cnt[0]:>2} / {cnt[1]} / {cnt[2]}   {members}")
    print(f"      {'total':<24} {tot['agree']:>2} / {tot['differ']} / {tot['nc']}")
    stray = {s["rel"].split("/")[0]: s["n_stray"] for s in solids
             if s["n_stray"]}
    print("    stray-WF warnings: " + ", ".join(
        f"{k} {v}" for k, v in sorted(stray.items(), key=lambda kv: -kv[1])))


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

    print("\n    structure facts read from the POSCAR files")
    st = parse_poscar(DATA / "k_hp/v3/POSCAR")
    print(f"      k_hp     {len(st.elements)} K atoms in the cell "
          f"(the paper's cell: 4 K + 2 non-nuclear attractors, p. 7079)")
    for tag in "abc":
        st = parse_poscar(DATA / f"hgf4_{tag}/v3/POSCAR")
        n_hg = st.elements.count("Hg")
        hf = dists(st, "Hg", "F")
        near = hf[hf < 2.6]
        ff = dists(st, "F", "F")
        print(f"      hgf4_{tag}   {len(st.elements)} atoms; Hg-F below 2.6 A: "
              f"{len(near) // n_hg} per Hg, {near.min():.2f}-{near.max():.2f} A;"
              f" shortest F-F {ff.min():.2f} A")


def cross_check_f(leos):
    """Compare with docs/RESULTS.md sec. 4 and results/bench74_golden.tsv."""
    print("\n    cross-checks of (f)")
    bad = 0
    for (case, el), ref in RESULTS_MD_MEANS.items():
        mean = float(element_stats(leos[f"{case}/v3"])[el][0])
        if abs(mean - ref) > TOL:
            bad += 1
            print(f"      MISMATCH {case} {el}: computed {mean:+.3f}, "
                  f"RESULTS.md {ref:+.3f}")
    print(f"      docs/RESULTS.md sec. 4 element means: "
          f"{len(RESULTS_MD_MEANS) - bad}/{len(RESULTS_MD_MEANS)} reproduced")
    got = {r.split("/")[0]: v["n_stray"] for r, v in leos.items()
           if INDEX[r]["role"] == "comas" and v["n_stray"]}
    print("      docs/RESULTS.md sec. 4 stray counts (12 of 30 cases): "
          + ("reproduced" if got == RESULTS_MD_STRAY else f"MISMATCH {got}"))
    tsv = {(r["batch"], r["case"]): r for r in read_rows(TSV["flow"])}
    for rel_dir, rec in leos.items():
        row = INDEX[rel_dir]
        if row["role"] != "bench74":
            continue
        case = row["ref_key"]
        q = tsv[(row["bench_set"], case)]
        ok = (int(q["n_stray"]) == rec["n_stray"]
              and int(q["n_ghost"]) == rec["ghost"])
        print(f"      bench74_golden.tsv {case:<10} status {q['status']:<4} "
              f"dev {q['dev']}  n_stray/n_ghost {q['n_stray']}/{q['n_ghost']} "
              f"{'= computed' if ok else '!= computed'}")


# --------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--recompute", action="store_true",
                    help="ignore the caches: re-run the evaluator and pyleos")
    ap.add_argument("--check-beta", action="store_true",
                    help="repeat the beta scan on the current data into "
                         "_cache/ (a few minutes; printed, never plotted)")
    args = ap.parse_args()

    tables = {m: read_tsv(p) for m, p in TSV.items()}
    if set(tables["flow"]) != set(tables["nearest"]):
        raise ValueError("the two benchmark tables list different cases")
    rel_rows = read_rows(TSV["flow"])
    beta_rows = read_rows(BETA_TSV)
    beta_chk = beta_check(args.check_beta)
    bench = {m: clean_bench(m, args.recompute)
             for m in ("flow", "nearest", "geo")}
    contacts = like_atom_contacts(bench["flow"], args.recompute)
    assert len(LIT) == 40 and len({r[1] for r in LIT}) == 40
    leos = leos_values(args.recompute)
    points, solids = build_points(leos)
    n_pass_now = sum(s == "OK" for s, _ in tables["flow"].values())

    # --- layout in cm from the left / top edges; the saved file is exactly
    # DOUBLE wide and height_cm tall
    pad = 0.02                                   # pad_inches of style.save
    top, h_ac, h_x, band = 0.95, 3.30, 0.92, 0.45
    left, right = 1.95, 16.90
    h_b = 2 * h_ac + h_x + band                  # (b) spans (a) and (c)
    top_c = top + h_ac + h_x + band
    top_2 = top + h_b + h_x + band               # second block: (d), (e), (f)
    side_f = 6.60                                # (f) is square
    h_e, h_xd = 1.95, 0.75
    h_d = side_f - h_e - h_xd - band
    height_cm = top_2 + side_f + 1.00       # room for the x labels
    fig = plt.figure(figsize=(style.DOUBLE - 2 * pad,
                              height_cm * style.CM - 2 * pad))
    # invisible frame: pins the tight bounding box to the full figure
    fig.add_artist(Rectangle((0, 0), 1, 1, transform=fig.transFigure,
                             fill=False, edgecolor="none", lw=0))
    ax_a = cm_axes(fig, left, top, 5.45, h_ac)
    ax_c1 = cm_axes(fig, left, top_c, 1.45, h_ac)
    ax_c2 = cm_axes(fig, 5.97, top_c, 1.15, h_ac)
    ax_b = cm_axes(fig, 12.10, top, right - 12.10, h_b)
    w_de = 6.85
    ax_d = cm_axes(fig, left, top_2, w_de, h_d)
    ax_e = cm_axes(fig, left, top_2 + side_f - h_e, w_de, h_e)
    ax_f = cm_axes(fig, right - side_f, top_2, side_f, side_f)

    prov_a = panel_a(ax_a, tables)
    prov_b, n_both, counts_b = panel_b(ax_b, tables)
    panel_c_summary(ax_c1, bench)
    rows_c = panel_c_cases(ax_c2, contacts)
    res_d = panel_d(ax_d, rel_rows)
    res_e = panel_e(ax_e, beta_rows, n_pass_now)
    align_ylabel(ax_e, ax_d)
    panel_f(ax_f, points)

    # panel labels on one grid: two columns (0.45 and 10.15 cm from the left
    # edge, right-aligned), each 0.10 cm above the top of its axes
    fw = fig.get_size_inches()[0] / style.CM
    for ax, letter, x_cm in ((ax_a, "a", 0.45), (ax_b, "b", 10.15),
                             (ax_c1, "c", 0.45), (ax_d, "d", 0.45),
                             (ax_e, "e", 0.45), (ax_f, "f", 10.15)):
        w, hh = ax_cm(ax)
        style.panel_label(ax, letter,
                          dx=(x_cm - ax.get_position().x0 * fw) / w,
                          dy=1 + 0.10 / hh)
    method_legend(fig, left, 0.27, ("flow", "nearest", "geo"))
    style.save(fig, "fig2_benchmark")

    # ------------------------------------------------------- provenance
    print("\nFigure 2 provenance")
    print("=" * 72)
    print_abc(prov_a, prov_b, n_both, counts_b, bench, rows_c)
    print_de(res_d, res_e, n_pass_now, beta_chk)
    print_f(points, solids, leos)
    structure_facts()
    cross_check_f(leos)
    print(f"\nwrote {HERE / 'fig2_benchmark.pdf'} and .png "
          f"({style.DOUBLE / style.CM:.1f} x {height_cm:.2f} cm)")


if __name__ == "__main__":
    main()
