#!/usr/bin/env python3
"""Figure 4 - Delithiation case study: who is oxidised when Li is removed.

(a) Li2O2 delithiation: mean oxidation state of O and Li.
(b) Where the hole goes in Li15O16: charge of every O2 unit against its O-O
    distance, for six calculations (PBE+U scan with U on O 2p, then PBE and
    HSE06 single points on the U = 8 eV geometry).
(c) Li-rich cathodes along 50 ps MD trajectories (20 frames each).

Layout: (a) above (b) in the left column, (c) over the full height of the
right column; the saved page is exactly 17.1 cm wide.

Run from this directory:   python3 fig4_delithiation.py [--recompute]

Every plotted number is read from, or computed with pyleos from, files under
cases/ and results/.  The baseline ("per-WF greedy rule", style key
"nearest") is the per-atom table that the VASP patch writes to each EOSCAR:
its OS(int) column for element means, its n_e(assigned) column for the charge
of an O2 unit.
pyleos results are cached in _cache/fig4_*.json together with the md5 of
every input file and of the pyleos sources; a changed input (or --recompute)
triggers recomputation.  Writes fig4_delithiation.pdf / .png and prints a
provenance table.
"""
import argparse
import csv
import hashlib
import json
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

import style
from style import C, MARKER, LABEL, formula

import matplotlib.pyplot as plt                     # noqa: E402
from matplotlib.lines import Line2D                 # noqa: E402
from matplotlib.patches import Rectangle            # noqa: E402
from matplotlib.ticker import FixedLocator, MultipleLocator  # noqa: E402

REPO = style.REPO
DATA = REPO / "cases"                       # <case>/<run>/ directories
CACHE = style.HERE / "_cache"
sys.path.insert(0, str(REPO))
from pyleos.parsers import (CaseData, Structure, WannierChannel,  # noqa: E402
                            load_case, parse_eoscar, parse_poscar)
from pyleos.transport import assign_case, reliability_margins  # noqa: E402

PYLEOS_SRC = [REPO / "pyleos" / f for f in
              ("parsers.py", "transport.py", "capacity.py")]
CACHE_VERSION = 2           # bump when the layout of the cached data changes
BASE = "nearest"            # style.py key of the per-WF greedy baseline

# ---------------------------------------------------------------- inputs --
LI2O2_RUN = "v31"           # cases/<name>/v31: Foppl cells (the v3 cells put 8 Li at O-O bond midpoints)
LI2O2_CASES = ["li16o16", "li15o16", "li14o16"]

# Li15O16 with one pre-shortened peroxide unit: run name and the Hubbard U
# on O 2p (eV) / hybrid flag that its INCAR must carry (checked on load).
# Run <name> is the directory cases/li15o16/polaron-<name> ("_" -> "-").
POLARON_DIR = DATA / "li15o16"
POLARON = [("u0", 0, False), ("u4", 4, False), ("u6", 6, False),
           ("u8", 8, False), ("sp_pbe", 0, False), ("sp_hse", 0, True)]
POLARON_GEOM = "u8"         # sp_pbe and sp_hse are single points on this one
PAIR_CUTOFF = 2.4           # Angstrom: longest O-O distance still paired
MD_DIR = DATA / "li15o16" / "v3-md2"    # 20 analysed steps (not plotted)

TS_DIR = REPO / "results" / "cathode_md"     # timeseries_<system>_x<x>.tsv
STRUCT_DIR = DATA                            # <system>/md-x<x>/POSCAR
X_LEVELS = ["0.00", "0.25", "0.50"]
# (system directory, element, label of the system)
TS_ROWS = [("li2mno3", "Mn", "Li2MnO3"), ("li2ruo3", "Ru", "Li2RuO3"),
           ("nmc", "Ni", "NMC"), ("nmc", "Mn", "NMC"), ("nmc", "Co", "NMC")]

PROV = []       # (panel, quantity, value, source)


def prov(panel, quantity, value, source):
    PROV.append((panel, quantity, str(value), str(source)))


def rel(p):
    return str(Path(p).relative_to(REPO))


def li2o2_dir(name):
    return DATA / name / LI2O2_RUN


def polaron_dir(name):
    return POLARON_DIR / ("polaron-" + name.replace("_", "-"))


# ------------------------------------------------------------ cache layer --
def md5(path):
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def cached(name, inputs, compute, force):
    """Return compute() from _cache/fig4_<name>.json when the md5 of every
    input (and of the pyleos sources) is unchanged; recompute otherwise."""
    path = CACHE / f"fig4_{name}.json"
    sig = {rel(p): md5(p) for p in list(inputs) + PYLEOS_SRC}
    sig["cache_version"] = CACHE_VERSION
    if path.is_file() and not force:
        blob = json.loads(path.read_text())
        if blob.get("inputs") == sig:
            return blob["data"]
    t0 = time.time()
    print(f"[fig4] computing {name} with pyleos ...", flush=True)
    data = compute()
    CACHE.mkdir(exist_ok=True)
    path.write_text(json.dumps(
        {"generated_by": "fig4_delithiation.py", "inputs": sig, "data": data},
        indent=1))
    print(f"[fig4] {name}: {time.time() - t0:.0f} s -> {rel(path)}")
    return data


# --------------------------------------------------------------- geometry --
SHIFTS = np.array([[a, b, c] for a in (-1, 0, 1) for b in (-1, 0, 1)
                   for c in (-1, 0, 1)])


def min_image_dist(st, i, j):
    lat = np.array(st.lattice)
    d = np.array(st.frac[j]) - np.array(st.frac[i])
    d -= np.round(d)
    return float(np.min(np.linalg.norm((d + SHIFTS) @ lat, axis=1)))


def idx(elements, el):
    return [i for i, e in enumerate(elements) if e == el]


def o2_units(st, ox, cutoff=PAIR_CUTOFF):
    """Pair O atoms by repeatedly taking the shortest remaining O-O distance
    (up to `cutoff`).  Returns (pairs [(i, j, d)], unpaired atoms)."""
    dist = sorted((min_image_dist(st, i, j), i, j)
                  for a, i in enumerate(ox) for j in ox[a + 1:])
    free, pairs = set(ox), []
    for d, i, j in dist:
        if d > cutoff:
            break
        if i in free and j in free:
            free -= {i, j}
            pairs.append((i, j, d))
    return pairs, sorted(free)


# ----------------------------------------------------- panel (a): Li2O2 ---
def compute_li2o2():
    out = {}
    for name in LI2O2_CASES:
        case = load_case(li2o2_dir(name))
        res = assign_case(case)                       # beta = 10, Allen
        st, el = case.structure, case.structure.elements
        li, ox = idx(el, "Li"), idx(el, "O")
        pairs, _ = o2_units(st, ox)
        g = case.eoscar.greedy_os                     # EOSCAR OS(int) column
        out[name] = {
            "n_li": len(li), "n_o": len(ox),
            "flow": {"Li": [res.os_int[i] for i in li],
                     "O": [res.os_int[i] for i in ox]},
            "greedy": {"Li": [g[i] for i in li], "O": [g[i] for i in ox]},
            "o2_unit_flow": [res.os_int[i] + res.os_int[j]
                             for i, j, _ in pairs],
            "ghost": res.ghost, "stray": len(res.stray),
            "dne": case.eoscar.dne,
        }
    return out


def load_li2o2(force):
    inputs = [li2o2_dir(n) / f for n in LI2O2_CASES
              for f in ("EOSCAR", "POSCAR")]
    raw = cached("li2o2", inputs, compute_li2o2, force)
    n0 = raw[LI2O2_CASES[0]]["n_li"]
    rows = []
    for name in LI2O2_CASES:
        d = raw[name]
        row = {"name": name, "removed": n0 - d["n_li"],
               # charge counting: Li = +1, charge-neutral cell
               "ref": {"Li": 1.0, "O": -d["n_li"] / d["n_o"]},
               "flow": {e: float(np.mean(d["flow"][e])) for e in ("Li", "O")},
               BASE: {e: float(np.mean(d["greedy"][e])) for e in ("Li", "O")}}
        rows.append(row)
        src = rel(li2o2_dir(name))
        prov("a", f"{name}: mean OS(O) ref / flow / greedy",
             f"{row['ref']['O']:+.4f} / {row['flow']['O']:+.4f} / "
             f"{row[BASE]['O']:+.4f}", src + "/{POSCAR,EOSCAR}")
        prov("a", f"{name}: mean OS(Li) ref / flow / greedy",
             f"{row['ref']['Li']:+.4f} / {row['flow']['Li']:+.4f} / "
             f"{row[BASE]['Li']:+.4f}  "
             f"(greedy: {d['greedy']['Li'].count(-1)} Li at -1)", src)
        prov("a", f"{name}: O2-unit charges, flow (not plotted)",
             sorted(d["o2_unit_flow"]), "pyleos, site-degenerate")
        prov("a", f"{name}: ghost e / stray WF / |dNe|",
             f"{d['ghost']} / {d['stray']} / {d['dne']:.2e}", src + "/EOSCAR")
    return rows


# ------------------------------------------ panel (b): hole in Li15O16 ----
def incar_tags(path):
    tags = {}
    for line in Path(path).read_text().splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            tags[k.strip().upper()] = v.split("#")[0].strip()
    return tags


def compute_polaron():
    out = {}
    for name, _, _ in POLARON:
        d = polaron_dir(name)
        case = load_case(d)
        res = reliability_margins(case, res=assign_case(case))
        st, el = case.structure, case.structure.elements
        ox, li = idx(el, "O"), idx(el, "Li")
        z = case.eoscar.zval["O"]
        ne, g = case.eoscar.greedy_ne, case.eoscar.greedy_os
        pairs, free = o2_units(st, ox)
        out[name] = {
            "units": [{"atoms": [i + 1, j + 1], "d": dist,
                       "flow": res.os_int[i] + res.os_int[j],
                       # greedy: electrons counted on the two atoms
                       "greedy": 2 * z - ne[i] - ne[j],
                       "greedy_int": g[i] + g[j]} for i, j, dist in pairs],
            "unpaired": [{"atom": i + 1, "flow": res.os_int[i],
                          "greedy": z - ne[i], "greedy_int": g[i]}
                         for i in free],
            "o_mean_flow": float(np.mean([res.os_int[i] for i in ox])),
            "o_mean_greedy": float(np.mean([g[i] for i in ox])),
            "li_flow": sorted({res.os_int[i] for i in li}),
            "li_greedy_not_plus1": [[i + 1, g[i], round(ne[i], 3)]
                                    for i in li if g[i] != 1],
            "ghost": res.ghost, "stray": len(res.stray),
            "dne": case.eoscar.dne,
            "n_wf": [len(ch.centres) for ch in case.channels],
            "r_site": res.case_conf, "r_elem": res.case_conf_elem,
            "min_site_margin": min(res.atom_margin[i] for i in ox),
        }
    return out


def polaron_label(name, tags, same_geom):
    if tags.get("LHFCALC", "").upper().startswith(".T"):
        level = "HSE06"
    else:
        level = "PBE"
    u = float(tags.get("LDAUU", "0 0").split()[-1])
    if same_geom and name != POLARON_GEOM:        # single point on U = 8
        return f"{level} // $U$ = 8"
    return "$U$ = 0" if u == 0 else f"$U$ = {u:g} eV"


def load_polaron(force):
    inputs = [polaron_dir(n) / f for n, _, _ in POLARON
              for f in ("EOSCAR", "POSCAR")]
    raw = cached("polaron", inputs, compute_polaron, force)
    geom = np.array(parse_poscar(polaron_dir(POLARON_GEOM) / "POSCAR").frac)
    out = []
    for name, u_expect, hybrid in POLARON:
        d = polaron_dir(name)
        tags = incar_tags(d / "INCAR")
        u = float(tags.get("LDAUU", "0 0").split()[-1])
        is_hyb = tags.get("LHFCALC", "").upper().startswith(".T")
        assert (u, is_hyb) == (u_expect, hybrid), (name, u, is_hyb)
        frac = np.array(parse_poscar(d / "POSCAR").frac)
        same = bool(np.abs(frac - geom).max() < 1e-6)   # same geometry
        assert same == (name in (POLARON_GEOM, "sp_pbe", "sp_hse")), name
        r = raw[name]
        units = sorted(r["units"], key=lambda q: q["d"])
        out.append({"name": name, "label": polaron_label(name, tags, same),
                    "units": units})
        src = rel(d)
        for q in units:
            prov("b", f"{name}: O{q['atoms'][0]}-O{q['atoms'][1]} "
                 "d (A) / flow / greedy n_e / greedy OS(int)",
                 f"{q['d']:.3f} / {q['flow']:+d} / {q['greedy']:+.2f} / "
                 f"{q['greedy_int']:+d}", src + "/{POSCAR,EOSCAR}")
        for q in r["unpaired"]:
            prov("b", f"{name}: unpaired O{q['atom']} (not plotted) "
                 "flow / greedy n_e / greedy OS(int)",
                 f"{q['flow']:+d} / {q['greedy']:+.2f} / "
                 f"{q['greedy_int']:+d}", src)
        short = units[0]
        hole = [q for q in units if q["flow"] > -2]
        prov("b", f"{name}: shortest unit d / flow / greedy; units with "
             "flow charge above -2",
             f"{short['d']:.3f} / {short['flow']:+d} / {short['greedy']:+.2f}"
             "; " + ", ".join(f"O{q['atoms'][0]}-O{q['atoms'][1]} "
                              f"{q['d']:.3f} A {q['flow']:+d}" for q in hole),
             src)
        prov("b", f"{name}: mean OS(O) flow / greedy OS(int); Li (flow); "
             "greedy Li not +1 [atom, OS, n_e]",
             f"{r['o_mean_flow']:+.4f} / {r['o_mean_greedy']:+.4f}; "
             f"{r['li_flow']}; {r['li_greedy_not_plus1']}", src + "/EOSCAR")
        prov("b", f"{name}: R_site / R_elem / smallest O site margin "
             "(cost units)",
             f"{r['r_site']:.2f} / {r['r_elem']:.1f} / "
             f"{r['min_site_margin']:.4f}", "pyleos reliability_margins")
        prov("b", f"{name}: N_WF per spin / ghost e / stray / |dNe|",
             f"{r['n_wf'][0]},{r['n_wf'][1]} / {r['ghost']} / {r['stray']} / "
             f"{r['dne']:.1e}", src + "/EOSCAR")
    return out


# ------------------------------- Li15O16 MD, 20 steps (provenance only) ---
def read_xdatcar(path):
    x = Path(path).read_text().splitlines()
    lat = np.array([[float(v) for v in x[i].split()] for i in (2, 3, 4)])
    lat *= float(x[1])
    els = [e for e, n in zip(x[5].split(), x[6].split()) for _ in range(int(n))]
    conf, i = [], 7
    while i < len(x):
        if x[i].startswith("Direct"):
            conf.append([[float(v) for v in x[i + 1 + k].split()[:3]]
                         for k in range(len(els))])
            i += 1 + len(els)
        else:
            i += 1
    return lat.tolist(), els, conf


def compute_md():
    """Flow per analysed step of the trajectory EOSCAR.  The step-to-geometry
    map is not recorded (analysis every second MD step), so both possible
    XDATCAR configurations are tried."""
    lines = (MD_DIR / "EOSCAR").read_text().splitlines()
    starts = [i for i, ln in enumerate(lines) if ln.startswith("# STEP")]
    lat, els, conf = read_xdatcar(MD_DIR / "XDATCAR")
    ox = idx(els, "O")
    out = []
    with tempfile.TemporaryDirectory() as tmp:
        for s, (a, b) in enumerate(zip(starts, starts[1:] + [len(lines)])):
            p = Path(tmp) / "EOSCAR"
            p.write_text("\n".join(lines[a:b]) + "\n")
            eo = parse_eoscar(p)
            chans = [WannierChannel([list(w.centre) for w in eo.wf_detail
                                     if w.spin == sp],
                                    [w.spread for w in eo.wf_detail
                                     if w.spin == sp]) for sp in (1, 2)]
            row = {"step": s + 1, "n_wf": [len(c.centres) for c in chans],
                   "greedy_o_mean": float(np.mean([eo.greedy_os[i]
                                                   for i in ox]))}
            for tag, k in (("odd", 2 * s), ("even", 2 * s + 1)):
                st = Structure(lat, list(els), conf[k])
                res = assign_case(CaseData(MD_DIR, st, chans, eo, 2))
                pairs, _ = o2_units(st, ox)
                row[tag] = {"o_mean": float(np.mean([res.os_int[i]
                                                     for i in ox])),
                            "units": [res.os_int[i] + res.os_int[j]
                                      for i, j, _ in pairs],
                            "ghost": res.ghost}
            out.append(row)
    bonds = [d for c in conf
             for _, _, d in o2_units(Structure(lat, list(els), c), ox)[0]]
    return {"n_conf": len(conf), "steps": out,
            "oo_range": [min(bonds), max(bonds)]}


def load_md(force):
    md = cached("li2o2_md", [MD_DIR / "EOSCAR", MD_DIR / "XDATCAR"],
                compute_md, force)
    st = md["steps"]
    means = sorted({round(r[t]["o_mean"], 4) for r in st
                    for t in ("odd", "even")})
    units = [u for r in st for u in r["odd"]["units"]]
    prov("b", "MD (not plotted): analysed steps / XDATCAR configurations / "
         "N_WF per spin / ghost e",
         f"{len(st)} / {md['n_conf']} / {st[0]['n_wf'][0]} / "
         f"{st[0]['odd']['ghost']}", rel(MD_DIR))
    prov("b", "MD (not plotted): flow mean OS(O), all steps, both geometry "
         "maps; greedy OS(int) mean min..max", f"{means}; "
         f"{min(r['greedy_o_mean'] for r in st):+.4f}.."
         f"{max(r['greedy_o_mean'] for r in st):+.4f}", rel(MD_DIR))
    prov("b", "MD (not plotted): flow O2-unit charges, min..max over "
         "20 steps x 8 units", f"{min(units):+d}..{max(units):+d}",
         rel(MD_DIR))
    prov("b", "MD (not plotted): O-O distance of the O2 units, min..max "
         "over all configurations (A)",
         f"{md['oo_range'][0]:.3f}..{md['oo_range'][1]:.3f}",
         rel(MD_DIR / "XDATCAR"))


# ----------------------------------------------- panel (c): time series ---
def read_tsv(path):
    lines = [ln for ln in Path(path).read_text().splitlines()
             if ln.strip() and not ln.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def poscar_counts(path):
    lines = Path(path).read_text().splitlines()
    return dict(zip(lines[5].split(), map(int, lines[6].split())))


def load_timeseries():
    out = {}
    for system in sorted({r[0] for r in TS_ROWS}):
        for x in X_LEVELS:
            f = TS_DIR / f"timeseries_{system}_x{x}.tsv"
            rows = read_tsv(f)
            counts = poscar_counts(STRUCT_DIR / system / f"md-x{x}" / "POSCAR")
            tm = [e for e in counts if e not in ("Li", "O")]
            n_tm = sum(counts[e] for e in tm)
            d = {"t": np.array([float(r["t_ps"]) for r in rows]),
                 "counts": counts,
                 # charge counting with Li+ and O2-, averaged over all TM atoms
                 "ref_tm": (2 * counts["O"] - counts["Li"]) / n_tm,
                 "stray": [int(r["n_stray"]) for r in rows],
                 "ghost": [int(r["n_ghost"]) for r in rows]}
            for col in rows[0]:
                if col.endswith(("_flow", "_greedy")):
                    d[col] = np.array([float(r[col]) for r in rows])
            d["tm_flow"] = sum(counts[e] * d[f"{e}_flow"] for e in tm) / n_tm
            out[(system, x)] = d
            comp = " ".join(f"{e}{n}" for e, n in counts.items())
            prov("c", f"{system} x={x}: cell / frames / t range (ps)",
                 f"{comp} / {len(rows)} / {d['t'][0]:.1f}-{d['t'][-1]:.1f}",
                 rel(f))
            prov("c", f"{system} x={x}: TM mean, charge counting / flow "
                 "(min..max over frames)",
                 f"{d['ref_tm']:.4f} / {d['tm_flow'].min():.4f}.."
                 f"{d['tm_flow'].max():.4f}",
                 rel(STRUCT_DIR / system / f"md-x{x}" / "POSCAR"))
            for e in tm:
                fl, gr = d[f"{e}_flow"], d[f"{e}_greedy"]
                prov("c", f"{system} x={x}: {e} flow mean+-sd [min,max] | "
                     "greedy mean+-sd [min,max]",
                     f"{fl.mean():.3f}+-{fl.std():.3f} [{fl.min():.2f},"
                     f"{fl.max():.2f}] | {gr.mean():.3f}+-{gr.std():.3f} "
                     f"[{gr.min():.2f},{gr.max():.2f}]", rel(f))
            lo = (d["Li_flow"].min(), d["Li_flow"].max(),
                  d["O_flow"].min(), d["O_flow"].max())
            prov("c", f"{system} x={x}: flow Li min/max, O min/max; greedy "
                 "O mean", f"{lo[0]:+.2f}/{lo[1]:+.2f}, {lo[2]:+.2f}/"
                 f"{lo[3]:+.2f}; {d['O_greedy'].mean():+.3f}", rel(f))
            prov("c", f"{system} x={x}: stray WF range / ghost e range",
                 f"{min(d['stray'])}-{max(d['stray'])} / "
                 f"{min(d['ghost'])}-{max(d['ghost'])}", rel(f))
    return out


# ================================================================ drawing ==
# All positions are in centimetres on the canvas.  The canvas is the page
# minus the pad that style.save adds on every side, so that the saved file
# is exactly PAGE_W x PAGE_H.
PAGE_W = style.DOUBLE / style.CM       # 17.1 cm
PAGE_H = 12.2                          # cm
PAD = 0.02 / style.CM                  # pad_inches of style.save, in cm
FIG_W = PAGE_W - 2 * PAD
FIG_H = PAGE_H - 2 * PAD
COL_L = (1.25, 5.00)                   # left column: axes left edge, width
COL_R = (7.75, 7.70)                   # right column: grid left edge, width
TOP = 0.92                             # top edge of (a) and (c)
A_H = (2.05, 1.30)                     # (a): heights of the O and Li axes
A_GAP = 0.15
B_TOP = 5.85                           # top edge of (b)
B_ROW, B_GAP = 1.70, 0.15              # (b): facet height, gap between facets
B_NCOL = 2
# (b): the O-O distance axis is broken into two windows drawn at one scale:
# bonded units and units elongated beyond 2.1 A (nothing lies in between)
B_SEG = [(1.34, 1.60), (2.11, 2.33)]   # Angstrom
B_BREAK = 0.10                         # cm, gap at the axis break
B_NROW = -(-len(POLARON) // B_NCOL)
BOTTOM = B_TOP + B_NROW * B_ROW + (B_NROW - 1) * B_GAP   # of (b) and (c)


def ax_cm(fig, left, top, width, height, **kw):
    """Axes placed in centimetres measured from the top-left corner."""
    return fig.add_axes([left / FIG_W, 1 - (top + height) / FIG_H,
                         width / FIG_W, height / FIG_H], **kw)


def label_panel(ax, letter):
    """style.panel_label at a fixed physical offset from the axes corner."""
    w, h = (ax.get_position().width * FIG_W, ax.get_position().height * FIG_H)
    style.panel_label(ax, letter, dx=-0.10 / w, dy=1 + 0.10 / h)


def tag(ax, text, x, y, ha="left", va="top"):
    ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va)


def ylabel_cm(fig, axes_left, top, height, text, offset=0.95):
    """Rotated y label centred on a block of stacked axes; its centre line
    sits `offset` cm left of the axes edge (same for every panel)."""
    fig.text((axes_left - offset) / FIG_W, 1 - (top + height / 2) / FIG_H,
             text, rotation=90, ha="center", va="center", fontsize=8)


REF_KW = dict(color=C["ref"], lw=0.7, ls=(0, (3, 2)))
SERIES = {  # zorder: reference < greedy < flow, sizes shrink accordingly
    "ref": dict(marker=MARKER["ref"], ms=6.0, mfc="none", mew=0.8, zorder=2,
                **REF_KW),
    BASE: dict(color=C[BASE], marker=MARKER[BASE], ms=4.0, mew=0, lw=1.0,
               zorder=3),
    "flow": dict(color=C["flow"], marker=MARKER["flow"], ms=3.0, mew=0,
                 lw=1.0, zorder=4),
}


def draw_a(fig, rows):
    ax_o = ax_cm(fig, COL_L[0], TOP, COL_L[1], A_H[0])
    ax_li = ax_cm(fig, COL_L[0], TOP + A_H[0] + A_GAP, COL_L[1], A_H[1],
                  sharex=ax_o)
    x = [r["removed"] for r in rows]
    for ax, el in ((ax_o, "O"), (ax_li, "Li")):
        for key in ("ref", BASE, "flow"):
            ax.plot(x, [r[key][el] for r in rows], **SERIES[key])
    ax_o.set_ylim(-1.13, 0.13)
    ax_o.yaxis.set_major_locator(MultipleLocator(0.5))
    ax_o.yaxis.set_minor_locator(MultipleLocator(0.25))
    ax_li.set_ylim(-0.17, 1.17)
    ax_li.yaxis.set_major_locator(MultipleLocator(0.5))
    ax_li.yaxis.set_minor_locator(MultipleLocator(0.25))
    ax_o.set_xlim(-0.35, 2.35)
    ax_o.xaxis.set_major_locator(FixedLocator(x))
    ax_o.tick_params(labelbottom=False)
    ax_li.set_xlabel("Li removed per " + formula("Li16O16") + " cell")
    tag(ax_o, "O", 0.04, 0.95)
    tag(ax_li, "Li", 0.04, 0.12, va="bottom")
    ylabel_cm(fig, COL_L[0], TOP, sum(A_H) + A_GAP, "Mean oxidation state")
    label_panel(ax_o, "a")


def draw_b(fig, pol):
    w = (COL_L[1] - B_GAP * (B_NCOL - 1)) / B_NCOL
    scale = (w - B_BREAK) / sum(b - a for a, b in B_SEG)   # cm per Angstrom
    ticks = [[1.4, 1.5], [2.2, 2.3]]
    first = None
    for k, calc in enumerate(pol):
        i, j = divmod(k, B_NCOL)
        d = np.array([q["d"] for q in calc["units"]])
        yg = [q["greedy"] for q in calc["units"]]
        yf = [q["flow"] for q in calc["units"]]
        assert all(any(a < x < b for a, b in B_SEG) for x in d), calc["name"]
        left = COL_L[0] + j * (w + B_GAP)
        seg = []
        for s, (a, b) in enumerate(B_SEG):
            ax = ax_cm(fig, left, B_TOP + i * (B_ROW + B_GAP),
                       (b - a) * scale, B_ROW)
            left += (b - a) * scale + B_BREAK
            # one thin connector per O2 unit: its two charges belong together
            ax.vlines(d, yg, yf, color=C["muted"], lw=0.5, zorder=2)
            ax.plot(d, yg, ls="none", color=C[BASE], marker=MARKER[BASE],
                    ms=3.3, mew=0, zorder=3)
            ax.plot(d, yf, ls="none", color=C["flow"], marker=MARKER["flow"],
                    ms=2.7, mew=0, zorder=4)
            ax.set_xlim(a, b)
            ax.set_ylim(-4.6, 2.0)
            ax.xaxis.set_major_locator(FixedLocator(ticks[s]))
            ax.xaxis.set_minor_locator(MultipleLocator(0.05))
            ax.yaxis.set_major_locator(FixedLocator([-4, -2, 0]))
            ax.yaxis.set_minor_locator(FixedLocator([-3, -1, 1]))
            ax.tick_params(labelbottom=(i == B_NROW - 1),
                           labelleft=(j == 0 and s == 0), pad=2)
            seg.append(ax)
        ax_l, ax_r = seg
        ax_l.spines["right"].set_visible(False)
        ax_r.spines["left"].set_visible(False)
        ax_l.tick_params(which="both", right=False)
        ax_r.tick_params(which="both", left=False)
        slash = dict(marker=[(-1, -2), (1, 2)], ms=4.5, ls="none", mew=0.6,
                     color=C["ink"], mec=C["ink"], clip_on=False, zorder=6)
        for ax, x in ((ax_l, 1), (ax_r, 0)):       # axis-break marks
            ax.plot([x, x], [0, 1], transform=ax.transAxes, **slash)
        ax_r.text(0.90, 0.93, calc["label"], transform=ax_r.transAxes,
                  ha="right", va="top")
        first = first or ax_l
    height = BOTTOM - B_TOP
    ylabel_cm(fig, COL_L[0], B_TOP, height, "Charge of " + formula("O2")
              + " unit")
    fig.text((COL_L[0] + COL_L[1] / 2) / FIG_W,
             1 - (BOTTOM + 0.56) / FIG_H, "O–O distance (Å)", ha="center",
             va="center", fontsize=8)
    label_panel(first, "b")


def draw_c(fig, ts):
    n_col = len(X_LEVELS)
    cgap, rgap = 0.15, 0.15
    cw = (COL_R[1] - cgap * (n_col - 1)) / n_col
    height = BOTTOM - TOP
    h_ord = 1.55                                    # ordered systems
    n_ord = sum(1 for r in TS_ROWS if r[0] != "nmc")
    n_nmc = len(TS_ROWS) - n_ord
    h_nmc = (height - rgap * (len(TS_ROWS) - 1) - h_ord * n_ord) / n_nmc
    top = TOP
    first = None
    for i, (system, el, name) in enumerate(TS_ROWS):
        ordered = system != "nmc"
        h = h_ord if ordered else h_nmc
        vals = np.concatenate([ts[(system, x)][f"{el}_{m}"] for x in X_LEVELS
                               for m in ("flow", "greedy")])
        if ordered:
            ylim, ticks = (3.62, 5.28), MultipleLocator(0.5)
        else:                                       # same scale for Ni/Mn/Co
            mid = 0.5 * (vals.min() + vals.max())
            ylim, ticks = (mid - 2.6, mid + 2.6), MultipleLocator(1)
        row_axes = []
        for j, x in enumerate(X_LEVELS):
            ax = ax_cm(fig, COL_R[0] + j * (cw + cgap), top, cw, h)
            d = ts[(system, x)]
            if ordered:     # on top: the flow series coincides with it
                ax.axhline(d["ref_tm"], zorder=5, **REF_KW)
            ax.plot(d["t"], d[f"{el}_greedy"], color=C[BASE], lw=0.6,
                    marker=MARKER[BASE], ms=3.0, mew=0, zorder=3)
            ax.plot(d["t"], d[f"{el}_flow"], color=C["flow"], lw=0.6,
                    marker=MARKER["flow"], ms=2.4, mew=0, zorder=4)
            ax.set_xlim(-3, 50.5)
            ax.set_ylim(*ylim)
            ax.xaxis.set_major_locator(MultipleLocator(20))
            ax.xaxis.set_minor_locator(MultipleLocator(10))
            ax.yaxis.set_major_locator(ticks)
            ax.tick_params(labelbottom=(i == len(TS_ROWS) - 1),
                           labelleft=(j == 0), pad=2)
            if i == 0:
                ax.set_title(rf"$x$ = {float(x):.2f}".replace("0.00", "0"),
                             fontsize=7.5, pad=3)
            row_axes.append(ax)
        first = first or row_axes[0]
        row_axes[-1].text(1.04, 0.5, f"{formula(name)}\n{el}",
                          transform=row_axes[-1].transAxes, ha="left",
                          va="center", linespacing=1.15)
        if i == len(TS_ROWS) - 1:
            row_axes[n_col // 2].set_xlabel("Time (ps)")
        top += h + rgap
    ylabel_cm(fig, COL_R[0], TOP, height, "Mean oxidation state", offset=0.85)
    label_panel(first, "c")


def method_legend(fig):
    handles = [
        Line2D([], [], label=LABEL["flow"], **SERIES["flow"]),
        Line2D([], [], label=LABEL[BASE], **SERIES[BASE]),
        Line2D([], [], label=LABEL["ref"] + " (charge counting)",
               **SERIES["ref"]),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=3,
               bbox_to_anchor=(0.5, 1.0), columnspacing=2.2,
               handlelength=2.2, borderaxespad=0.1)


def print_provenance():
    w = [max(len(r[i]) for r in PROV) for i in range(3)]
    w[1], w[2] = min(w[1], 78), min(w[2], 60)
    print("\nPROVENANCE (panel | quantity | value | source)")
    print("-" * 150)
    for p, q, v, s in PROV:
        print(f"({p}) {q:<{w[1]}}  {v:<{w[2]}}  {s}")
    print("-" * 150)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--recompute", action="store_true",
                    help="ignore _cache/fig4_*.json and rerun pyleos")
    args = ap.parse_args()

    li = load_li2o2(args.recompute)
    pol = load_polaron(args.recompute)
    load_md(args.recompute)
    ts = load_timeseries()

    fig = plt.figure(figsize=(FIG_W * style.CM, FIG_H * style.CM))
    # invisible frame: pins the tight bounding box of style.save to the canvas
    fig.add_artist(Rectangle((0, 0), 1, 1, transform=fig.transFigure,
                             fill=False, edgecolor="none", lw=0))
    method_legend(fig)
    draw_a(fig, li)
    draw_b(fig, pol)
    draw_c(fig, ts)
    style.save(fig, "fig4_delithiation")
    print_provenance()
    print(f"wrote {style.HERE / 'fig4_delithiation.pdf'} and .png "
          f"(page {PAGE_W:.1f} x {PAGE_H:.1f} cm)")


if __name__ == "__main__":
    main()
