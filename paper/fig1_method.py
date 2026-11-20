#!/usr/bin/env python3
"""Figure 1 - Method.

(a) Wannier-function (WF) centres of one H2S2 molecule (orthographic
    projection), classed geometrically as bond pairs or lone pairs;
(b) fractional position f of every X-H bond-pair WF centre along its bond
    (0 at X, 1 at H) in the scored benchmark cases that contain hydrogen,
    grouped by the heavy atom X and ordered by the Allen electronegativity
    difference chi(X) - chi(H); f = 0.5 is where a nearest-atom rule hands
    the pair to H;
(c) the assignment of panel (a) written as a flow network: source -> WFs ->
    atoms -> sink, with the edges that carry flow in the optimal solution.

Run from this directory:  python3 fig1_method.py
Needs numpy + matplotlib + stdlib; pyleos is imported from the repository
root.  Every plotted number is computed from files under cases/ through
pyleos (load_case, capacity, assign_case); a provenance table is printed.
Nothing is typed in by hand except display names of compounds (labels only)
and the geometric thresholds below, which are stated in the caption.
"""
import itertools
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import style  # noqa: E402
from style import C, CM, REPO, formula  # noqa: E402

sys.path.insert(0, str(REPO))
from pyleos.bench import collect  # noqa: E402
from pyleos.capacity import capacity_bounds, eneg  # noqa: E402
from pyleos.parsers import load_case  # noqa: E402
from pyleos.transport import assign_case  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle  # noqa: E402

# ---------------------------------------------------------------- Data sources
EXAMPLE = REPO / "cases/h2s2/v3"                   # panels (a) and (c)
CLEAN_JSON = REPO / "pyleos/sets/clean_set.json"    # flags, printed only

# ------------------------------------------------------- Rules (see caption)
BETA = 10.0            # weight of the electronegativity term (default)
RHO_MAX = 0.10         # A; a bond-pair centre lies this close to the axis
TERMINAL_RATIO = 1.25  # H is terminal if its 2nd non-H neighbour is farther
                       # than this multiple of the nearest one
TIE_TOL = 0.01         # A^2; edges this close in cost count as degenerate
SHADE_DISAGREEMENT = True   # panel (b): grey where nearest atom != IUPAC

# Display names (labels and provenance only).
NAME = {
    "B6H10": "B6H10", "BH6N": "NH3BH3", "CH2F2": "CH2F2",
    "AuClB_C6H5_3_P_C6H5_3_2": "AuCl[B(C6H5)3][P(C6H5)3]2",
    "hno3": "HNO3", "tmao": "(CH3)3NO", "h2s2": "H2S2", "h2s3": "H2S3",
    "h2s5": "H2S5", "h3po3": "H3PO3", "hi3": "HI.I2 ('hi3')",
    "rhcl_co__p_c6h5_3_2": "RhCl(CO)(SO2)[P(C6H5)3]2",
    "rhh": "RhH", "srbe_oh_4": "SrBe(OH)4",
}

# ------------------------------------------------------------------- Geometry
_IMG = np.array(list(itertools.product(range(-2, 3), repeat=3)), dtype=float)


def image_vectors(dfrac, lattice):
    """Cartesian vectors of all periodic images (|n| <= 2) of fractional
    differences; shape (..., 125, 3)."""
    df = np.asarray(dfrac, dtype=float)
    df = df - np.round(df)
    return (df[..., None, :] + _IMG) @ lattice


def min_image(dfrac, lattice):
    """Shortest Cartesian image of each fractional difference."""
    cand = image_vectors(dfrac, lattice)
    idx = np.argmin((cand ** 2).sum(-1), axis=-1)
    return np.take_along_axis(cand, idx[..., None, None], axis=-2)[..., 0, :]


def axis_position(v, b):
    """Fractional position f of point v along the vector b (both measured
    from the same origin) and its perpendicular offset rho."""
    f = float(v @ b / (b @ b))
    return f, float(np.linalg.norm(v - f * b))


# --------------------------------------------------------- Panel (b): survey
def benchmark_cases():
    """(batch, name, directory, has_reference) of the 74-case benchmark,
    resolved exactly as in `pyleos bench full` (pyleos.bench.collect)."""
    for batch, name, cdir, ref in collect():
        yield batch, name, cdir, bool(ref)


def xh_bond_wfs(case):
    """All X-H bond-pair WFs of one case.

    For every H atom: X is its nearest non-H atom (all periodic images
    count); the H is terminal if the next non-H neighbour is farther than
    TERMINAL_RATIO times d(X-H).  The bond-pair WF is the WF centre nearest
    to that H, accepted if it lies within RHO_MAX of the X-H axis with
    0 < f < 1.  Returns (records, rejected); each record is a dict."""
    st = case.structure
    els = st.elements
    latt = np.array(st.lattice)
    linv = np.linalg.inv(latt)
    frac = np.array(st.frac)
    heavy = [j for j, e in enumerate(els) if e != "H"]
    records, rejected = [], defaultdict(int)
    for ih, el in enumerate(els):
        if el != "H":
            continue
        cand = image_vectors(frac[heavy] - frac[ih], latt)   # H -> heavy
        dist = np.linalg.norm(cand, axis=-1).ravel()
        order = np.argsort(dist)
        d1, d2 = dist[order[0]], dist[order[1]]
        jx = heavy[order[0] // len(_IMG)]
        if d2 <= TERMINAL_RATIO * d1:
            rejected[f"{els[jx]}-H: H not terminal"] += 1
            continue
        bond = -cand.reshape(-1, 3)[order[0]]                  # X -> H
        for ch in case.channels:
            wfrac = np.array(ch.centres) @ linv
            h_to_w = min_image(wfrac - frac[ih], latt)
            k = int(np.argmin(np.linalg.norm(h_to_w, axis=1)))
            x_to_w = h_to_w[k] + bond
            f, rho = axis_position(x_to_w, bond)
            if rho > RHO_MAX or not 0.0 < f < 1.0:
                rejected[f"{els[jx]}-H: no WF on the bond axis"] += 1
                continue
            records.append(dict(X=els[jx], f=f, rho=rho, d=float(d1),
                                ratio=float(d2 / d1),
                                nearer_H=bool(np.linalg.norm(h_to_w[k])
                                              < np.linalg.norm(x_to_w))))
    return records, dict(rejected)


def hydrogen_os(case):
    """Oxidation state of every H atom under (i) a pure nearest-atom rule
    (each WF to its nearest atom), (ii) the EOSCAR in-situ baseline column
    and (iii) the min-cost flow (beta = BETA).  Printed, not plotted."""
    res = assign_case(case, beta=BETA)
    els = case.structure.elements
    ne = [0] * len(els)
    for w in res.wfs:
        ne[w.nearest] += w.supply
    zval = case.eoscar.zval
    near = [round(zval[e]) - n for e, n in zip(els, ne)]
    hs = [i for i, e in enumerate(els) if e == "H"]
    return ([near[i] for i in hs], [case.eoscar.greedy_os[i] for i in hs],
            [res.os_int[i] for i in hs])


def survey():
    """Collect X-H bond-pair WFs over the scored H-containing cases."""
    flags = {(e["batch"], e["case"]): ",".join(e["reason"])
             for e in json.loads(CLEAN_JSON.read_text())["excluded"]}
    groups = defaultdict(list)          # X -> [(case name, f), ...]
    rows, accepted = [], []
    hyd_tot = np.zeros(4, dtype=int)    # H atoms; H(-1) by three rules
    for batch, name, cdir, has_ref in benchmark_cases():
        if not (cdir / "EOSCAR").is_file():
            continue
        case = load_case(cdir)
        if "H" not in case.structure.elements:
            continue
        rel = cdir.relative_to(REPO)
        flag = flags.get((batch, name), "-")
        if not has_ref:
            rows.append((name, str(rel), "skipped: no reference", flag))
            continue
        recs, rej = xh_bond_wfs(case)
        accepted += recs
        near, base, flow = hydrogen_os(case)
        hyd_tot += [len(near)] + [sum(v == -1 for v in vals)
                                  for vals in (near, base, flow)]
        by_x = defaultdict(list)
        for r in recs:
            by_x[r["X"]].append(r)
            groups[r["X"]].append((name, r["f"]))
        parts = []
        for x, rs in by_x.items():
            fs = [r["f"] for r in rs]
            parts.append(f"{x}-H n={len(rs)} f={min(fs):.3f}..{max(fs):.3f} "
                         f"(nearer H: {sum(r['nearer_H'] for r in rs)})")
        for reason, n in rej.items():
            parts.append(f"[excluded {n}: {reason}]")
        hyd = (f"H(-1): nearest-atom {sum(v == -1 for v in near)}/{len(near)}"
               f", EOSCAR baseline {sum(v == -1 for v in base)}/{len(base)}"
               f", flow {sum(v == -1 for v in flow)}/{len(flow)}")
        rows.append((name, str(rel), "; ".join(parts) + "; " + hyd, flag))
    margins = dict(rho_max=max(r["rho"] for r in accepted),
                   ratio_min=min(r["ratio"] for r in accepted),
                   d_min=min(r["d"] for r in accepted),
                   d_max=max(r["d"] for r in accepted),
                   hyd=[int(v) for v in hyd_tot])
    return groups, rows, margins


# ------------------------------------------------- Panels (a), (c): molecule
def load_molecule(case_dir):
    """One molecule (cell holds a single molecule): unwrapped atom and WF
    coordinates, WF classes, bonds, and the optimal assignment."""
    case = load_case(case_dir)
    st = case.structure
    els = st.elements
    latt = np.array(st.lattice)
    frac = np.array(st.frac)
    atoms = min_image(frac - frac[0], latt)
    wfs = np.vstack([min_image(np.array(ch.centres) @ np.linalg.inv(latt)
                               - frac[0], latt) for ch in case.channels])
    origin = atoms.mean(axis=0)
    atoms, wfs = atoms - origin, wfs - origin

    # Class of each WF from its two nearest atoms A, B: bond pair if the
    # centre lies within RHO_MAX of the A-B axis and between the nuclei.
    kinds, bonds, geom = [], set(), []
    for w in wfs:
        dist = np.linalg.norm(atoms - w, axis=1)
        a, b = (int(i) for i in np.argsort(dist)[:2])
        f, rho = axis_position(w - atoms[a], atoms[b] - atoms[a])
        is_bond = rho <= RHO_MAX and 0.0 < f < 1.0
        kinds.append("bond" if is_bond else "lone")
        geom.append((a, b, dist[a], dist[b], f, rho))
        if is_bond:
            bonds.add((min(a, b), max(a, b)))

    count = defaultdict(int)
    labels = []
    for e in els:                       # H1, H2, S1, S2: POSCAR order
        count[e] += 1
        labels.append(f"{e}{count[e]}")

    res = assign_case(case, beta=BETA)
    chi_max = max(eneg(e) for e in set(els))
    d2 = ((wfs[:, None, :] - atoms[None, :, :]) ** 2).sum(-1)
    cost = d2 + BETA * (1.0 - np.array([eneg(e) for e in els]) / chi_max)
    bounds = [capacity_bounds(e, case.eoscar.zval[e]) for e in els]
    sent = np.zeros_like(cost)
    for i, w in enumerate(res.wfs):
        for a, n in w.dests:
            sent[i, a] += n
    # The cost recomputed here from the published formula must reproduce
    # the solver's own objective (valid when no n_min reward is active).
    if all(lo == 0 for lo, _ in bounds):
        assert abs((sent * cost).sum() - res.transport_cost) < 1e-6
    # Degenerate alternatives: an unused edge whose cost equals that of the
    # used one within TIE_TOL and whose use violates no capacity.
    ties = []
    for i in range(len(wfs)):
        for a in np.flatnonzero(sent[i] > 0):
            for a2 in range(len(els)):
                if (a2 != a and sent[i, a2] == 0
                        and abs(cost[i, a2] - cost[i, a]) < TIE_TOL
                        and res.ne[a2] + sent[i, a] <= bounds[a2][1]
                        and res.ne[a] - sent[i, a] >= bounds[a][0]):
                    ties.append((i, int(a2), float(cost[i, a2] - cost[i, a])))
    formula_str = "".join(f"{e}{n if n > 1 else ''}"
                          for e, n in count.items())
    return dict(case=case, els=els, labels=labels, atoms=atoms, wfs=wfs,
                formula=formula_str,
                kinds=kinds, bonds=sorted(bonds), geom=geom, res=res,
                cost=cost, sent=sent, bounds=bounds, ties=ties,
                supply=[w.supply for w in res.wfs],
                nelect=case.eoscar.nelect_wann)


def choose_view(atoms, wfs, bonds, n_dir=20000):
    """Orthographic view: among the directions that keep every pair of drawn
    objects apart (atom-atom 0.70 A, atom-WF 0.30 A, WF-WF 0.25 A, each with
    a 20 % margin), take the one whose most foreshortened bond is longest.
    Returns 2D coordinates (atoms, WFs), principal axis vertical."""
    pts = np.vstack([atoms, wfs])
    na, n = len(atoms), len(pts)
    req = np.full((n, n), 0.25)
    req[:na, :] = req[:, :na] = 0.30
    req[:na, :na] = 0.70
    iu = np.triu_indices(n, 1)
    blen = [np.linalg.norm(atoms[a] - atoms[b]) for a, b in bonds]
    best = None
    for k in range(n_dir):
        z = 1.0 - (k + 0.5) / n_dir                 # Fibonacci hemisphere
        r, phi = np.sqrt(1.0 - z * z), np.pi * (1.0 + 5.0 ** 0.5) * k
        nv = np.array([r * np.cos(phi), r * np.sin(phi), z])
        t = np.array([1.0, 0, 0]) if abs(nv[0]) < 0.9 else np.array([0, 1.0, 0])
        u = np.cross(nv, t)
        u /= np.linalg.norm(u)
        q = pts @ np.stack([u, np.cross(nv, u)]).T
        clear = (np.linalg.norm(q[iu[0]] - q[iu[1]], axis=1) / req[iu]).min()
        short = min(np.linalg.norm(q[a] - q[b]) / l
                    for (a, b), l in zip(bonds, blen)) if bonds else 1.0
        key = (clear >= 1.2, short if clear >= 1.2 else clear)
        if best is None or key > best[0]:
            best = (key, q, clear, short)
    _, q, clear, short = best
    qa = q[:na] - q[:na].mean(axis=0)
    _, vecs = np.linalg.eigh(qa.T @ qa)
    ey = vecs[:, -1]                                 # principal axis -> y
    ex = np.array([ey[1], -ey[0]])
    q = (q - q[:na].mean(axis=0)) @ np.stack([ex, ey]).T
    if q[0, 1] < 0:                                  # first atom at the top
        q[:, 1] *= -1
    if q[0, 0] > 0:                                  # ... and on the left
        q[:, 0] *= -1
    return q[:na], q[na:], clear, short


# -------------------------------------------------------------------- Drawing
WF_STYLE = {   # same marker for a WF class in panels (a), (b) and (c)
    "bond": dict(marker="o", mfc=C["ink"], mec=C["ink"]),
    "lone": dict(marker="o", mfc="white", mec=C["ink"]),
}
WF_LABEL = {"bond": "Bond pair", "lone": "Lone pair"}


def cm_axes(fig, left, bottom, width, height):
    w, h = fig.get_size_inches() / CM
    return fig.add_axes([left / w, bottom / h, width / w, height / h])


def label_at(fig, ax, letter, x_cm, y_cm):
    """style.panel_label with its anchor at a given figure position (cm)."""
    w, h = fig.get_size_inches() / CM
    box = ax.get_position()
    style.panel_label(ax, letter, dx=(x_cm / w - box.x0) / box.width,
                      dy=(y_cm / h - box.y0) / box.height)


def panel_a(ax, mol, qa, qw, w_cm, h_cm):
    """Atoms, bonds and WF centres; data units are Angstrom."""
    pts = np.vstack([qa, qw])
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    r_cm = 0.21                                     # atom circle radius
    reserve = 1.30                                  # legend + scale bar, cm
    head = 0.75                                     # panel label + formula
    scale = min((w_cm - 2 * r_cm - 0.3) / (hi[0] - lo[0]),
                (h_cm - reserve - head - 2 * r_cm - 0.2) / (hi[1] - lo[1]))
    xc = 0.5 * (lo[0] + hi[0])
    ax.set_xlim(xc - 0.5 * w_cm / scale, xc + 0.5 * w_cm / scale)
    # centre the molecule between the header and the legend
    yc = 0.5 * (lo[1] + hi[1])
    mid = reserve + 0.5 * (h_cm - reserve - head)   # cm above the bottom
    ax.set_ylim(yc - mid / scale, yc + (h_cm - mid) / scale)
    ax.set_aspect("equal")
    ax.axis("off")
    rad = r_cm / scale
    for a, b in mol["bonds"]:
        ax.plot(*zip(qa[a], qa[b]), color=C["ink"], lw=0.7, zorder=1,
                solid_capstyle="butt")
    for (x, y), lab in zip(qa, mol["labels"]):
        ax.add_patch(Circle((x, y), rad, fc="white", ec=C["ink"], lw=0.7,
                            zorder=2))
        ax.text(x, y, lab, ha="center", va="center_baseline", fontsize=7,
                zorder=3)
    for kind in ("bond", "lone"):
        sel = [i for i, k in enumerate(mol["kinds"]) if k == kind]
        ax.plot(qw[sel, 0], qw[sel, 1], ls="none", ms=4.2, mew=0.8,
                zorder=4, **WF_STYLE[kind])
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    ax.text(x1 - 0.10 / scale, y1 - 0.42 / scale, formula(mol["formula"]),
            ha="right", va="top", fontsize=8)
    # 1 Angstrom scale bar, bottom left
    bx, by = x0 + 0.25 / scale, y0 + 0.22 / scale
    tick = 0.07 / scale
    ax.plot([bx, bx + 1.0], [by, by], color=C["ink"], lw=1.0,
            solid_capstyle="butt")
    for x in (bx, bx + 1.0):
        ax.plot([x, x], [by - tick, by + tick], color=C["ink"], lw=0.7)
    ax.text(bx + 0.5, by + 0.12 / scale, "1 Å", ha="center", va="bottom",
            fontsize=7.5)
    handles = [Line2D([], [], ls="none", ms=4.2, mew=0.8, label=WF_LABEL[k],
                      **WF_STYLE[k]) for k in ("bond", "lone")]
    ax.legend(handles=handles, loc="lower right", bbox_to_anchor=(1.0, 0.0),
              borderaxespad=0.0, borderpad=0.1, labelspacing=0.3,
              handletextpad=0.2, handlelength=1.0, title="WF centre",
              title_fontsize=7.5, alignment="left")
    return scale


def panel_b(ax, groups):
    """Strip plot of f, one row per bond type X-H."""
    order = sorted(groups, key=lambda x: -(eneg(x) - eneg("H")))
    dchi = [eneg(x) - eneg("H") for x in order]
    n_rows = len(order)
    tr = ax.get_yaxis_transform()            # x: axes fraction, y: data
    # Where a nearest-atom rule contradicts the electronegativity rule:
    # f > 0.5 for chi(X) > chi(H), f < 0.5 for chi(X) < chi(H).
    n_pos = sum(d > 0 for d in dchi)
    if SHADE_DISAGREEMENT:
        for x0, y0, height in ((0.5, -0.5, n_pos),
                               (0.0, n_pos - 0.5, n_rows - n_pos)):
            ax.add_patch(Rectangle((x0, y0), 0.5, height, fc=C["light"],
                                   ec="none", alpha=0.55, zorder=0))
    if 0 < n_pos < n_rows:
        ax.axhline(n_pos - 0.5, color=C["ink"], lw=0.5, zorder=1)
    ax.axvline(0.5, color=C["nearest"], lw=1.0, ls=(0, (4, 2)), zorder=2)
    table = []
    for row, x in enumerate(order):
        by_case = defaultdict(list)
        for name, f in groups[x]:
            by_case[name].append(f)
        lanes = (np.linspace(-0.22, 0.22, len(by_case))
                 if len(by_case) > 1 else [0.0])
        for off, (name, fs) in zip(lanes, by_case.items()):
            fs = np.sort(fs)
            spread = (np.linspace(-0.055, 0.055, len(fs)) if len(fs) > 1
                      else np.zeros(1))
            ax.plot(fs, row + off + spread, ls="none", ms=3.0, mew=0.0,
                    marker="o", mfc=C["ink"], zorder=3)
            table.append((f"{x}-H", f"{dchi[row]:+.3f}", name, len(fs),
                          fs.min(), np.median(fs), fs.max(),
                          int((fs > 0.5).sum())))
        ax.text(-0.455, row, f"{x}–H", transform=tr, ha="left",
                va="center", fontsize=7.5)
        ax.text(-0.04, row,
                f"({dchi[row]:+.2f})".replace("-", "−"),
                transform=tr, ha="right", va="center", fontsize=7.5)
        ax.text(1.035, row, f"$n$ = {len(groups[x])}", transform=tr,
                ha="left", va="center", fontsize=7.5)
    head = dict(transform=ax.transAxes, ha="center", va="bottom",
                fontsize=7.5)
    ax.text(0.0, 1.02, "X", **head)
    ax.text(1.0, 1.02, "H", **head)
    ax.set_xlim(0, 1)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"])
    ax.set_yticks([])
    ax.set_xlabel("Position of WF centre along X–H, $f$")
    ax.text(-0.56, 0.5, "Bond X–H  "
            r"($\chi_{\mathrm{X}} - \chi_{\mathrm{H}}$)",
            transform=ax.transAxes, rotation=90, ha="center", va="center",
            fontsize=8)
    return order, table


def panel_c(ax, mol, qa, qw, w_cm, h_cm):
    """Flow network of the molecule; data units are centimetres."""
    ax.set_xlim(0, w_cm)
    ax.set_ylim(0, h_cm)
    ax.axis("off")
    pt = 1.0 / (72 * CM)                            # cm per point
    xs, xw, xa, xt = 0.28, 1.55, 3.80, w_cm - 0.28
    y_top, y_bot = h_cm - 1.12, 2.05
    r_node, r_atom = 0.17, 0.21
    r_wf = 0.5 * 4.2 * pt + 0.02
    # Top to bottom here = top to bottom in panel (a).  WFs are ordered by
    # the mean rank of their two nearest atoms (bond pair) or the rank of
    # their host (lone pair), which keeps edge crossings low.
    a_order = [int(a) for a in np.argsort(-qa[:, 1], kind="stable")]
    rank = {a: r for r, a in enumerate(a_order)}
    key = []
    for i, (a, b, *_rest) in enumerate(mol["geom"]):
        r = (0.5 * (rank[a] + rank[b]) if mol["kinds"][i] == "bond"
             else float(rank[a]))
        key.append((r, -qw[i, 1], i))
    w_order = [i for _, _, i in sorted(key)]
    yw = dict(zip(w_order, np.linspace(y_top, y_bot, len(w_order))))
    inset = 0.5 * (y_top - y_bot) / len(w_order)
    ya = dict(zip(a_order, np.linspace(y_top - inset, y_bot + inset,
                                       len(a_order))))
    ym = 0.5 * (y_top + y_bot)
    sent, ne = mol["sent"], mol["res"].ne

    def edge(p, q, used, r_p, r_q, ls="-"):
        if used:
            ax.add_patch(FancyArrowPatch(
                p, q, arrowstyle="-|>", mutation_scale=5.0, lw=1.0, ls=ls,
                color=C["flow"], shrinkA=r_p / pt, shrinkB=r_q / pt,
                zorder=3, joinstyle="miter", capstyle="butt"))
        else:
            ax.plot(*zip(p, q), color=C["light"], lw=0.6, zorder=1)

    tied = {(i, a) for i, a, _ in mol["ties"]}
    for i in w_order:
        edge((xs, ym), (xw, yw[i]), True, r_node, r_wf)
        for a in a_order:
            if (i, a) in tied:
                edge((xw, yw[i]), (xa, ya[a]), True, r_wf, r_atom,
                     ls=(0, (2.5, 1.5)))
            else:
                edge((xw, yw[i]), (xa, ya[a]), sent[i, a] > 0, r_wf, r_atom)
    for a in a_order:
        edge((xa, ya[a]), (xt, ym), ne[a] > 0, r_atom, r_node)
        lo, hi = mol["bounds"][a]
        ax.text(xa + r_atom + 0.06, ya[a] + (0.21 if ya[a] >= ym else -0.21),
                f"[{lo}, {hi}]", ha="left", va="center", fontsize=7)
    for x, lab in ((xs, "s"), (xt, "t")):
        ax.add_patch(Circle((x, ym), r_node, fc="white", ec=C["ink"], lw=0.7,
                            zorder=4))
        ax.text(x, ym, f"${lab}$", ha="center", va="center_baseline",
                fontsize=7.5, zorder=5)
    for i in w_order:
        ax.plot([xw], [yw[i]], ls="none", ms=4.2, mew=0.8, zorder=5,
                **WF_STYLE[mol["kinds"][i]])
    for a in a_order:
        ax.add_patch(Circle((xa, ya[a]), r_atom, fc="white", ec=C["ink"],
                            lw=0.7, zorder=4))
        ax.text(xa, ya[a], mol["labels"][a], ha="center",
                va="center_baseline", fontsize=7, zorder=5)
    head = dict(ha="center", va="bottom", fontsize=7.5)
    ax.text(xw, y_top + 0.27, "WF $i$", **head)
    ax.text(xa, y_top + 0.27, "atom $j$", **head)
    supply = sorted(set(mol["supply"]))
    note = dict(ha="center", va="top", fontsize=7.5, linespacing=1.3)
    y_note = y_bot - 0.36
    ax.text(0.5 * (xs + xw), y_note,
            "supply\n" + ", ".join(f"{s}$e$" for s in supply), **note)
    ax.text(0.5 * (xw + xa), y_note, "cost\n$c_{ij}$", **note)
    ax.text(0.5 * (xa + xt) - 0.05, y_note,
            "capacity\n" r"[$n_j^{\min}$, $n_j^{\max}$]", **note)
    ax.text(0.5 * w_cm, 0.05,
            r"$c_{ij} = d_{ij}^{\,2} + \beta\,(1 - \chi_j/\chi_{\max})$"
            r"$,\quad \Sigma_j\, n_j \leq N_{\mathrm{elect}}$",
            ha="center", va="bottom", fontsize=7.5)
    return w_order, a_order


# ----------------------------------------------------------------------- Main
def main():
    groups, rows, margins = survey()
    mol = load_molecule(EXAMPLE)
    qa, qw, clear, short = choose_view(mol["atoms"], mol["wfs"], mol["bonds"])

    w_cm, h_cm = style.DOUBLE / CM, 7.0
    fig = plt.figure(figsize=(w_cm * CM, h_cm * CM))
    y_lab = h_cm - 0.42                 # baseline of the panel labels
    wa, wb, wc = 3.75, 3.70, 5.75       # widths: (a), axes of (b), (c)
    xb = wa + 0.35                      # left edge of the block of panel (b)
    ax_a = cm_axes(fig, 0.0, 0.0, wa, h_cm)
    ax_b = cm_axes(fig, xb + 2.30, 0.95, wb, h_cm - 0.95 - 0.95)
    ax_c = cm_axes(fig, w_cm - wc, 0.0, wc, h_cm)

    panel_a(ax_a, mol, qa, qw, wa, h_cm)
    order, table = panel_b(ax_b, groups)
    w_order, a_order = panel_c(ax_c, mol, qa, qw, wc, h_cm)
    label_at(fig, ax_a, "a", 0.50, y_lab)
    label_at(fig, ax_b, "b", xb + 0.50, y_lab)
    label_at(fig, ax_c, "c", w_cm - wc + 0.50, y_lab)
    style.save(fig, "fig1_method")

    # ------------------------------------------------------------ provenance
    rel = EXAMPLE.relative_to(REPO)
    print("Panel (a): Wannier centres of one molecule")
    print(f"  source: {rel}/POSCAR, {rel}/EOSCAR (#WF block)")
    print(f"  view: clearance {clear:.2f} x required, most foreshortened "
          f"bond drawn at {short:.2f} of its length")
    for i, (a, b, da, db, f, rho) in enumerate(mol["geom"]):
        print(f"  WF {i + 1}: {WF_LABEL[mol['kinds'][i]].lower():9s} nearest "
              f"{mol['labels'][a]} {da:.3f} A, {mol['labels'][b]} {db:.3f} A,"
              f" f = {f:+.3f}, offset from axis {rho:.3f} A")
    print("\nPanel (b): X-H bond-pair WFs, f = 0 at X and 1 at H")
    print("  cases surveyed (scored benchmark cases that contain H):")
    for name, path, text, flag in rows:
        print(f"  {NAME.get(name, name):28s} {path}")
        print(f"      {text}")
        if flag != "-":
            print(f"      clean_set.json flag: {flag}")
    print(f"  accepted bonds: d(X-H) = {margins['d_min']:.3f}.."
          f"{margins['d_max']:.3f} A; largest offset of a WF centre from "
          f"its bond axis {margins['rho_max']:.3f} A (limit {RHO_MAX}); "
          f"smallest 2nd/1st neighbour distance ratio "
          f"{margins['ratio_min']:.2f} (limit {TERMINAL_RATIO})")
    print("  H atoms of these cases assigned -1: nearest-atom rule "
          "{1} of {0}, EOSCAR baseline {2}, min-cost flow {3}"
          .format(*margins["hyd"]))
    print("  plotted groups (top to bottom; lanes in the order listed):")
    print(f"  {'bond':5s} {'dchi':>7s} {'case':22s} {'n':>3s} {'f_min':>6s} "
          f"{'f_med':>6s} {'f_max':>6s} {'f>0.5':>5s}")
    for bond, d, name, n, lo, med, hi, nh in table:
        print(f"  {bond:5s} {d:>7s} {name:22s} {n:3d} {lo:6.3f} {med:6.3f} "
              f"{hi:6.3f} {nh:5d}")
    for x in order:
        fs = np.array([f for _, f in groups[x]])
        print(f"  {x}-H total n = {len(fs)}, median f = {np.median(fs):.3f}, "
              f"f > 0.5 in {int((fs > 0.5).sum())}")
    pos = np.array([f for x in order if eneg(x) > eneg("H")
                    for _, f in groups[x]])
    neg = np.array([f for x in order if eneg(x) < eneg("H")
                    for _, f in groups[x]])
    print(f"  chi(X) > chi(H): {int((pos > 0.5).sum())} of {len(pos)} bond "
          f"pairs lie nearer to H;  chi(X) < chi(H): "
          f"{int((neg > 0.5).sum())} of {len(neg)}")
    print("\nPanel (c): flow network of the molecule in panel (a)")
    print(f"  beta = {BETA:g}, Allen electronegativity (pyleos.capacity), "
          f"N_elect = {mol['nelect']:g} (EOSCAR header), "
          f"total cost {mol['res'].transport_cost:.4f}")
    print("  atoms top to bottom: " + ", ".join(
        f"{mol['labels'][a]} [{mol['bounds'][a][0]}, {mol['bounds'][a][1]}]"
        f" n = {mol['res'].ne[a]} OS = {mol['res'].os_int[a]:+d}"
        .replace("+0", "0")
        for a in a_order))
    d2 = ((mol["wfs"][:, None, :] - mol["atoms"][None, :, :]) ** 2).sum(-1)
    for i in w_order:
        used = ", ".join(f"{mol['labels'][a]} ({int(mol['sent'][i, a])} e, "
                         f"c = {mol['cost'][i, a]:.3f})"
                         for a in np.flatnonzero(mol["sent"][i] > 0))
        near = int(np.argmin(d2[i]))
        alt = "".join(f"; equal-cost alternative {mol['labels'][a]} "
                      f"(dc = {dc:+.4f})"
                      for j, a, dc in mol["ties"] if j == i)
        print(f"  WF {i + 1} ({WF_LABEL[mol['kinds'][i]].lower()}): -> {used}"
              f"; nearest atom {mol['labels'][near]} "
              f"(d = {d2[i, near] ** 0.5:.3f} A, "
              f"c = {mol['cost'][i, near]:.3f}){alt}")
    print(f"\nwrote {style.HERE / 'fig1_method.pdf'}")
    print(f"wrote {style.HERE / 'fig1_method.png'}")


if __name__ == "__main__":
    main()
