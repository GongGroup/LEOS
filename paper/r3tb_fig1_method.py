#!/usr/bin/env python3
"""Method figure for the in-situ electronegativity assignment (r3-T-B).

Counterpart of fig1_method.py for the plugin
pyleos/methods/frozen/r3-T-B-cap-repair.py, on the clean38 benchmark
(`python3 -m pyleos bench clean --split dev`).

(a) Wannier-function (WF) centres of one H2S2 molecule, drawn by the classes
    of the method: contested WFs (bond pairs) and hybrid own WFs (lone pairs);
(b) position f of the WF centre along every terminal X-H bond of clean38
    (0 at X, 1 at H) with two dividing points: the bond midpoint (nearest-atom
    rule) and the force-balance point of the method,
    f* = sqrt(Z_eff,X) / (sqrt(Z_eff,X) + sqrt(Z_eff,H));
(c) the decision rule: own WFs of an atom (dipole-free shell, lone pairs), a
    contested WF at distances d_A and d_B, the two criteria
    chi = Z_eff / rho^2 and Z_eff / d^2, and two worked examples.

Run from this directory:  python3 r3tb_fig1_method.py
Writes r3tb_fig1_method.pdf / .png and prints every plotted number with its
source; the printout is kept in _cache/r3tb/fig1_method.run.log.

Every number comes from the plugin itself (its `_analyse`, `_Shells` and
`assign`) applied to the case as the clean evaluator shows it
(pyleos.bench.scrub).  The script only adds geometry (f, f*) and checks: the
plugin's distances and valence radii are recomputed with numpy, the dividing
point is checked against the plugin's verdict for every WF, and the oxidation
states are checked against the per-atom reference and against
results/clean38_r3-T-B.tsv.

Each panel is a function draw_<name>(ax) that takes nothing but the axes; text
and line weights scale with the axes, so a panel can be redrawn alone at
poster size (design sizes: WA x H, WB x HB, WC x H below).
"""
import contextlib
import csv
import functools
import itertools
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import style  # noqa: E402
from style import C, CM, REPO, formula  # noqa: E402

sys.path.insert(0, str(REPO))
from pyleos import bench  # noqa: E402
from pyleos.parsers import load_case  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Arc, Circle, FancyArrowPatch, Rectangle  # noqa: E402

# ---------------------------------------------------------------- Data sources
PLUGIN = REPO / "pyleos/methods/frozen/r3-T-B-cap-repair.py"
RESULTS = REPO / "results/clean38_r3-T-B.tsv"      # per-case score of the plugin
EXAMPLE = ("wave2", "h2s2")                         # panel (a); force example
CHI_EXAMPLE = ("wave2", "as4s4")                    # chi example of panel (c)
LOG = style.HERE / "_cache" / "r3tb" / "fig1_method.run.log"
NAME = "r3tb_fig1_method"

# ------------------------------------------------------- Rules (see caption)
RHO_MAX = 0.10         # A; the WF of an X-H bond lies this close to the axis
TERMINAL_RATIO = 1.25  # H is terminal if its 2nd non-H neighbour is farther
                       # than this multiple of the nearest one
THIS = C["aux1"]       # colour of the method (team brief)
MID = C["ref"]         # bond midpoint: a neutral guide line, not a method colour

# Design sizes in cm (figure 17.1 cm wide); panels scale from these.
W, H = style.DOUBLE / CM, 7.0
WA, WB, HB, WC = 3.60, 3.75, 5.10, 6.40

# Display names (labels and provenance only).
DISPLAY = {"B6H10": "B6H10", "BH6N": "NH3BH3", "CH2F2": "CH2F2", "hno3": "HNO3",
           "tmao": "(CH3)3NO", "h2s2": "H2S2", "h2s3": "H2S3", "h2s5": "H2S5",
           "srbe_oh_4": "SrBe(OH)4", "as4s4": "As4S4"}


# ------------------------------------------------------------------- Geometry
# image_vectors, min_image, axis_position and choose_view are those of
# fig1_method.py, so that panel (a) is the same projection as fig. 1(a).
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
    from the same origin) and its perpendicular offset."""
    f = float(v @ b / (b @ b))
    return f, float(np.linalg.norm(v - f * b))


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


# ------------------------------------------------------------ Method (plugin)
@functools.lru_cache(maxsize=None)
def plugin():
    return bench.load_plugin(PLUGIN)


@functools.lru_cache(maxsize=None)
def clean38():
    """(batch, name, dir, reference, site spec) of the 38 cases, resolved as
    in `python3 -m pyleos bench clean --split dev` (bench.main_clean)."""
    wanted = {(b, n): spec for b, cases in bench.CLEAN["cases"].items()
              for n, spec in cases.items()}
    held_out = {tuple(x) for x in bench.CLEAN["splits"]["test_regression"]}
    return tuple((b, n, d, ref, wanted[(b, n)])
                 for b, n, d, ref in bench.collect_clean()
                 if (b, n) in wanted and (b, n) not in held_out)


@functools.lru_cache(maxsize=None)
def analyse(case_dir):
    """The plugin's own analysis of one case.

    Returns the scrubbed case, the dict of plugin._analyse (per WF: class,
    hosts, vectors to all atoms) and, per atom, the plugin's _Shells object
    built from the own WFs exactly as _analyse builds it.  Every contested WF
    that has a heteronuclear partner gets its decision record g['note']
    (plugin 'notes': rule, Z_eff and d of winner and loser) and the two atoms
    g['top'], g['low'] the record refers to."""
    P = plugin()
    case = bench.scrub(load_case(case_dir))
    an = P._analyse(case)
    els, geo, na = an["els"], an["geo"], an["na"]
    rows = [[] for _ in range(na)]
    for i, g in enumerate(geo):
        g["index"] = i                  # position in the flattened WF list
        if g["own"]:
            rows[g["near"]].append((math.sqrt(g["d2"][g["near"]]), g["rad"],
                                    g["e"], g["pp"]))
    shells = [P._Shells(r) for r in rows]
    # The notes are appended in WF order; each records the two distances
    # sqrt(d2[top]), sqrt(d2[low]) of its WF, which identify the WF and the
    # two atoms exactly (same floats).
    notes, k = an["notes"], 0
    for g in geo:
        g["note"] = None
        if not g["contested"] or k == len(notes):
            continue
        n = notes[k]
        top = [a for a in g["hosts"] if els[a] == n["winner"]
               and math.sqrt(g["d2"][a]) == n["d"][0]]
        low = [b for b in range(na) if els[b] == n["loser"]
               and math.sqrt(g["d2"][b]) == n["d"][1]
               and top and P._in_bond_region(g["vecs"][top[0]], g["vecs"][b])]
        if top and low:
            g["note"], g["top"], g["low"] = n, top[0], low[0]
            k += 1
            # Z_eff of the record = ZVAL - inner electrons - 5/16, with the
            # inner electrons from the plugin's shells at these distances
            for a, z in zip((top[0], low[0]), n["zeff"]):
                inner, _ = shells[a].seen_from(math.sqrt(g["d2"][a]))
                assert z == an["zval"][els[a]] - inner - P.PAIR_SCREEN
    assert k == len(notes), (case_dir, k, len(notes))
    return case, an, shells


def wf_class(g):
    """Class of a WF in the method: bond (contested), lone (hybrid own WF),
    shell (dipole-free own WF)."""
    if g["contested"]:
        return "bond"
    if g["own"]:
        return "shell" if g["pp"] else "lone"
    return "other"


def ends(P, shells, an, g):
    """What the plugin compares for the contested WF g: per end (winner
    first) element, atom index, d, inner electrons, Z_eff, rho^2, Z_eff/d^2
    and chi = Z_eff/rho^2 (None without lone pairs)."""
    out = []
    for a, z, d in zip((g["top"], g["low"]), g["note"]["zeff"], g["note"]["d"]):
        inner, rho2 = shells[a].seen_from(d)
        out.append(dict(el=an["els"][a], atom=a, d=d, inner=inner, zeff=z,
                        zval=an["zval"][an["els"][a]], rho2=rho2,
                        force=z / (d * d) if d > 0 else math.inf,
                        chi=None if rho2 is None else z / rho2))
    return out


def numpy_geometry(case):
    """Own arithmetic, from POSCAR and the WF centres/spreads only: vectors
    and distances from every WF centre to every atom (minimum image) and the
    rms radius of every WF, in the order of the plugin's WF list."""
    st = case.structure
    latt = np.array(st.lattice)
    frac = np.array(st.frac)
    cen = np.vstack([np.array(ch.centres) for ch in case.channels])
    rms = np.sqrt(np.concatenate([np.array(ch.spreads) for ch in case.channels]))
    vec = min_image(frac[None, :, :] - (cen @ np.linalg.inv(latt))[:, None, :],
                    latt)                                         # WF -> atom
    return vec, np.linalg.norm(vec, axis=-1), rms


def independent_check(case, g, e):
    """Check one end e of the contested WF g against own arithmetic: the
    distance d, and rho^2 as the mean of offset^2 + spread over the lone
    pairs of that atom, i.e. the WFs whose nearest atom it is, that are not
    nearer a bond midpoint than either nucleus, lie within one rms radius of
    the nucleus, are displaced by at least PARITY rms radii and are not
    wholly inside the sphere of radius d (all WFs of these cases hold the
    same charge)."""
    P = plugin()
    vec, dist, rms = numpy_geometry(case)
    a = e["atom"]
    d = dist[g["index"], a]
    assert abs(d - e["d"]) < 1e-9, (d, e["d"])
    if e["rho2"] is not None:
        mine = []
        for i in np.flatnonzero(dist.argmin(axis=1) == a):
            off = dist[i, a]
            mid = np.linalg.norm(0.5 * (vec[i] + vec[i, a]), axis=-1)
            mid[a] = np.inf
            contested = bool(np.any((mid < off) & (mid < dist[i])))
            if (not contested and P.PARITY * rms[i] <= off <= P.REACH_SIGMA * rms[i]
                    and off + P.REACH_SIGMA * rms[i] >= d):
                mine.append(off ** 2 + rms[i] ** 2)
        assert abs(np.mean(mine) - e["rho2"]) < 1e-9, (np.mean(mine), e["rho2"])
    return True


# --------------------------------------------------------- Panel (b): survey
def xh_bonds(case, an):
    """One record per terminal X-H bond of a case.

    Selection as in fig1_method.py: X is the non-H atom nearest to the H (all
    periodic images); the H is terminal if the next non-H neighbour is
    farther than TERMINAL_RATIO x d(X-H); the WF of the bond is the centre
    nearest to that H, accepted within RHO_MAX of the X-H axis with
    0 < f < 1.  Class, verdict, Z_eff and d come from the plugin."""
    P = plugin()
    st = case.structure
    els, zval = an["els"], an["zval"]
    latt = np.array(st.lattice)
    frac = np.array(st.frac)
    _vec, wf_dist, _rms = numpy_geometry(case)
    heavy = [j for j, e in enumerate(els) if e != "H"]
    records, rejected = [], Counter()
    for ih, el in enumerate(els):
        if el != "H":
            continue
        cand = image_vectors(frac[heavy] - frac[ih], latt)      # H -> heavy
        dist = np.linalg.norm(cand, axis=-1).ravel()
        order = np.argsort(dist)
        d1, d2 = dist[order[0]], dist[order[1]]
        jx = heavy[order[0] // len(_IMG)]
        if d2 <= TERMINAL_RATIO * d1:
            rejected[f"{els[jx]}-H: H not terminal"] += 1
            continue
        for ch in sorted({g["ch"] for g in an["geo"]}):
            g = min((g for g in an["geo"] if g["ch"] == ch),
                    key=lambda g: g["d2"][ih])
            v_x, v_h = np.array(g["vecs"][jx]), np.array(g["vecs"][ih])
            bond = v_h - v_x                                    # X -> H
            assert abs(np.linalg.norm(bond) - d1) < 1e-6
            f, off = axis_position(-v_x, bond)
            if off > RHO_MAX or not 0.0 < f < 1.0:
                rejected[f"{els[jx]}-H: no WF on the bond axis"] += 1
                continue
            rec = dict(X=els[jx], f=f, off=off, L=float(d1), ratio=float(d2 / d1),
                       nearer_H=g["near"] == ih, cls=wf_class(g),
                       index=g["index"], fstar=None, margin=None, rule=None)
            hosts = g["hosts"]
            assert len(hosts) == 1 and next(iter(hosts)) in (jx, ih)
            rec["to"] = "H" if ih in hosts else "X"
            if g["contested"]:
                n = g["note"]
                assert {g["top"], g["low"]} == {jx, ih} and g["top"] in hosts
                z = dict(zip((g["top"], g["low"]), n["zeff"]))
                d = dict(zip((g["top"], g["low"]), n["d"]))
                rec.update(rule=n["rule"], zX=z[jx], zH=z[ih], dX=d[jx], dH=d[ih],
                           margin=n["force_margin"],
                           fstar=math.sqrt(z[jx])
                           / (math.sqrt(z[jx]) + math.sqrt(z[ih])))
                # the dividing point reproduces the plugin's verdict
                assert rec["rule"] == "force"
                assert (rec["f"] < rec["fstar"]) == (rec["to"] == "X")
                # own arithmetic: distances with numpy; no inner electrons,
                # so Z_eff = ZVAL - 5/16; ratio of the two forces
                assert abs(wf_dist[g["index"], jx] - rec["dX"]) < 1e-9
                assert abs(wf_dist[g["index"], ih] - rec["dH"]) < 1e-9
                assert rec["zX"] == zval[els[jx]] - P.PAIR_SCREEN
                assert rec["zH"] == zval["H"] - P.PAIR_SCREEN
                ratio = (rec["zX"] / rec["dX"] ** 2) / (rec["zH"] / rec["dH"] ** 2)
                assert abs((ratio if rec["to"] == "X" else 1.0 / ratio)
                           / rec["margin"] - 1.0) < 1e-9
            else:
                assert rec["cls"] == "lone" and rec["to"] == "H"
            records.append(rec)
    return records, dict(rejected)


@functools.lru_cache(maxsize=None)
def survey():
    """Run the plugin on the 38 cases; collect the X-H bonds, the class
    counts, the chi decisions and the checks against the reference."""
    P = plugin()
    with RESULTS.open() as fh:
        tsv = {(r["batch"], r["case"]): r for r in csv.DictReader(fh, delimiter="\t")}
    groups = defaultdict(list)             # X -> [(case, record), ...]
    checks = []
    tot = Counter()
    for batch, name, cdir, ref, spec in clean38():
        case, an, shells = analyse(cdir)
        geo, els = an["geo"], an["els"]
        tot.update(cases=1, atoms=an["na"], wfs=len(geo))
        tot.update(wf_class(g) for g in geo)
        tot.update(f"rule:{g['note']['rule']}" for g in geo if g["note"])
        tot["contested, same-element ends only"] += sum(
            g["contested"] and g["note"] is None for g in geo)
        # oxidation states against the per-atom reference and the result file
        os_ = P.assign(case)
        labels = bench.site_labels(name, case.structure, ref, spec)
        n_lab = sum(lab is not None for lab in labels)
        n_ok = sum(lab is not None and lab == o for lab, o in zip(labels, os_))
        row = tsv[(batch, name)]
        assert (int(row["atoms_ok"]), int(row["atoms_labelled"])) == (n_ok, n_lab), name
        tot.update(labelled=n_lab, correct=n_ok, full=n_ok == n_lab)
        hs = [i for i, e in enumerate(els) if e == "H"]
        if not hs:
            continue
        ne = Counter()                                  # pure nearest-atom rule
        for g in geo:
            ne[g["near"]] += an["n_pair"]
        near_os = [round(an["zval"]["H"]) - ne[i] for i in hs]
        recs, rej = xh_bonds(case, an)
        for r in recs:
            groups[r["X"]].append((name, r))
        checks.append(dict(
            name=name, dir=str(cdir.relative_to(REPO)), n_h=len(hs),
            h_ref=Counter(labels[i] for i in hs),
            h_ok=sum(labels[i] == os_[i] for i in hs),
            h_near_minus=sum(v == -1 for v in near_os),
            atoms=f"{n_ok}/{n_lab}", recs=recs, rej=rej))
    assert tot["cases"] == len(tsv) == 38
    return dict(groups=dict(groups), checks=checks, tot=tot)


def row_order(groups):
    """Bond types ordered by the dividing point f* (largest first), then by
    bond length; f* is one number per bond type here (asserted)."""
    fstar = {}
    for x, items in groups.items():
        vals = [r["fstar"] for _, r in items if r["fstar"] is not None]
        assert max(vals) - min(vals) < 1e-12
        fstar[x] = vals[0]
    order = sorted(groups, key=lambda x: (-round(fstar[x], 9),
                                          np.mean([r["L"] for _, r in groups[x]])))
    return order, fstar


# ------------------------------------------------- Panels (a), (c): examples
@functools.lru_cache(maxsize=None)
def molecule():
    """H2S2 (one molecule per cell): unwrapped atoms and WF centres, the
    class of every WF from the plugin, and the projection."""
    cdir = next(d for b, n, d, _r, _s in clean38() if (b, n) == EXAMPLE)
    case, an, shells = analyse(cdir)
    st = case.structure
    els = st.elements
    latt = np.array(st.lattice)
    frac = np.array(st.frac)
    atoms = min_image(frac - frac[0], latt)
    wfs = np.vstack([min_image(np.array(ch.centres) @ np.linalg.inv(latt)
                               - frac[0], latt) for ch in case.channels])
    origin = atoms.mean(axis=0)
    atoms, wfs = atoms - origin, wfs - origin
    kinds = [wf_class(g) for g in an["geo"]]
    # a line joins the two atoms nearest to each contested WF
    bonds = sorted({tuple(sorted(int(i) for i in np.argsort(g["d2"])[:2]))
                    for g in an["geo"] if g["contested"]})
    count, labels = defaultdict(int), []
    for e in els:                       # H1, H2, S1, S2: POSCAR order
        count[e] += 1
        labels.append(f"{e}{count[e]}")
    qa, qw, clear, short = choose_view(atoms, wfs, bonds)
    return dict(dir=cdir, case=case, an=an, shells=shells, els=els, labels=labels,
                kinds=kinds, bonds=bonds, qa=qa, qw=qw, clear=clear, short=short,
                formula="".join(f"{e}{n if n > 1 else ''}" for e, n in count.items()),
                os=plugin().assign(case))


@functools.lru_cache(maxsize=None)
def examples():
    """Worked examples of panel (c), each the least decisive WF of its kind:
    'force' - the S-H bond WF of H2S2 with the smallest Z_eff/d^2 ratio;
    'chi'   - the As-S bond WF of As4S4 with the smallest chi ratio."""
    P = plugin()
    out = {}
    for key, target in (("force", EXAMPLE), ("chi", CHI_EXAMPLE)):
        cdir = next(d for b, n, d, _r, _s in clean38() if (b, n) == target)
        case, an, shells = analyse(cdir)
        if key == "force":
            hs = {i for i, e in enumerate(an["els"]) if e == "H"}
            cand = [g for g in an["geo"] if g["note"] and g["note"]["rule"] == "force"
                    and g["near"] in hs]
            g = min(cand, key=lambda g: g["note"]["force_margin"])
        else:
            cand = [g for g in an["geo"] if g["note"] and g["note"]["rule"] == "chi"]
            g = min(cand, key=lambda g: g["note"]["chi_margin"])
        e = ends(P, shells, an, g)
        for end in e:
            independent_check(case, g, end)
        if key == "chi":
            assert abs(e[0]["chi"] / e[1]["chi"] - g["note"]["chi_margin"]) < 1e-12
        assert abs(e[0]["force"] / e[1]["force"] - g["note"]["force_margin"]) < 1e-9
        assert list(g["hosts"]) == [e[0]["atom"]]
        out[key] = dict(case=target[1], dir=str(cdir.relative_to(REPO)),
                        wf=g["index"] + 1, n_same=len(cand), ends=e,
                        rule=g["note"]["rule"],
                        margin=g["note"][f"{key}_margin"],
                        force_margin=g["note"]["force_margin"],
                        os=[plugin().assign(case)[x["atom"]] for x in e])
    return out


# -------------------------------------------------------------------- Drawing
WF_STYLE = {   # same marker for a WF class in panels (a), (b) and (c)
    "bond": dict(marker="o", mfc=C["ink"], mec=C["ink"]),
    "lone": dict(marker="o", mfc="white", mec=C["ink"]),
}
WF_LABEL = {"bond": "Bond pair (contested)", "lone": "Lone pair (own)"}


def axes_cm(ax):
    """Width and height of an axes in cm."""
    w, h = ax.figure.get_size_inches() / CM
    box = ax.get_position()
    return box.width * w, box.height * h


def cm_axes(fig, left, bottom, width, height):
    w, h = fig.get_size_inches() / CM
    return fig.add_axes([left / w, bottom / h, width / w, height / h])


def label_at(fig, ax, letter, x_cm, y_cm):
    """style.panel_label with its anchor at a given figure position (cm)."""
    w, h = fig.get_size_inches() / CM
    box = ax.get_position()
    style.panel_label(ax, letter, dx=(x_cm / w - box.x0) / box.width,
                      dy=(y_cm / h - box.y0) / box.height)


def draw_molecule(ax):
    """Panel (a): atoms, bonds and WF centres of H2S2; data units are
    Angstrom.  Design size WA x H cm."""
    mol = molecule()
    qa, qw = mol["qa"], mol["qw"]
    w_cm, h_cm = axes_cm(ax)
    k = min(w_cm / WA, h_cm / H)
    pts = np.vstack([qa, qw])
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    r_cm = 0.21 * k                                 # atom circle radius
    reserve = 1.66 * k                              # legend + scale bar
    head = 0.55 * k                                 # panel label + formula
    scale = min((w_cm - 2 * r_cm - 0.3 * k) / (hi[0] - lo[0]),
                (h_cm - reserve - head - 2 * r_cm - 0.2 * k) / (hi[1] - lo[1]))
    xc = 0.5 * (lo[0] + hi[0])
    ax.set_xlim(xc - 0.5 * w_cm / scale, xc + 0.5 * w_cm / scale)
    yc = 0.5 * (lo[1] + hi[1])
    mid = reserve + 0.5 * (h_cm - reserve - head)   # cm above the bottom
    ax.set_ylim(yc - mid / scale, yc + (h_cm - mid) / scale)
    ax.set_aspect("equal")
    ax.axis("off")
    rad = r_cm / scale
    for a, b in mol["bonds"]:
        ax.plot(*zip(qa[a], qa[b]), color=C["ink"], lw=0.7 * k, zorder=1,
                solid_capstyle="butt")
    for (x, y), lab in zip(qa, mol["labels"]):
        ax.add_patch(Circle((x, y), rad, fc="white", ec=C["ink"], lw=0.7 * k,
                            zorder=2))
        ax.text(x, y, lab, ha="center", va="center_baseline", fontsize=7 * k,
                zorder=3)
    for kind in ("bond", "lone"):
        sel = [i for i, c in enumerate(mol["kinds"]) if c == kind]
        ax.plot(qw[sel, 0], qw[sel, 1], ls="none", ms=4.2 * k, mew=0.8 * k,
                zorder=4, **WF_STYLE[kind])
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    ax.text(x1 - 0.10 * k / scale, y1 - 0.42 * k / scale, formula(mol["formula"]),
            ha="right", va="top", fontsize=8 * k)
    # 1 Angstrom scale bar, above the legend on the left
    bx, by = x0 + 0.25 * k / scale, y0 + 1.05 * k / scale
    tick = 0.07 * k / scale
    ax.plot([bx, bx + 1.0], [by, by], color=C["ink"], lw=1.0 * k,
            solid_capstyle="butt")
    for x in (bx, bx + 1.0):
        ax.plot([x, x], [by - tick, by + tick], color=C["ink"], lw=0.7 * k)
    ax.text(bx + 0.5, by + 0.12 * k / scale, "1 Å", ha="center", va="bottom",
            fontsize=7.5 * k)
    handles = [Line2D([], [], ls="none", ms=4.2 * k, mew=0.8 * k,
                      label=WF_LABEL[c], **WF_STYLE[c]) for c in ("bond", "lone")]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.05, 0.0),
              borderaxespad=0.0, borderpad=0.1, labelspacing=0.3,
              handletextpad=0.2, handlelength=1.0, fontsize=7.5 * k,
              alignment="left")
    return scale


def draw_dividing_point(ax):
    """Panel (b): strip plot of f, one row per bond type X-H, with the bond
    midpoint and the force-balance point f* of the method.  Design size of
    the axes WB x HB cm; the row labels lie outside it (keep 0.62 axes widths
    free on the left, 0.33 on the right and room for one text line above)."""
    groups = survey()["groups"]
    order, fstar = row_order(groups)
    k = axes_cm(ax)[0] / WB
    n_rows = len(order)
    tr = ax.get_yaxis_transform()            # x: axes fraction, y: data
    # right of f* the pair goes to H
    stairs = []
    for row, x in enumerate(order):
        ax.add_patch(Rectangle((fstar[x], row - 0.5), 1.0 - fstar[x], 1.0,
                               fc=THIS, ec="none", alpha=0.22, zorder=0))
        stairs += [(fstar[x], row - 0.5), (fstar[x], row + 0.5)]
    ax.plot(*zip(*stairs), color=THIS, lw=1.2 * k, zorder=2,
            solid_joinstyle="miter", solid_capstyle="butt")
    ax.axvline(0.5, color=MID, lw=1.0 * k, ls=(0, (4, 2)), zorder=2)
    table = []
    for row, x in enumerate(order):
        by_case = defaultdict(list)
        for name, r in groups[x]:
            by_case[name].append(r)
        lanes = (np.linspace(-0.22, 0.22, len(by_case))
                 if len(by_case) > 1 else [0.0])
        for off, (name, rs) in zip(lanes, by_case.items()):
            rs = sorted(rs, key=lambda r: r["f"])
            fs = np.array([r["f"] for r in rs])
            cls = {r["cls"] for r in rs}
            assert len(cls) == 1
            cls = cls.pop()
            # filled dots may merge; open ones are stepped so each ring shows
            half = 0.055 if cls == "bond" else min(0.14, 0.02 * (len(fs) - 1))
            spread = (np.linspace(-half, half, len(fs)) if len(fs) > 1
                      else np.zeros(1))
            ax.plot(fs, row + off + spread, ls="none", zorder=3,
                    ms=(3.0 if cls == "bond" else 2.9) * k,
                    mew=(0.0 if cls == "bond" else 0.5) * k, **WF_STYLE[cls])
            table.append(dict(bond=f"{x}-H", case=name, n=len(fs), cls=cls,
                              f_min=fs.min(), f_max=fs.max(),
                              nearer_H=sum(r["nearer_H"] for r in rs),
                              to=Counter(r["to"] for r in rs),
                              zX=rs[0].get("zX"), zH=rs[0].get("zH"),
                              fstar=rs[0]["fstar"],
                              margin=(min(r["margin"] for r in rs)
                                      if cls == "bond" else None)))
        z_x = next(r["zX"] for _, r in groups[x] if r["fstar"] is not None)
        ax.text(-0.455, row, f"{x}–H", transform=tr, ha="left", va="center",
                fontsize=7.5 * k)
        ax.text(-0.04, row, f"({z_x:.2f})", transform=tr, ha="right",
                va="center", fontsize=7.5 * k)
        ax.text(1.035, row, f"$n$ = {len(groups[x])}", transform=tr,
                ha="left", va="center", fontsize=7.5 * k)
    head = dict(transform=ax.transAxes, va="bottom", fontsize=7.5 * k)
    ax.text(0.0, 1.02, "X", ha="center", **head)
    ax.text(1.0, 1.02, "H", ha="center", **head)
    ax.text(0.5, 1.02, "midpoint", ha="center", color=MID, **head)
    ax.text(fstar[order[0]], 1.02, "$f\\,$*", ha="center", color=THIS, **head)
    note = dict(ha="center", va="center", fontsize=7 * k, linespacing=1.25,
                zorder=4)
    y_mid = 0.5 * (n_rows - 1) - 0.5
    ax.text(0.235, y_mid, "pair to X\nH(+1)", **note)
    ax.text(0.5 * (1.0 + fstar[order[0]]), y_mid, "pair\nto H\nH(−1)", **note)
    ax.set_xlim(0, 1)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"], fontsize=7 * k)
    ax.set_yticks([])
    ax.tick_params(axis="x", length=2.5 * k, width=0.6 * k, pad=3.5 * k)
    for sp in ax.spines.values():
        sp.set_linewidth(0.6 * k)
    ax.set_xlabel("Position of WF centre along X–H, $f$", fontsize=8 * k,
                  labelpad=4 * k)
    ax.text(-0.56, 0.5, "Bond X–H  ($Z_{\\mathrm{eff}}$ of X)",
            transform=ax.transAxes, rotation=90, ha="center", va="center",
            fontsize=8 * k)
    return order, fstar, table


def draw_rule(ax):
    """Panel (c): the decision rule; data units are centimetres of the
    design size WC x H."""
    ex = examples()
    w_cm, h_cm = axes_cm(ax)
    k = min(w_cm / WC, h_cm / H)
    ax.set_xlim(0, WC)
    ax.set_ylim(0, H)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    ink, muted = C["ink"], C["muted"]

    # ---- schematic: atoms A and B, own WFs, one contested WF
    yb, xa, xb = 5.38, 1.42, 4.92
    xw = xa + 0.60 * (xb - xa)
    r_atom, r_shell = 0.21, 0.36
    rho_a, rho_b = 0.88, 0.76
    ax.plot([xa, xb], [yb, yb], color=ink, lw=0.7 * k, zorder=1)
    ax.add_patch(Circle((xa, yb), r_shell, fc=C["light"], ec=muted,
                        lw=0.5 * k, zorder=2))
    for x, lab in ((xa, "A"), (xb, "B")):
        ax.add_patch(Circle((x, yb), r_atom, fc="white", ec=ink, lw=0.7 * k,
                            zorder=3))
        ax.text(x, yb, lab, ha="center", va="center_baseline", fontsize=7 * k,
                zorder=4)
    lone = []
    for x, sign, off, rho, lab in ((xa, -1, 0.60, rho_a, "A"),
                                   (xb, +1, 0.54, rho_b, "B")):
        for ang in (62.0, -62.0):
            t = math.radians(ang if sign > 0 else 180.0 - ang)
            lone.append((x + off * math.cos(t), yb + off * math.sin(t)))
        mid = 0.0 if sign > 0 else 180.0
        ax.add_patch(Arc((x, yb), 2 * rho, 2 * rho, theta1=mid - 42,
                         theta2=mid + 42, color=muted, lw=0.6 * k,
                         ls=(0, (3, 2)), zorder=1))
        start = r_shell if lab == "A" else r_atom
        ax.add_patch(FancyArrowPatch(
            (x + sign * start, yb), (x + sign * rho, yb), arrowstyle="-|>",
            mutation_scale=4.5 * k, lw=0.6 * k, color=ink, shrinkA=0,
            shrinkB=0, zorder=2))
        ax.text(x + sign * 0.5 * (start + rho) - sign * 0.02, yb + 0.07,
                f"$\\rho_{{\\mathrm{{{lab}}}}}$", ha="center", va="bottom",
                fontsize=7.5 * k)
    ax.plot(*zip(*lone), ls="none", ms=4.2 * k, mew=0.8 * k, zorder=5,
            **WF_STYLE["lone"])
    ax.plot([xw], [yb], ls="none", ms=4.2 * k, mew=0.8 * k, zorder=5,
            **WF_STYLE["bond"])
    # distances d_A, d_B
    yd, tick = yb - 0.60, 0.07
    for x0, x1, lab in ((xa, xw, "A"), (xw, xb, "B")):
        ax.plot([x0, x1], [yd, yd], color=ink, lw=0.6 * k, zorder=2)
        ax.text(0.5 * (x0 + x1), yd - 0.09, f"$d_{{\\mathrm{{{lab}}}}}$",
                ha="center", va="top", fontsize=7.5 * k)
    for x, top in ((xa, yb - r_shell), (xw, yb - 0.10), (xb, yb - r_atom)):
        ax.plot([x, x], [yd - tick, yd + tick], color=ink, lw=0.6 * k, zorder=2)
        ax.plot([x, x], [yd + tick, top - 0.03], color=muted, lw=0.4 * k,
                ls=(0, (1.5, 1.5)), zorder=1)
    # labels of the three WF classes
    lab = dict(fontsize=7 * k, va="bottom")
    ax.text(xw, yb + 0.16, "bond pair", ha="center", **lab)
    ax.text(lone[2][0] + 0.10, lone[2][1] + 0.15, "lone pair", ha="center", **lab)
    xs, ys = xa + 0.62, yb + 0.80
    ax.text(xs + 0.05, ys, "dipole-free shell", ha="left", va="center",
            fontsize=7 * k)
    ax.plot([xs, xa + r_shell * math.cos(math.radians(52)) + 0.01],
            [ys - 0.02, yb + r_shell * math.sin(math.radians(52)) + 0.01],
            color=muted, lw=0.5 * k, zorder=2)

    # ---- the rule
    xm = 0.5 * WC
    ax.text(xm, 4.06, "$Z_{\\mathrm{eff}} = Z_{\\mathrm{val}} - "
            "n_{\\mathrm{inner}} - 5/16$", ha="center", va="center",
            fontsize=8 * k)
    ax.text(xm, 3.68, "$n_{\\mathrm{inner}}$: electrons of own WFs inside "
            "radius $d$", ha="center", va="center", fontsize=7 * k)
    ax.text(xm, 3.20, "The pair goes to the end with the larger",
            ha="center", va="center", fontsize=7.5 * k)
    gap = 0.24
    wcol = 0.5 * (WC - gap) - 0.06
    blocks = (
        (0.06, "chi", "$\\chi = Z_{\\mathrm{eff}}\\,/\\,\\rho^{2}$",
         "lone pairs on both ends",
         ("$Z_{\\mathrm{eff}}$", "$\\rho^{2}$ (Å$^{2}$)", "$\\chi$"),
         lambda e: (e["zeff"], e["rho2"], e["chi"]), ("{:.2f}", "{:.2f}", "{:.2f}")),
        (0.06 + wcol + gap, "force", "$Z_{\\mathrm{eff}}\\,/\\,d^{\\,2}$",
         "otherwise",
         ("$Z_{\\mathrm{eff}}$", "$d$ (Å)", "$Z_{\\mathrm{eff}}/d^{\\,2}$"),
         lambda e: (e["zeff"], e["d"], e["force"]), ("{:.2f}", "{:.3f}", "{:.2f}")),
    )
    y_f, y_c = 2.76, 2.38
    y_top, y_head, y_sep, y_bot = 2.10, 1.87, 1.64, 0.80
    y_rows = (1.43, 1.02)
    base = 0.095                 # baseline below the centre of a 7.5 pt line
    for x0, key, crit, cond, heads, vals, fmts in blocks:
        e = ex[key]
        xc = x0 + 0.5 * wcol
        ax.text(xc, y_f - 0.11, crit, ha="center", va="baseline",
                fontsize=8.5 * k)
        ax.text(xc, y_c, cond, ha="center", va="center", fontsize=7.5 * k)
        cols = (x0 + 0.06, x0 + 0.30 * wcol, x0 + 0.56 * wcol, x0 + 0.85 * wcol)
        for y, lw in ((y_top, 0.7), (y_sep, 0.4), (y_bot, 0.7)):
            ax.plot([x0, x0 + wcol], [y, y], color=ink, lw=lw * k,
                    solid_capstyle="butt", zorder=2)
        for x, h in zip(cols[1:], heads):
            ax.text(x, y_head - base, h, ha="center", va="baseline",
                    fontsize=7 * k)
        for row, (y, end) in enumerate(zip(y_rows, e["ends"])):
            if row == 0:                             # the end that wins
                ax.add_patch(Rectangle((x0, y - 0.19), wcol, 0.38, fc=THIS,
                                       ec="none", alpha=0.30, zorder=0))
            ax.text(cols[0], y - base, end["el"], ha="left", va="baseline",
                    fontsize=7.5 * k)
            for j, (x, v, fmt) in enumerate(zip(cols[1:], vals(end), fmts)):
                ax.text(x, y - base, fmt.format(v), ha="center", va="baseline",
                        fontsize=7.5 * k,
                        fontweight="bold" if row == 0 and j == 2 else "normal")
        pair = "–".join(sorted((x["el"] for x in e["ends"]),
                               key=lambda el: (el == "H", el)))
        ax.text(xc, y_bot - 0.30,
                f"{pair} pair in {formula(DISPLAY[e['case']])} "
                f"→ {e['ends'][0]['el']}",
                ha="center", va="center", fontsize=7.5 * k)
    return ex


# ----------------------------------------------------------------------- Main
def build():
    fig = plt.figure(figsize=(W * CM, H * CM))
    y_lab = H - 0.42                    # baseline of the panel labels
    xb = WA + 0.10                      # left edge of the block of panel (b)
    ax_a = cm_axes(fig, 0.0, 0.0, WA, H)
    ax_b = cm_axes(fig, xb + 2.20, 0.95, WB, HB)
    ax_c = cm_axes(fig, W - WC, 0.0, WC, H)
    draw_molecule(ax_a)
    order, fstar, table = draw_dividing_point(ax_b)
    ex = draw_rule(ax_c)
    label_at(fig, ax_a, "a", 0.50, y_lab)
    label_at(fig, ax_b, "b", xb + 0.50, y_lab)
    label_at(fig, ax_c, "c", W - WC + 0.50, y_lab)
    return fig, order, fstar, table, ex


def report(order, fstar, table, ex):
    P = plugin()
    sv, mol = survey(), molecule()
    tot = sv["tot"]
    rel = lambda p: str(Path(p).relative_to(REPO))
    print(f"Method: {rel(PLUGIN)}")
    print(f"  constants: PAIR_SCREEN = {P.PAIR_SCREEN} (5/16), REACH_SIGMA = "
          f"{P.REACH_SIGMA}, PARITY = {P.PARITY}, FLOAT_EPS = {P.FLOAT_EPS}")
    print(f"Benchmark: clean38 = `python3 -m pyleos bench clean --split dev`; "
          f"{tot['cases']} cases, {tot['atoms']} atoms, {tot['wfs']} WFs")
    print(f"  WF classes (plugin _analyse): own {tot['lone'] + tot['shell']} "
          f"(hybrid {tot['lone']}, dipole-free {tot['shell']}), contested "
          f"{tot['bond']}, other {tot['other']}")
    print(f"  contested WFs with ends of two elements: "
          f"{tot['rule:force'] + tot['rule:chi']} (decided by Z_eff/d^2: "
          f"{tot['rule:force']}, by chi: {tot['rule:chi']}); same-element ends "
          f"only: {tot['contested, same-element ends only']}")
    print(f"  atoms correct {tot['correct']}/{tot['labelled']}, cases fully "
          f"correct {tot['full']}/{tot['cases']} - recomputed here and equal, "
          f"case by case, to {rel(RESULTS)}")

    print("\nPanel (a): Wannier centres of one molecule, classes of the method")
    print(f"  source: {rel(mol['dir'])}/POSCAR, EOSCAR (#WF block); classes "
          f"from plugin _analyse")
    print(f"  view: clearance {mol['clear']:.2f} x required, most foreshortened "
          f"bond drawn at {mol['short']:.2f} of its length")
    for i, g in enumerate(mol["an"]["geo"]):
        a = g["near"]
        off, rad = math.sqrt(g["d2"][a]), g["rad"]
        hosts = ", ".join(f"{mol['labels'][h]} {s:.3f}"
                          for h, s in sorted(g["hosts"].items()))
        extra = (f"offset/rms = {off / rad:.3f} (dipole-free below "
                 f"{P.PARITY})" if g["own"] else
                 "nearer a bond midpoint than either nucleus")
        print(f"  WF {i + 1}: {wf_class(g):5s} nearest {mol['labels'][a]} "
              f"{off:.3f} A, rms radius {rad:.3f} A, {extra}; hosts: {hosts}")
    print("  oxidation states (plugin assign): " + ", ".join(
        f"{lab} {o:+d}" for lab, o in zip(mol["labels"], mol["os"])))

    print("\nPanel (b): terminal X-H bonds of clean38, f = 0 at X and 1 at H")
    print("  cases (all clean38 cases that contain H):")
    for c in sv["checks"]:
        parts = []
        by_x = defaultdict(list)
        for r in c["recs"]:
            by_x[r["X"]].append(r)
        for x, rs in by_x.items():
            fs = [r["f"] for r in rs]
            parts.append(f"{x}-H n={len(rs)} f={min(fs):.3f}..{max(fs):.3f}")
        parts += [f"[excluded {n}: {why}]" for why, n in c["rej"].items()]
        ref = ", ".join(f"{v:+d} x{n}" for v, n in sorted(c["h_ref"].items()))
        print(f"  {DISPLAY[c['name']]:10s} {c['dir']}: " + "; ".join(parts))
        print(f"      H atoms {c['n_h']}: reference {ref}; plugin correct "
              f"{c['h_ok']}/{c['n_h']}; nearest-atom rule gives -1 to "
              f"{c['h_near_minus']}; all atoms {c['atoms']} correct")
    n_h = sum(c["n_h"] for c in sv["checks"])
    ref_minus = sum(c["h_ref"][-1] for c in sv["checks"])
    print(f"  H atoms in clean38: {n_h}; reference -1: {ref_minus}, +1: "
          f"{sum(c['h_ref'][1] for c in sv['checks'])}; plugin correct: "
          f"{sum(c['h_ok'] for c in sv['checks'])}/{n_h}; pure nearest-atom "
          f"rule gives -1 to {sum(c['h_near_minus'] for c in sv['checks'])}")
    recs = [r for x in order for _, r in sv["groups"][x]]
    print(f"  accepted bonds: {len(recs)}; d(X-H) = "
          f"{min(r['L'] for r in recs):.3f}..{max(r['L'] for r in recs):.3f} A; "
          f"largest offset of a WF centre from its bond axis "
          f"{max(r['off'] for r in recs):.3f} A (limit {RHO_MAX}); smallest "
          f"2nd/1st neighbour distance ratio {min(r['ratio'] for r in recs):.2f}"
          f" (limit {TERMINAL_RATIO})")
    print("  plotted rows (top to bottom; lanes in the order listed); Z_eff and"
          " verdict from the plugin, f* = sqrt(Z_X)/(sqrt(Z_X)+sqrt(Z_H)):")
    print(f"  {'bond':5s} {'case':10s} {'n':>3s} {'class':5s} {'f_min':>6s} "
          f"{'f_max':>6s} {'Z_X':>6s} {'Z_H':>6s} {'f*':>6s} {'nearer H':>8s} "
          f"{'to X':>4s} {'to H':>4s} {'min ratio':>9s}")
    for t in table:
        z = (f"{t['zX']:6.4f} {t['zH']:6.4f} {t['fstar']:6.4f}"
             if t["cls"] == "bond" else f"{'-':>6s} {'-':>6s} {'-':>6s}")
        m = f"{t['margin']:9.3f}" if t["margin"] is not None else f"{'-':>9s}"
        print(f"  {t['bond']:5s} {DISPLAY[t['case']]:10s} {t['n']:3d} "
              f"{t['cls']:5s} {t['f_min']:6.3f} {t['f_max']:6.3f} {z} "
              f"{t['nearer_H']:8d} {t['to']['X']:4d} {t['to']['H']:4d} {m}")
    for x in order:
        rs = [r for _, r in sv["groups"][x]]
        print(f"  {x}-H: n = {len(rs)}, f* = {fstar[x]:.4f}, nearer H than X "
              f"{sum(r['nearer_H'] for r in rs)}, pair to X "
              f"{sum(r['to'] == 'X' for r in rs)}, pair to H "
              f"{sum(r['to'] == 'H' for r in rs)}")
    n_cont = sum(r["cls"] == "bond" for r in recs)
    print(f"  all rows: {len(recs)} X-H bond WFs; nearer H than X: "
          f"{sum(r['nearer_H'] for r in recs)} (a nearest-atom rule gives these "
          f"to H); in-situ rule: pair to X {sum(r['to'] == 'X' for r in recs)}, "
          f"pair to H {sum(r['to'] == 'H' for r in recs)}")
    print("  Z_eff = ZVAL - 5/16 at both ends of every contested X-H WF (no "
          "inner electrons); distances and force ratios reproduced with numpy")
    print(f"  {n_cont} of them are contested and decided by Z_eff/d^2 (f < f* "
          f"<=> pair to X, checked for each); {len(recs) - n_cont} are "
          f"one-centre WFs of H (open circles; centre nearer H than the bond "
          f"midpoint, f > 0.75) and stay on H")

    print("\nPanel (c): decision rule and worked examples (plugin numbers)")
    for key in ("chi", "force"):
        e = ex[key]
        print(f"  {key}: {e['dir']}, WF {e['wf']} (least decisive of "
              f"{e['n_same']} WFs of this kind in the case); rule used by the "
              f"plugin: {e['rule']}")
        for end, o in zip(e["ends"], e["os"]):
            rho = ("-" if end["rho2"] is None else
                   f"{end['rho2']:.4f} A^2 (rho = {math.sqrt(end['rho2']):.4f} A)")
            chi = "-" if end["chi"] is None else f"{end['chi']:.4f}"
            print(f"    {end['el']}{end['atom'] + 1}: ZVAL {end['zval']:g}, inner "
                  f"{end['inner']:.4f}, Z_eff {end['zeff']:.4f}, d "
                  f"{end['d']:.4f} A, Z_eff/d^2 {end['force']:.4f}, rho^2 {rho}, "
                  f"chi {chi}; OS {o:+d}")
        print(f"    ratio winner/loser: {key} {e['margin']:.4f}"
              + (f" (Z_eff/d^2 {e['force_margin']:.4f})" if key == "chi" else "")
              + f"; pair to {e['ends'][0]['el']}; distances and rho^2 "
              f"reproduced with numpy from POSCAR and WF centres/spreads")


class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, s):
        for st in self.streams:
            st.write(s)

    def flush(self):
        for st in self.streams:
            st.flush()


def main():
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("w") as fh, contextlib.redirect_stdout(Tee(sys.stdout, fh)):
        fig, order, fstar, table, ex = build()
        style.save(fig, NAME)
        report(order, fstar, table, ex)
        print(f"\nwrote {style.HERE / (NAME + '.pdf')}")
        print(f"wrote {style.HERE / (NAME + '.png')}")
        print(f"log   {LOG}")


if __name__ == "__main__":
    main()
