#!/usr/bin/env python3
"""Figure 4 - When can the assignment be trusted.

(a) reliability index R_elem of the 72 scored benchmark cases, grouped by
    agreement with the reference oxidation states;
(b) number of agreeing cases as a function of the electronegativity weight
    beta;
(c) excess total Wannier spread of five transition-metal oxides (two spin
    channels each): archived wannier90 run (100 iterations), gpuw90 run from
    the same starting gauge, gpuw90 runs from 16 random starting gauges.

Run from this directory:  python3 figS2_reliability.py
Needs numpy + matplotlib + stdlib.  Every plotted number is read from the
files listed under "Data sources"; a provenance table is printed.

Panel (c) uses PROVISIONAL data that are not in the repository (held-out
test inputs).  Only Wannier-quality numbers are read from there (total
spread, convergence flags, stray counts) - never an oxidation state.  If the
directory is missing the panel is skipped with a notice.

No ghost panel: results/ghost_ab_summary.tsv does not support one cleanly;
ghost_check() prints the reason at every run.
"""
import csv
import json
import re
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True      # leave no __pycache__ in the repository
sys.path.insert(0, str(Path(__file__).resolve().parent))
import style  # noqa: E402
from style import C, CM, MARKER, REPO, formula  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.legend_handler import HandlerTuple  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from matplotlib.ticker import FixedLocator, MultipleLocator  # noqa: E402

# ---------------------------------------------------------------- Data sources
REL_TSV = REPO / "results/bench74_golden.tsv"
BETA_TSV = REPO / "results/scan_beta.tsv"
GHOST_TSV = REPO / "results/ghost_ab_summary.tsv"
REF_JSON = REPO / "pyleos/sets/reference_os.json"    # ghost_check() only
# Provisional, outside the repository; the only place this path appears.
HOLDOUT_A1 = Path("/private/tmp/claude-501/-Users-mastreina-Desktop-EOS/"
                  "f1ef2bf1-0cdf-4500-9d1a-d3a47241909d/scratchpad/holdout/a1")

# Quoted values - used for the printed cross-check only, never plotted.
# docs/RESULTS.md sections 2, 7 and 8; the stray totals are those stated in
# reports/multistart.json and review/verdict.json of the hold-out directory.
DOC = {"n_ok": 64, "n_not": 8, "med_ok": 83.0, "med_not": 25.55,
       "beta10_scan": 62, "stray_arch": 13, "stray_std": 2, "stray_best": 0}

# Display names (labels only).  "(legacy)": entry of the historical data set
# (run v1), which shares its formula with a newer entry of the benchmark.
NAME = {
    ("v1_mol", "As4S4"): "As4S4 (legacy)",
    ("v1_mol", "CSe2"): "CSe2 (legacy)",
    ("v1_mol", "CH2F2"): "CH2F2",
    ("wave1", "ca3aun"): "Ca3AuN",
    ("wave1", "cse2"): "CSe2",
    ("wave1", "tmao"): "(CH3)3NO",
    ("wave2", "hi3"): "HI3",
    ("wave2", "rhcl_co__p_c6h5_3_2"): "RhCl(CO)(PPh3)2",
    ("wave2", "wau12"): "WAu12",
    ("wave2", "wcl5"): "WCl5",
    "coo": "CoO", "feo": "FeO", "mno": "MnO", "nio": "NiO", "reo3": "ReO3",
}
# Cases whose reference follows another convention than the electronegativity
# rule (docs/RESULTS.md, sections 3 and 8); they get an asterisk.
CONVENTION = {("wave1", "ca3aun"), ("wave2", "wau12")}

ANG2 = "\u00c5\u00b2"   # Angstrom squared
GEQ, DAGGER = "\u2265", "\u2020"
PROV = []          # rows of the provenance table


def prov(panel, what, value, source):
    PROV.append((panel, what, str(value), source))


def rel(path):
    try:
        return str(Path(path).relative_to(REPO))
    except ValueError:
        return str(path)


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ------------------------------------------------------------------- geometry
def axes_cm(fig, left, bottom, width, height):
    w, h = fig.get_size_inches() / CM
    return fig.add_axes([left / w, bottom / h, width / w, height / h])


def axes_points(ax):
    """Width and height of an axes in points."""
    fig = ax.figure
    box = ax.get_position()
    w, h = fig.get_size_inches() * 72
    return box.width * w, box.height * h


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


# -------------------------------------------------------------------- panel a
def panel_a(ax, rows):
    scored = [r for r in rows if r["status"] != "run"]
    groups = [[r for r in scored if r["status"] == "OK"],
              [r for r in scored if r["status"] != "OK"]]
    ticklabels = ["Agrees with\nreference", "Deviates from\nreference"]
    med = [float(np.median([float(r["R_elem"]) for r in g])) for g in groups]

    ax.set_xlim(-0.55, 2.32)
    ax.set_ylim(-4, 104)
    w_pt, h_pt = axes_points(ax)
    x_per_pt = (ax.get_xlim()[1] - ax.get_xlim()[0]) / w_pt
    pt_per_y = h_pt / (ax.get_ylim()[1] - ax.get_ylim()[0])

    ms = 3.3
    for gi, grp in enumerate(groups):
        r_el = np.array([float(r["R_elem"]) for r in grp])
        stray = np.array([int(r["n_stray"]) > 0 for r in grp])
        q1, q2, q3 = np.percentile(r_el, [25, 50, 75])
        ax.add_patch(Rectangle((gi - 0.27, q1), 0.54, q3 - q1, fc=C["light"],
                               ec="none", alpha=0.75, zorder=1))
        ax.plot([gi - 0.31, gi + 0.31], [q2, q2], color=C["ink"], lw=1.0,
                solid_capstyle="butt", zorder=2)
        x = gi + swarm(r_el * pt_per_y, dist=ms + 0.5) * x_per_pt
        ax.plot(x[~stray], r_el[~stray], ls="none", marker=MARKER["flow"],
                ms=ms, mfc=C["flow"], mec="white", mew=0.3, zorder=4)
        ax.plot(x[stray], r_el[stray], ls="none", marker=MARKER["aux3"],
                ms=ms + 1.0, mfc=C["flow"], mec="white", mew=0.3, zorder=4)
        for r, xi in zip(grp, x):
            r["_x"] = xi

    # names of the deviating cases, with leader lines
    bad = groups[1]
    y_pts = np.array([float(r["R_elem"]) for r in bad])
    sep = 8.6 / pt_per_y
    y_lab = spread_labels(y_pts, sep, lo=0.5, hi=99.5)
    x_lab = 1.36
    for r, yl in zip(bad, y_lab):
        key = (r["batch"], r["case"])
        txt = formula(NAME.get(key, r["case"])) + ("*" if key in CONVENTION else "")
        x0, y0, x1 = r["_x"], float(r["R_elem"]), x_lab - 0.035
        dx, dy = (x1 - x0) / x_per_pt, (yl - y0) * pt_per_y     # points
        gap = 2.6 / np.hypot(dx, dy)
        ax.plot([x0 + gap * (x1 - x0), x1], [y0 + gap * (yl - y0), yl],
                color=C["muted"], lw=0.4, zorder=3, solid_capstyle="round")
        ax.text(x_lab, yl, txt, ha="left", va="center", fontsize=7)
    # agreeing cases that lie below the median of the deviating group
    low_ok = [r for r in groups[0] if float(r["R_elem"]) < med[1]]
    for r in low_ok:
        ax.annotate(formula(NAME.get((r["batch"], r["case"]), r["case"])),
                    (r["_x"], float(r["R_elem"])), xytext=(4.5, 0),
                    textcoords="offset points", ha="left", va="center",
                    fontsize=7)

    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"{t}\n($n$ = {len(g)})"
                        for t, g in zip(ticklabels, groups)])
    ax.tick_params(axis="x", which="both", length=0, pad=4)
    ax.yaxis.set_major_locator(MultipleLocator(20))
    ax.yaxis.set_minor_locator(MultipleLocator(10))
    ax.set_ylabel(r"Reliability index $R_{\mathrm{elem}}$")
    handles = [
        Line2D([], [], ls="none", marker=MARKER["flow"], ms=ms, mfc=C["flow"],
               mec="white", mew=0.3),
        Line2D([], [], ls="none", marker=MARKER["aux3"], ms=ms + 1.0,
               mfc=C["flow"], mec="white", mew=0.3)]
    ax.legend(handles, ["no stray WFs", GEQ + " 1 stray WF"], loc="upper right",
              handletextpad=0.2, borderaxespad=0.5, labelspacing=0.35)

    n_stray = [sum(int(r["n_stray"]) > 0 for r in g) for g in groups]
    prov("a", "scored cases (status != run)", len(scored), rel(REL_TSV))
    prov("a", "agreeing (status OK) / deviating (~, FAIL)",
         f"{len(groups[0])} / {len(groups[1])}", rel(REL_TSV))
    prov("a", "median R_elem, agreeing / deviating",
         f"{med[0]:.2f} / {med[1]:.2f}", rel(REL_TSV))
    prov("a", "cases with n_stray > 0, agreeing / deviating",
         f"{n_stray[0]} / {n_stray[1]}", rel(REL_TSV))
    for r in sorted(bad, key=lambda r: float(r["R_elem"])):
        prov("a", f"deviating: {r['batch']}/{r['case']}",
             f"R_elem {r['R_elem']}, dev {r['dev']}, n_stray {r['n_stray']}",
             rel(REL_TSV))
    for r in low_ok:
        prov("a", f"agreeing, R_elem < {med[1]:.2f}: {r['batch']}/{r['case']}",
             f"R_elem {r['R_elem']}", rel(REL_TSV))
    return {"n_ok": len(groups[0]), "n_not": len(groups[1]),
            "med_ok": med[0], "med_not": med[1]}


# -------------------------------------------------------------------- panel b
def panel_b(ax, rows):
    beta = np.array([float(r["beta"]) for r in rows])
    n_ok = np.array([int(r["n_OK"]) for r in rows])
    n_all = {int(r["n_OK"]) + int(r["n_mid"]) + int(r["n_FAIL"]) for r in rows}
    assert len(n_all) == 1, "beta scan: number of scored cases varies"
    n_all = n_all.pop()
    # the reference (default) beta is encoded in the column name
    col = [k for k in rows[0] if k.startswith("changed_vs_beta")][0]
    beta0 = float(re.search(r"beta([\d.]+)$", col).group(1))

    ax.axvline(beta0, color=C["muted"], lw=0.5, ls=(0, (4, 3)), zorder=1)
    ax.plot(beta, n_ok, color=C["flow"], lw=1.0, marker=MARKER["flow"], ms=3.8,
            mfc=C["flow"], mec="white", mew=0.4, zorder=3)
    lo, hi = n_ok.min(), n_ok.max()
    ax.set_xlim(1.5, 21.5)
    ax.set_ylim(lo - 1.6, hi + 1.6)
    ax.text(beta0 + 0.45, lo - 1.25, "default", ha="left", va="bottom",
            fontsize=7)
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.set_xlabel(r"Electronegativity weight $\beta$ (" + ANG2 + ")")
    ax.set_ylabel(f"Agreeing cases (of {n_all})")

    at0 = int(n_ok[np.argmin(np.abs(beta - beta0))])
    prov("b", "beta values", " ".join(f"{b:g}" for b in beta), rel(BETA_TSV))
    prov("b", "n_OK", " ".join(str(n) for n in n_ok), rel(BETA_TSV))
    prov("b", "scored cases per beta", n_all, rel(BETA_TSV))
    prov("b", "default beta (column name) / n_OK there", f"{beta0:g} / {at0}",
         rel(BETA_TSV))
    for r in rows:
        if r[col] != "-":
            prov("b", f"changes at beta = {r['beta']}", r[col], rel(BETA_TSV))
    try:                                   # cross-check against the code default
        import inspect
        sys.path.insert(0, str(REPO))
        from pyleos.transport import assign_case
        prov("b", "default beta in pyleos.transport.assign_case",
             inspect.signature(assign_case).parameters["beta"].default,
             "pyleos/transport.py")
    except Exception as exc:               # pragma: no cover
        prov("b", "default beta in pyleos", f"not checked ({exc})", "-")
    return {"beta10_scan": at0}


# -------------------------------------------------------------------- panel c
def load_holdout():
    """Wannier-quality numbers of the held-out oxide runs (no oxidation
    states are read).  Returns None when the directory is absent."""
    if not HOLDOUT_A1.is_dir():
        return None
    multi = json.loads((HOLDOUT_A1 / "multistart_summary.json").read_text())
    std = json.loads((HOLDOUT_A1 / "compare_standard.json").read_text())
    std = {(r["case"], r["spin"]): r for r in std["rows"]}
    flat = {(r["case"], int(r["spin"])): r
            for r in read_tsv(HOLDOUT_A1 / "multistart_summary.tsv")}
    out = []
    for r in multi["rows"]:
        key = (r["case"], r["spin"])
        s, t = std[key], flat[key]
        ch = {
            "case": r["case"], "spin": r["spin"],
            "arch": r["archived_omega"], "std": r["standard_omega"],
            "std_conv": bool(s["wan_converged"]),
            "seed": np.array([x["omega"] for x in r["seeds"]]),
            "seed_conv": np.array([bool(x["converged"]) for x in r["seeds"]]),
            "seed_stray": np.array([x["n_stray"] for x in r["seeds"]]),
            "stray_arch": r["archived_n_stray"],
            "stray_std": r["standard_n_stray"],
            "stray_best": r["best_n_stray"],
            "best_conv": bool(r["best_converged"]),
        }
        # the three files must tell the same story
        assert abs(ch["arch"] - s["omega_arch"]["total"]) < 1e-8, key
        assert abs(ch["std"] - s["omega_gpu"]["total"]) < 1e-8, key
        assert abs(ch["arch"] - float(t["archived_omega"])) < 1e-8, key
        assert abs(ch["std"] - float(t["standard_omega"])) < 1e-8, key
        assert abs(ch["seed"].min() - float(t["best_omega"])) < 1e-8, key
        assert ch["stray_arch"] == len(s["stray_arch"]) == int(t["archived_n_stray"])
        assert ch["stray_std"] == len(s["stray_gpu"]) == int(t["standard_n_stray"])
        assert int(ch["seed_conv"].sum()) == int(t["n_converged"]), key
        out.append(ch)
    return out


def panel_c(ax, chans):
    lin = 1e-3                                   # linear below, log above
    ax.set_yscale("symlog", linthresh=lin, linscale=0.9)
    cases = []
    for ch in chans:
        if ch["case"] not in cases:
            cases.append(ch["case"])
    ms_a, ms_s, ms_r = 3.8, 4.6, 3.0
    col_a, col_s, col_r = C["ref"], C["aux2"], C["aux1"]
    xticks, xlabels, ymax = [], [], 0.0
    for i, ch in enumerate(chans):
        ref = min(ch["arch"], ch["std"], ch["seed"].min())
        d_arch, d_std = ch["arch"] - ref, ch["std"] - ref
        d_seed = ch["seed"] - ref
        ymax = max(ymax, d_arch, d_std, d_seed.max())
        if not ch["best_conv"]:                  # no start converged
            ax.axvspan(i - 0.5, i + 0.5, fc=C["light"], ec="none", alpha=0.5,
                       zorder=0)
        # archived wannier90 run: stopped after 100 iterations (open symbol)
        xa, xs = i - 0.36, i - 0.23
        ax.plot(xa, d_arch, ls="none", marker=MARKER["ref"], ms=ms_a,
                mfc="white", mec=col_a, mew=0.9, zorder=4)
        ax.plot(xs, d_std, ls="none", marker=MARKER["aux2"], ms=ms_s,
                mfc=col_s if ch["std_conv"] else "white", mec=col_s, mew=0.7,
                zorder=4)
        order = np.argsort(d_seed, kind="stable")
        xr = i - 0.10 + 0.51 * np.arange(len(order)) / (len(order) - 1)
        for conv in (True, False):
            sel = ch["seed_conv"][order] == conv
            ax.plot(xr[sel], d_seed[order][sel], ls="none",
                    marker=MARKER["aux1"], ms=ms_r,
                    mfc=col_r if conv else "none", mec=col_r, mew=0.6,
                    zorder=2 if conv else 3)
        for x, y, n in ((xa, d_arch, ch["stray_arch"]),
                        (xs, d_std, ch["stray_std"])):
            if n:
                ax.annotate(str(n), (x, y), xytext=(0, 3.4),
                            textcoords="offset points", ha="center",
                            va="bottom", fontsize=7)
        xticks.append(i)
        arrow = r"$\uparrow$" if ch["spin"] == 1 else r"$\downarrow$"
        xlabels.append(arrow if ch["best_conv"] else arrow + DAGGER)

        tag = f"{ch['case']} spin {ch['spin']}"
        prov("c", f"{tag}: Omega archived / same start / lowest (A^2)",
             f"{ch['arch']:.4f} / {ch['std']:.4f} / {ref:.4f}",
             "multistart_summary.json")
        prov("c", f"{tag}: excess of 16 random starts, min / median / max",
             f"{d_seed.min():.1e} / {np.median(d_seed):.1e} / {d_seed.max():.1e}",
             "multistart_summary.json")
        prov("c", f"{tag}: converged, same start / random starts",
             f"{ch['std_conv']} / {int(ch['seed_conv'].sum())} of {len(d_seed)}",
             "compare_standard.json, multistart_summary.json")
        prov("c", f"{tag}: stray WFs archived / same start / best random",
             f"{ch['stray_arch']} / {ch['stray_std']} / {ch['stray_best']}",
             "multistart_summary.json")

    ax.axhline(lin, color=C["muted"], lw=0.4, ls=(0, (1, 2.5)), zorder=1)
    for k in range(1, len(cases)):               # separators between compounds
        ax.axvline(2 * k - 0.5, color=C["muted"], lw=0.4, zorder=1)
    ax.set_xlim(-0.5, len(chans) - 0.5)
    ax.set_ylim(-2.6e-4, 6 * ymax)
    decades = [10.0 ** k for k in range(-3, int(np.floor(np.log10(6 * ymax))) + 1)]
    ax.yaxis.set_major_locator(FixedLocator([0.0] + decades))
    ax.set_yticklabels(["0"] + [f"$10^{{{int(np.log10(d))}}}$" for d in decades])
    ax.yaxis.set_minor_locator(FixedLocator(
        [m * d for d in decades for m in range(2, 10)] +
        [0.25 * lin, 0.5 * lin, 0.75 * lin]))
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=8)
    ax.tick_params(axis="x", which="both", length=0, pad=3)
    for k, name in enumerate(cases):
        ax.text(2 * k + 0.5, -0.135, formula(NAME.get(name, name)),
                transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=8)
    ax.set_ylabel(r"Excess spread $\Omega-\Omega_{\mathrm{min}}$ (" + ANG2 + ")")

    h_arch = Line2D([], [], ls="none", marker=MARKER["ref"], ms=ms_a,
                    mfc="white", mec=col_a, mew=0.9)
    h_std = Line2D([], [], ls="none", marker=MARKER["aux2"], ms=ms_s,
                   mfc=col_s, mec=col_s, mew=0.7)
    h_rnd = Line2D([], [], ls="none", marker=MARKER["aux1"], ms=ms_r,
                   mfc=col_r, mec=col_r, mew=0.6)
    h_std_o = Line2D([], [], ls="none", marker=MARKER["aux2"], ms=ms_s,
                     mfc="white", mec=col_s, mew=0.7)
    h_rnd_o = Line2D([], [], ls="none", marker=MARKER["aux1"], ms=ms_r,
                     mfc="none", mec=col_r, mew=0.6)
    n_seed = len(chans[0]["seed"])
    ax.legend(
        [h_arch, h_std, h_rnd, (h_std_o, h_rnd_o)],
        ["wannier90, 100 iterations", "gpuw90, same starting gauge",
         f"gpuw90, {n_seed} random starting gauges", "not converged"],
        handler_map={tuple: HandlerTuple(ndivide=None, pad=0.6)},
        loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=4,
        columnspacing=1.6, handletextpad=0.4, borderaxespad=0.25)

    tot = {k: sum(ch[k] for ch in chans)
           for k in ("stray_arch", "stray_std", "stray_best")}
    prov("c", "stray WFs, total: archived / same start / best random",
         f"{tot['stray_arch']} / {tot['stray_std']} / {tot['stray_best']}",
         "multistart_summary.json")
    prov("c", "stray WFs summed over all random-start end points",
         int(sum(ch["seed_stray"].sum() for ch in chans)),
         "multistart_summary.json")
    prov("c", "channels in which no start converged",
         ", ".join(f"{ch['case']} spin {ch['spin']}" for ch in chans
                   if not ch["best_conv"]) or "none",
         "multistart_summary.json (best_converged)")
    return tot


# ---------------------------------------------- ghost data: check, not a panel
def ghost_check():
    """Why there is no ghost panel: print what the file does and does not
    support.  Nothing from here is drawn."""
    rows = read_tsv(GHOST_TSV)
    ref = json.loads(REF_JSON.read_text())
    ref = {name: (v.get("ref") or v.get("ref_mean") or {})
           for sect in ref.values() if isinstance(sect, dict)
           for name, v in sect.items() if isinstance(v, dict)}
    bench = {r["case"]: r for r in read_tsv(REL_TSV)
             if r["batch"] in ("wave1", "wave2")}
    changed_ghost, dev_changed, unscored, not_bench = [], [], [], []
    for r in rows:
        if r["ghost_v3"] != r["ghost_v31"]:
            changed_ghost.append(f"{r['case']} {r['ghost_v3']}->{r['ghost_v31']}")
        means = {m.group(1): float(m.group(2)) for m in
                 re.finditer(r"([A-Z][a-z]?)([+-][\d.]+)", r["flow_elem_avg_v31"])}
        target = ref.get(r["case"], {})
        if not target or r["old_dev"] == "-":
            unscored.append(r["case"])
        else:
            dev = max(abs(means[el] - t) for el, t in target.items())
            if abs(dev - float(r["old_dev"])) > 0.01:
                dev_changed.append(f"{r['case']} {r['old_dev']}->{dev:.2f}")
        if r["case"] not in bench or bench[r["case"]]["n_ghost"] != r["ghost_v3"]:
            not_bench.append(r["case"])
    prov("-", "ghost A/B: systems / with a changed ghost count",
         f"{len(rows)} / {len(changed_ghost)} ({', '.join(changed_ghost)})",
         rel(GHOST_TSV))
    prov("-", "ghost A/B: 'before' element means in the file", "none (only "
         "old_status, old_dev and the v3.1 means)", rel(GHOST_TSV))
    prov("-", "ghost A/B: dev of the v3.1 means differs from old_dev",
         ", ".join(dev_changed) or "none",
         f"{rel(GHOST_TSV)} + {rel(REF_JSON)}")
    prov("-", "ghost A/B: no reference, change not testable",
         ", ".join(unscored) or "none", rel(GHOST_TSV))
    prov("-", "ghost A/B: 'before' run is not the current benchmark entry",
         ", ".join(not_bench) or "none", f"{rel(GHOST_TSV)} vs {rel(REL_TSV)}")
    return dev_changed, not_bench


# ----------------------------------------------------------------------- main
def main():
    rel_rows = read_tsv(REL_TSV)
    beta_rows = read_tsv(BETA_TSV)
    chans = load_holdout()
    if chans is None:
        print(f"NOTICE: {HOLDOUT_A1} not found - panel (c) skipped "
              "(provisional data outside the repository).")

    w_cm = style.DOUBLE / CM
    top_h, bot_h = 4.30, 3.95                    # axes heights, cm
    h_cm = 0.55 + top_h + (1.95 + bot_h + 1.15 if chans else 1.25)
    fig = plt.figure(figsize=(w_cm * CM, h_cm * CM))
    top_y = h_cm - 0.55 - top_h
    ax_a = axes_cm(fig, 1.24, top_y, 8.30, top_h)
    ax_b = axes_cm(fig, 11.50, top_y, 5.56, top_h)

    got = panel_a(ax_a, rel_rows)
    got.update(panel_b(ax_b, beta_rows))
    style.panel_label(ax_a, "a", dx=-0.075)
    style.panel_label(ax_b, "b", dx=-0.115)
    if chans:
        ax_c = axes_cm(fig, 1.24, 1.15, 15.82, bot_h)
        got.update(panel_c(ax_c, chans))
        style.panel_label(ax_c, "c", dx=-0.039)
    ghost_check()

    style.save(fig, "figS2_reliability")

    # ---- provenance table
    print("\nProvenance (panel | quantity | value | source)")
    if chans:
        print(f"  panel (c) files are read from the PROVISIONAL directory\n"
              f"  {HOLDOUT_A1}")
    wq = max(len(p[1]) for p in PROV)
    for panel, what, value, source in PROV:
        print(f"  {panel:1s} | {what:<{wq}s} | {value} | {source}")
    print("\nCross-check against docs/RESULTS.md (stray totals: hold-out reports)")
    for key, label in (("n_ok", "agreeing cases"), ("n_not", "deviating cases"),
                       ("med_ok", "median R_elem, agreeing"),
                       ("med_not", "median R_elem, deviating"),
                       ("beta10_scan", "n_OK at default beta (scan)"),
                       ("stray_arch", "stray WFs, archived"),
                       ("stray_std", "stray WFs, same start"),
                       ("stray_best", "stray WFs, best random start")):
        if key not in got:
            print(f"  {label:<32s} not evaluated")
            continue
        ok = abs(float(got[key]) - float(DOC[key])) < 1e-6
        print(f"  {label:<32s} data {got[key]:g}   quoted {DOC[key]:g}   "
              f"{'OK' if ok else 'MISMATCH'}")
    print("\nWrote figS2_reliability.pdf / .png in", style.HERE)


if __name__ == "__main__":
    main()
