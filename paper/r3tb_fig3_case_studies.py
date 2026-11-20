#!/usr/bin/env python3
"""r3-T-B figure 3 - Case studies of the in-situ electronegativity assignment.

(a) Delithiation of Li2O2 (Li16O16 -> Li15O16 -> Li14O16): mean oxidation
    state of O and of Li against charge counting, and the oxidation state of
    every O atom.
(b) O vacancy in CeO2 (Ce4O8 -> Ce4O7): oxidation state of each Ce atom, and
    the four Ce atoms around the vacancy.
(c) The 13 held-out compounds with homonuclear bonds: atoms that receive the
    reference oxidation state.

Every calculation in the figure is free of unlocalised (stray) Wannier
functions; main() asserts it.

Layout: (a) over the full height on the left, (b) above (c) on the right;
the saved page is exactly 17.1 cm wide.

Run from this directory:   python3 r3tb_fig3_case_studies.py [--recompute]

"This work" is the frozen plugin pyleos/methods/frozen/r3-T-B-cap-repair.py,
called as the clean benchmark calls it (bench.load_plugin, bench.scrub,
plugin.assign).  All numbers come from r3tb_applications.py (same directory),
which runs the plugin, the min-cost flow and reads the greedy column of every
EOSCAR, writes results/r3tb_applications.tsv and caches the per-atom results
in _cache/r3tb/applications.json (md5 of every input, of the plugin and of
the pyleos sources; --recompute ignores the cache).  No number in the figure
is typed by hand.

Every panel is a function draw_<name>(ax) that fills the rectangle of the
axes it is given (stacked rows are inset axes), so that a single panel can be
redrawn elsewhere:

    import r3tb_fig3_case_studies as f3
    f3.draw_delithiation(ax); f3.draw_vacancy(ax); f3.draw_heldout(ax)

Marker and font sizes follow matplotlib's rcParams at the time of the call.

draw_hydrides(ax) (the eight MH4 hydrides against Comas-Vila et al., with
hydride_handles() for its legend) stays in the module but is not part of the
figure: four of the eight calculations carry stray WFs.  Its literature
values are apps.LIT_MH4.
Writes r3tb_fig3_case_studies.pdf / .png, prints every plotted number with
its source and keeps the printout in _cache/r3tb/fig3_case_studies.run.log.
"""
import argparse
import functools
import itertools
import sys
from collections import Counter

import numpy as np

import style
from style import C, MARKER, formula

import matplotlib as mpl                            # noqa: E402
import matplotlib.pyplot as plt                     # noqa: E402
from matplotlib.lines import Line2D                 # noqa: E402
from matplotlib.patches import Circle, Patch, Rectangle  # noqa: E402
from matplotlib.ticker import (FixedLocator, FuncFormatter,  # noqa: E402
                               MultipleLocator, NullLocator)

import r3tb_applications as apps                    # noqa: E402
from pyleos.parsers import parse_poscar             # noqa: E402

REPO = style.REPO
NAME = "r3tb_fig3_case_studies"
LOG = style.HERE / "_cache" / "r3tb" / "fig3_case_studies.run.log"

# ----------------------------------------------------- methods and labels --
KEY, BASE = "aux1", "nearest"          # style.py keys: this work, greedy
THIS, GREEDY, FLOW, REFC = C[KEY], C[BASE], C["flow"], C["ref"]
LABEL_THIS = "In-situ electronegativity (this work)"
LABEL_GREEDY = style.LABEL[BASE]
LABEL_FLOW = "Min-cost flow"
LIT_LABEL = "Comas-Vilà et al."
ONE_SITE_LABEL = "cell with one H site"
PLUGIN = "pyleos/methods/frozen/r3-T-B-cap-repair.py"

# ---------------------------------------------------------------- inputs --
LI2O2 = [("li16o16", "v31"), ("li15o16", "v31"), ("li14o16", "v31")]
MH4 = apps.MH4                          # cases/<name>/v3
CEO2 = [("ce4o8", "Ce4O8"), ("ce4o7", "Ce4O7")]      # cases/<name>/v3
HELDOUT_LABEL = {"s8": "S8"}            # formula of test_set.json is "S"

FORCE = False                           # set by --recompute
PROV = []                               # (panel, quantity, value, source)


def prov(panel, quantity, value, source):
    PROV.append((panel, quantity, str(value), str(source)))


@functools.lru_cache(maxsize=None)
def results():
    """The cached per-atom results of r3tb_applications.py."""
    return apps.load(FORCE)


def signed(v, _pos=None):
    v = int(round(v))
    return "0" if v == 0 else (f"+{v}" if v > 0 else f"−{-v}")


def minus(text):
    return text.replace("-", "−")


# ------------------------------------------------- panel (a): Li2O2 data ---
@functools.lru_cache(maxsize=None)
def delithiation_rows():
    runs = results()["runs"]
    n0 = None
    rows = []
    for name, run in LI2O2:
        rec = runs[f"{name}/{run}"]
        v = apps.li2o2_view(rec)
        n0 = v["n_li"] if n0 is None else n0
        row = {"name": name, "removed": n0 - v["n_li"], "n_o": v["n_o"],
               "ref": {"o": v["ref_o_mean"], "li": 1.0},
               "stray": rec["stray"], "n_wf": sum(rec["n_wf"])}
        for m, k in (("r3tb", "this"), ("greedy", "greedy"), ("flow", "flow")):
            row[k] = {"o": float(np.mean(v[m]["o"])),
                      "li": float(np.mean(v[m]["li"])),
                      "li_set": sorted(set(v[m]["li"])),
                      "o_atoms": dict(sorted(Counter(v[m]["o"]).items())),
                      "units": dict(sorted(Counter(v[m]["units"]).items()))}
        rows.append(row)
        src = f"{rec['dir']} (POSCAR, EOSCAR, wout)"
        prov("a", f"{name}: Li removed; mean OS(O) charge counting / this "
             "work / greedy", f"{row['removed']}; {row['ref']['o']:+.4f} / "
             f"{row['this']['o']:+.4f} / {row['greedy']['o']:+.4f}",
             src + "; reference -n(Li)/n(O)")
        prov("a", f"{name}: mean OS(Li) charge counting / this work / greedy",
             f"{row['ref']['li']:+.4f} / {row['this']['li']:+.4f} / "
             f"{row['greedy']['li']:+.4f}  (this work: all Li at "
             f"{row['this']['li_set']})", src)
        prov("a", f"{name}: O atoms by OS {{OS: atoms}}, this work | greedy",
             f"{row['this']['o_atoms']} | {row['greedy']['o_atoms']}", src)
        prov("a", f"{name}: (not drawn) O atoms by OS, flow; O2-unit charges "
             "{charge: units} this work | flow | greedy",
             f"{row['flow']['o_atoms']}; {row['this']['units']} | "
             f"{row['flow']['units']} | {row['greedy']['units']}", src)
        prov("a", f"{name}: (not drawn) stray WFs / N_WF", f"{row['stray']} "
             f"/ {row['n_wf']}", "assign_case(...).stray")
    return rows


# ----------------------------- hydride data (not in the default figure) ---
@functools.lru_cache(maxsize=None)
def hydride_rows():
    runs = results()["runs"]
    rows = []
    for name, metal in MH4:
        rec = runs[f"{name}/v3"]
        v = apps.mh4_view(rec, metal)
        v["name"] = name
        rows.append(v)
        lit, t, g = v["lit"], v["r3tb"], v["greedy"]
        if v["two_sites"]:
            h = ("Hb " + ",".join(f"{x:+d}" for x in t["hb"][:2]) + "; Ha-Ha "
                 + ",".join(f"{x:+d}" for x in t["ha_pairs"][0]))
            dist = f"Ha-Ha {v['d_ha']:.3f} A"
        else:
            h = "pairs " + " ".join("(" + ",".join(f"{x:+d}" for x in p) + ")"
                                    for p in t["pairs"])
            dist = f"all H {v['d_h']:.3f} A"
        prov("-", f"{name}: {metal} this work / lit / greedy; H atoms of one "
             "formula unit, this work; lit Hb, Ha; shortest H-H",
             f"{sorted(set(t['m']))} / {lit[0]:+d} / {sorted(set(g['m']))}; "
             f"{h}; {lit[1]:+d}, {lit[2]:+g}; {dist}",
             f"{rec['dir']}; notes N{lit[3]}, p. {lit[4]}")
        prov("-", f"{name}: (not drawn) mean OS(H) this work / lit / flow; "
             "stray WFs / N_WF; |dNe|",
             f"{t['h_mean']:+.2f} / {(2 * lit[1] + 2 * lit[2]) / 4:+.2f} / "
             f"{v['flow']['h_mean']:+.2f}; {v['stray']} / {v['n_wf']}; "
             f"{v['dne']:.1e}", rec["dir"])
    return rows


# ---------------------------------------------------- panel (b): CeO2 data --
SHIFTS = np.array(list(itertools.product((-1, 0, 1), repeat=3)), float)


def min_image_vec(st, frac_from, frac_to):
    """Shortest periodic vector (A) between two fractional positions."""
    d = np.array(frac_to, float) - np.array(frac_from, float)
    d -= np.round(d)
    v = (d + SHIFTS) @ np.array(st.lattice)
    return v[np.argmin(np.linalg.norm(v, axis=1))]


@functools.lru_cache(maxsize=None)
def vacancy_data():
    runs = results()["runs"]
    out = {}
    for name, _ in CEO2:
        rec = runs[f"{name}/v3"]
        ce = apps.idx(rec["elements"], "Ce")
        ox = apps.idx(rec["elements"], "O")
        out[name] = {"ce_atoms": [i + 1 for i in ce],
                     "this": [rec["r3tb"][i] for i in ce],
                     "flow": [rec["flow"][i] for i in ce],
                     "greedy": [rec["greedy"][i] for i in ce],
                     "o_this": sorted({rec["r3tb"][i] for i in ox}),
                     "stray": rec["stray"], "n_wf": sum(rec["n_wf"])}
        d = out[name]
        prov("b", f"{name}: Ce1-Ce4 this work | (not drawn) flow | greedy; "
             "O, this work; stray WFs / N_WF",
             f"{d['this']} | {d['flow']} | {d['greedy']}; {d['o_this']}; "
             f"{d['stray']} / {d['n_wf']}", rec["dir"])
    # vacancy = the O of ce4o8 that has no counterpart in ce4o7
    full = parse_poscar(REPO / "cases" / "ce4o8" / "v3" / "POSCAR")
    vac = parse_poscar(REPO / "cases" / "ce4o7" / "v3" / "POSCAR")
    missing = [full.frac[i] for i in apps.idx(full.elements, "O")
               if all(np.linalg.norm(min_image_vec(full, full.frac[i],
                                                   vac.frac[j])) > 0.2
                      for j in apps.idx(vac.elements, "O"))]
    assert len(missing) == 1, missing
    vec = [min_image_vec(vac, missing[0], vac.frac[i])
           for i in apps.idx(vac.elements, "Ce")]
    out["vec_ce"] = [list(map(float, x)) for x in vec]
    out["d_ce"] = [float(np.linalg.norm(x)) for x in vec]
    prov("b", "Ce-vacancy distances (A); projected length in the sketch (A)",
         f"{[round(x, 3) for x in out['d_ce']]}; "
         f"{sorted({round(float(np.hypot(x[0], x[1])), 2) for x in vec})}",
         "POSCAR of cases/ce4o8/v3 and cases/ce4o7/v3")
    return out


# ------------------------------------------- panel (c): held-out set data --
@functools.lru_cache(maxsize=None)
def heldout_rows():
    homo, runs = results()["homonuclear"], results()["runs"]
    rows, tot = [], {m: [0, 0, 0] for m in ("this", "flow", "greedy")}
    for name, h in homo.items():
        test = apps.bench.TEST["test_homonuclear"]["cases"][name]
        row = {"name": name, "label": HELDOUT_LABEL.get(name, test["formula"]),
               "n": h["r3tb"]["atoms_labelled"],
               "stray": runs[h["dir"].split("cases/", 1)[1]]["stray"]}
        for m, k in (("r3tb", "this"), ("flow", "flow"), ("greedy", "greedy")):
            row[k] = h[m]["atoms_ok"]
            tot[k][0] += h[m]["atoms_ok"]
            tot[k][1] += h[m]["atoms_labelled"]
            tot[k][2] += h[m]["atoms_ok"] == h[m]["atoms_labelled"]
        rows.append(row)
        prov("c", f"{name}: atoms correct this work / flow / greedy, of "
             "labelled; wrong sites of this work; stray WFs",
             f"{row['this']} / {row['flow']} / {row['greedy']}, of "
             f"{row['n']}; {h['r3tb']['wrong'] or '-'}; {row['stray']}",
             f"{h['dir']}; labels pyleos/sets/test_set.json")
    for k, lab in (("this", "this work"), ("flow", "flow"),
                   ("greedy", "greedy")):
        prov("c", f"total, {lab}: atoms correct / labelled; compounds with "
             "every atom correct", f"{tot[k][0]}/{tot[k][1]}; {tot[k][2]}/"
             f"{len(rows)}", "bench.score_case over the 13 cases"
             if k != "greedy" else "OS(int) of the EOSCAR against "
             "bench.site_labels")
    return rows, tot


# ================================================================ drawing ==
def ms(k=1.0):
    """Marker size relative to rcParams['lines.markersize']."""
    return mpl.rcParams["lines.markersize"] * k


def small():
    return mpl.rcParams["xtick.labelsize"]


def size_cm(ax):
    """Width and height of an axes on the page, in cm."""
    p, fig = ax.get_position(), ax.figure
    return (p.width * fig.get_figwidth() / style.CM,
            p.height * fig.get_figheight() / style.CM)


def ylabel(ax, text, y=0.5, offset=0.80):
    """Rotated y label, its centre `offset` cm left of the axes rectangle (at
    7 pt tick labels; the distance scales with the tick-label size) and at
    the height y (axes fraction); works for a parent with its axis off."""
    w, _ = size_cm(ax)
    offset *= mpl.rcParams["ytick.labelsize"] / 7.0
    ax.text(-offset / w, y, text, transform=ax.transAxes, rotation=90,
            ha="center", va="center", fontsize=mpl.rcParams["axes.labelsize"])


def this_kw(**kw):
    out = dict(color=THIS, marker=MARKER[KEY], ms=ms(1.05), mew=0, lw=0,
               zorder=5)
    out.update(kw)
    return out


def greedy_kw(**kw):
    out = dict(color=GREEDY, marker=MARKER[BASE], ms=ms(1.1), mew=0, lw=0,
               zorder=3)
    out.update(kw)
    return out


def ref_kw(**kw):
    out = dict(color=REFC, marker=MARKER["ref"], ms=ms(1.75), mfc="none",
               mew=0.8, lw=0, zorder=2)
    out.update(kw)
    return out


# --- panel (a) ---------------------------------------------------------------
A_YLABEL = 1.22          # cm from the axes to the centre of the y label


def draw_delithiation(ax):
    """Li2O2 delithiation: mean OS of O and Li, and OS of every O atom."""
    rows = delithiation_rows()
    ax.axis("off")
    f_o, f_li, f_at, gap = 0.31, 0.15, 0.48, 0.03       # heights, top down
    ax_o = ax.inset_axes([0, 1 - f_o, 1, f_o])
    ax_li = ax.inset_axes([0, 1 - f_o - gap - f_li, 1, f_li], sharex=ax_o)
    ax_at = ax.inset_axes([0, 0, 1, f_at], sharex=ax_o)
    x = np.array([r["removed"] for r in rows], float)

    # mean oxidation states: charge counting, greedy, this work (nested)
    for a, el in ((ax_o, "o"), (ax_li, "li")):
        a.plot(x, [r["ref"][el] for r in rows],
               **ref_kw(ms=ms(2.2), lw=0.8, ls=(0, (4, 2))))
        a.plot(x, [r["greedy"][el] for r in rows], **greedy_kw(ms=ms(1.45)))
        a.plot(x, [r["this"][el] for r in rows], **this_kw(lw=1.0))
    o_ref = [r["ref"]["o"] for r in rows]
    ax_o.set_ylim(min(o_ref) - 0.045, max(o_ref) + 0.045)
    ax_o.yaxis.set_major_locator(FixedLocator(o_ref))
    ax_o.yaxis.set_minor_locator(NullLocator())
    ax_o.yaxis.set_major_formatter(FuncFormatter(
        lambda v, _p: minus(f"{v:.4f}".rstrip("0").rstrip("."))))
    ax_o.text(0.035, 0.90, "O, mean", transform=ax_o.transAxes, ha="left",
              va="top")
    ax_o.legend(handles=[Line2D([], [], label="Charge counting",
                                **ref_kw(ms=ms(1.75), lw=0.8,
                                         ls=(0, (4, 2))))],
                loc="lower right", handlelength=2.6, borderaxespad=0.5)
    ax_li.set_ylim(-0.45, 1.75)
    ax_li.yaxis.set_major_locator(FixedLocator([0, 1]))
    ax_li.yaxis.set_minor_locator(NullLocator())
    ax_li.yaxis.set_major_formatter(FuncFormatter(signed))
    ax_li.text(0.035, 0.14, "Li, mean", transform=ax_li.transAxes, ha="left",
               va="bottom")

    # every O atom: one symbol per oxidation state, the number of atoms
    # printed next to it; the symbol grows with that number
    for r in rows:
        for k, kw, side, dx in (("this", this_kw, -1, -0.13),
                                ("greedy", greedy_kw, +1, 0.20)):
            for os_, n in r[k]["o_atoms"].items():
                size = ms(0.62 + 0.34 * np.sqrt(n))
                ax_at.plot([r["removed"] + dx], [os_], **kw(ms=size))
                ax_at.annotate(str(n), (r["removed"] + dx, os_),
                               xytext=(side * (0.5 * size
                                               + 0.31 * small()), 0),
                               textcoords="offset points",
                               ha="right" if side < 0 else "left",
                               va="center_baseline", fontsize=small())
    levels = sorted({os_ for r in rows for k in ("this", "greedy")
                     for os_ in r[k]["o_atoms"]})
    ax_at.set_ylim(levels[0] - 0.75, levels[-1] + 1.75)   # room for the text
    ax_at.yaxis.set_major_locator(FixedLocator(
        range(levels[0], levels[-1] + 1)))
    ax_at.yaxis.set_minor_locator(NullLocator())
    ax_at.yaxis.set_major_formatter(FuncFormatter(signed))
    ax_at.text(0.035, 0.955, "O atoms (number at each value)",
               transform=ax_at.transAxes, ha="left", va="top")

    ax_o.set_xlim(-0.62, 2.62)
    ax_o.set_xticks(x)
    for a in (ax_o, ax_li, ax_at):
        a.xaxis.set_minor_locator(NullLocator())
    ax_o.tick_params(labelbottom=False)
    ax_li.tick_params(labelbottom=False)
    ax_at.set_xticklabels([f"{int(v)}" for v in x])
    ax_at.set_xlabel("Li removed per " + formula("Li16O16") + " cell")
    ylabel(ax, "Mean oxidation state", y=1 - 0.5 * (f_o + gap + f_li),
           offset=A_YLABEL)
    ylabel(ax, "Oxidation state", y=0.5 * f_at, offset=A_YLABEL)


# --- hydrides (not in the default figure) ------------------------------------
def hydride_handles():
    """Legend handles that draw_hydrides needs besides the methods."""
    return [Line2D([], [], label=LIT_LABEL, **ref_kw()),
            Patch(facecolor=C["light"], alpha=0.45, lw=0,
                  label=ONE_SITE_LABEL)]


def draw_hydrides(ax):
    """MH4 hydrides: metal, H atoms of one formula unit, shortest H-H.
    Dagger: stray WFs present; double dagger: |dNe| above the warning level
    of the EOSCAR."""
    rows = hydride_rows()
    ax.axis("off")
    f_m, f_h, f_d, gap = 0.42, 0.37, 0.17, 0.02
    ax_m = ax.inset_axes([0, 1 - f_m, 1, f_m])
    ax_h = ax.inset_axes([0, f_d + gap, 1, f_h], sharex=ax_m)
    ax_d = ax.inset_axes([0, 0, 1, f_d], sharex=ax_m)
    lit = ref_kw()
    xb, xa, dx = -0.22, 0.22, 0.078       # H_b and H_a slots, atom offset
    for n, r in enumerate(rows):
        t, g, lv = r["r3tb"], r["greedy"], r["lit"]
        if not r["two_sites"]:         # cell is not the literature structure
            for a in (ax_m, ax_h, ax_d):
                a.axvspan(n - 0.5, n + 0.5, color=C["light"], alpha=0.45,
                          lw=0, zorder=0)
        # metal
        ax_m.plot([n], [lv[0]], **lit)
        ax_m.plot([n] * len(set(t["m"])), sorted(set(t["m"])), **this_kw())
        if set(g["m"]) != set(t["m"]):
            # one square per distinct value among the metal atoms of the cell
            gm = sorted(set(g["m"]))
            ax_m.plot([n] * len(gm), gm, **greedy_kw())
        # hydrogen: literature per site, this work per atom of a formula unit
        ax_h.plot([n + xb, n + xa], [lv[1], lv[2]], **lit)
        if r["two_sites"]:
            ax_h.plot([n + xb - dx, n + xb + dx], t["hb"][:2], **this_kw())
            ax_h.plot([n + xa - dx, n + xa + dx], t["ha_pairs"][0],
                      **this_kw(lw=0.8 if r["h2"] else 0))
        else:
            ax_h.plot(n + np.array([xb - dx, xb + dx, xa - dx, xa + dx]),
                      t["pairs"][0] + t["pairs"][1], **this_kw())
        ax_d.plot([n], [r["d_min"]], color=C["ink"], marker="o", ms=ms(0.95),
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
    ax_h.text(0.015, 0.92, "H atoms", transform=ax_h.transAxes, ha="left",
              va="top")
    # name the two H sites once, above a compound in which both are at -1
    k = next(k for k, r in enumerate(rows) if r["two_sites"]
             and r["lit"][2] == -1)
    for off, name in ((xb, "H$_b$"), (xa, "H$_a$")):
        ax_h.annotate(name, (k + off, -1), xytext=(0, 0.86 * small()),
                      textcoords="offset points", ha="center", va="bottom",
                      fontsize=small())
    ax_d.set_ylim(0.35, 2.25)
    ax_d.yaxis.set_major_locator(MultipleLocator(1))
    ax_d.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax_m.set_xlim(-0.55, len(rows) - 0.45)
    ax_m.set_xticks(range(len(rows)))
    ax_d.set_xticklabels([formula(r["metal"] + "H4")
                          + ("$^{\\dagger}$" if r["stray"] else "")
                          + ("$^{\\ddagger}$" if r["lost"] else "")
                          for r in rows])
    for a in (ax_m, ax_h, ax_d):
        a.xaxis.set_minor_locator(NullLocator())
        a.tick_params(axis="x", which="both", top=False)
    ax_m.tick_params(labelbottom=False)
    ax_h.tick_params(labelbottom=False)
    ax_d.tick_params(axis="x", pad=1.5 * small())
    for lab in ax_d.get_xticklabels():
        lab.set_verticalalignment("baseline")
    ylabel(ax, "Oxidation state", y=1 - 0.5 * (f_m + gap + f_h))
    ylabel(ax, r"$d_{\mathrm{H–H}}$ (Å)", y=0.5 * f_d)


# --- panel (b) ---------------------------------------------------------------
def draw_vacancy(ax):
    """O vacancy in CeO2: OS of each Ce atom, and the vacancy with its four
    Ce neighbours seen along [001]."""
    d = vacancy_data()
    ax.axis("off")
    f_main, gap = 0.625, 0.03
    ax_m = ax.inset_axes([0, 0, f_main, 1])
    ax_s = ax.inset_axes([f_main + gap, 0, 1 - f_main - gap, 1])
    xs = {"ce4o8": np.arange(4), "ce4o7": np.arange(4) + 4.6}
    for name, lab in CEO2:
        x = xs[name]
        ax_m.plot(x, d[name]["this"], **this_kw(ms=ms(1.2)))
        n3 = d[name]["this"].count(3)
        ax_m.text(float(np.mean(x)), 4.56, formula(lab), ha="center",
                  va="center")
        ax_m.text(float(np.mean(x)), 3.5,
                  ("no" if n3 == 0 else str(n3)) + r" Ce$^{3{+}}$",
                  ha="center", va="center", fontsize=small())
    ax_m.axvline(3.8, color=C["ink"], lw=0.6)
    ax_m.set_xlim(-0.6, 8.2)
    ax_m.set_xticks(list(xs["ce4o8"]) + list(xs["ce4o7"]))
    ax_m.set_xticklabels([str(n) for n in d["ce4o8"]["ce_atoms"]]
                         + [str(n) for n in d["ce4o7"]["ce_atoms"]])
    ax_m.set_xlabel("Ce atom")
    ax_m.xaxis.set_minor_locator(NullLocator())
    ax_m.tick_params(axis="x", which="both", top=False)
    ax_m.set_ylim(2.65, 4.85)
    ax_m.yaxis.set_major_locator(FixedLocator([3, 4]))
    ax_m.yaxis.set_minor_locator(NullLocator())
    ax_m.yaxis.set_major_formatter(FuncFormatter(signed))
    ylabel(ax, "Oxidation state")

    # sketch: the vacancy and its four Ce neighbours, seen along [001]
    w, h = size_cm(ax_s)
    reach = max(float(np.hypot(v[0], v[1])) for v in d["vec_ce"]) + 0.70
    scale = min(w / (2 * reach + 0.1), h / (2 * reach + 1.0))   # cm per A
    for n, vec in enumerate(d["vec_ce"]):
        ax_s.plot([0, vec[0]], [0, vec[1]], color=C["ink"], lw=0.6,
                  ls=(0, (2, 1.5)), zorder=1)
        ax_s.add_patch(Circle(vec[:2], 0.70, facecolor=C["light"],
                              edgecolor=C["ink"], lw=0.5, zorder=2))
        ax_s.text(vec[0], vec[1], f"Ce{n + 1}", ha="center",
                  va="center_baseline", fontsize=small(), zorder=3)
    ax_s.add_patch(Circle((0, 0), 0.60, facecolor="white",
                          edgecolor=C["ink"], lw=0.6, ls=(0, (2, 1.5)),
                          zorder=2))
    ax_s.text(0, 0, r"V$_{\mathrm{O}}$", ha="center", va="center_baseline",
              fontsize=small(), zorder=3)
    dist = sorted(set(round(v, 2) for v in d["d_ce"]))
    assert len(dist) == 1, dist
    ymax = 0.5 * h / scale
    ax_s.text(0, -ymax + 0.10, f"4 × {dist[0]:.2f} Å", ha="center",
              va="bottom", fontsize=small())
    ax_s.set_xlim(-0.5 * w / scale, 0.5 * w / scale)
    ax_s.set_ylim(-ymax - 0.25, ymax - 0.25)
    ax_s.set_aspect("equal")
    ax_s.axis("off")


# --- panel (c) ---------------------------------------------------------------
def draw_heldout(ax):
    """Held-out compounds with homonuclear bonds: atoms that receive the
    reference oxidation state, per compound, for three methods."""
    rows, tot = heldout_rows()
    x = np.arange(len(rows), dtype=float)
    series = (("this", THIS, MARKER[KEY], -0.25, 5),
              ("flow", FLOW, MARKER["flow"], 0.0, 4),
              ("greedy", GREEDY, MARKER[BASE], 0.25, 3))
    handles = []
    for k, col, mk, off, z in series:
        y = np.array([100.0 * r[k] / r["n"] for r in rows])
        ax.vlines(x + off, 0, y, color=col, lw=0.8, zorder=z - 0.5)
        ax.plot(x + off, y, color=col, marker=mk, ms=ms(1.0), mew=0, lw=0,
                zorder=z, clip_on=False)
        handles.append(Line2D([], [], color=col, marker=mk, ms=ms(1.0),
                              mew=0, lw=0,
                              label=f"{tot[k][0]}/{tot[k][1]}"))
    ax.set_ylim(0, 112)
    ax.yaxis.set_major_locator(FixedLocator([0, 50, 100]))
    ax.yaxis.set_minor_locator(FixedLocator([25, 75]))
    ax.set_xlim(-0.6, len(rows) - 0.4)
    ax.set_xticks(x)
    ax.set_xticklabels([formula(r["label"]) for r in rows], rotation=45,
                       ha="right", rotation_mode="anchor")
    ax.xaxis.set_minor_locator(NullLocator())
    ax.tick_params(axis="x", which="both", top=False)
    ylabel(ax, "Atoms correct (%)")
    leg = ax.legend(handles=handles, loc="lower right",
                    bbox_to_anchor=(1.0, 1.0), ncol=3, columnspacing=1.2,
                    handletextpad=0.2, borderaxespad=0.25,
                    title="Atoms correct, 13 compounds:")
    leg._legend_box.align = "right"
    leg.get_title().set_fontsize(mpl.rcParams["legend.fontsize"])


# ----------------------------------------------------------------- page ----
# All positions are in centimetres on the page, measured from its top-left
# corner.  The canvas is the page minus the pad that style.save adds on every
# side, so that the saved file is exactly PAGE_W wide.
PAGE_W = style.DOUBLE / style.CM       # 17.1 cm
PAD = 0.02 / style.CM                  # pad_inches of style.save, in cm
FIG_W = PAGE_W - 2 * PAD

TOP = 0.95                             # top edge of (a) and (b)
A_L, A_W = 1.80, 5.70                  # left column: (a)
B_L, B_W = 9.45, 7.45                  # right column: (b) above (c)
B_H = 2.35                             # height of (b)
C_TOP = TOP + B_H + 1.80               # labels of (b), then the legend of (c)
C_H = 2.30                             # height of (c)
A_H = C_TOP + C_H - TOP                # (a) ends where (c) ends
PAGE_H = C_TOP + C_H + 1.02            # room for the rotated labels of (c)
FIG_H = PAGE_H - 2 * PAD
LABEL_DX = -0.10                       # panel label: right end from axes edge


def ax_cm(fig, left, top, width, height, **kw):
    """Axes placed in centimetres measured from the top-left corner."""
    return fig.add_axes([(left - PAD) / FIG_W,
                         1 - (top - PAD + height) / FIG_H,
                         width / FIG_W, height / FIG_H], **kw)


def label_panel(ax, letter, x_cm=LABEL_DX, y_cm=0.10):
    """style.panel_label with its right end x_cm from the left axes edge and
    its baseline y_cm above the top edge."""
    w, h = size_cm(ax)
    style.panel_label(ax, letter, dx=x_cm / w, dy=1 + y_cm / h)


def method_legend(fig):
    handles = [
        Line2D([], [], label=LABEL_THIS, **this_kw()),
        Line2D([], [], label=LABEL_GREEDY, **greedy_kw()),
        Line2D([], [], color=FLOW, marker=MARKER["flow"], ms=ms(1.0), mew=0,
               lw=0, label=LABEL_FLOW),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=len(handles),
               bbox_to_anchor=(0.5, 1.0), columnspacing=2.4,
               handletextpad=0.4, borderaxespad=0.15, handlelength=1.3,
               handleheight=0.8)


def print_provenance():
    w = [max(len(r[i]) for r in PROV) for i in range(3)]
    w[1], w[2] = min(w[1], 78), min(w[2], 70)
    print("\nPROVENANCE OF THE PLOTTED VALUES (panel | quantity | value | "
          "source)")
    print("-" * 170)
    for p, q, v, s in PROV:
        print(f"({p}) {q:<{w[1]}}  {v:<{w[2]}}  {s}")
    print("-" * 170)


def main():
    global FORCE
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--recompute", action="store_true",
                    help="ignore _cache/r3tb/applications.json and rerun "
                         "r3-T-B, the flow and the greedy read-out")
    args = ap.parse_args()
    FORCE = args.recompute
    tee = apps.Tee(sys.stdout)
    sys.stdout = tee
    try:
        data = results()
        print(f"this work = {PLUGIN} via bench.load_plugin + bench.scrub + "
              f"plugin.assign; tables read by the plugin: "
              f"{data['tables_used_r3tb'] or 'none'}")
        fig = plt.figure(figsize=(FIG_W * style.CM, FIG_H * style.CM))
        # invisible frame: pins the tight bounding box of style.save
        fig.add_artist(Rectangle((0, 0), 1, 1, transform=fig.transFigure,
                                 fill=False, edgecolor="none", lw=0))
        method_legend(fig)
        panels = (("a", draw_delithiation, A_L, TOP, A_W, A_H),
                  ("b", draw_vacancy, B_L, TOP, B_W, B_H),
                  ("c", draw_heldout, B_L, C_TOP, B_W, C_H))
        for letter, draw, left, top, width, height in panels:
            ax = ax_cm(fig, left, top, width, height)
            draw(ax)
            label_panel(ax, letter)
        style.save(fig, NAME)
        # one standard for everything drawn: no unlocalised WFs
        stray = ([(r["name"], r["stray"]) for r in delithiation_rows()]
                 + [(n, vacancy_data()[n]["stray"]) for n, _ in CEO2]
                 + [(r["name"], r["stray"]) for r in heldout_rows()[0]])
        assert all(n == 0 for _, n in stray), stray
        prov("a-c", "stray WFs in the calculations drawn",
             f"0 in each of the {len(stray)} calculations",
             "len(assign_case(load_case(dir)).stray)")
        print_provenance()
        print(f"wrote {NAME}.pdf / {NAME}.png (page {PAGE_W:.1f} x "
              f"{PAGE_H:.2f} cm)")
    finally:
        sys.stdout = tee.stream
        LOG.parent.mkdir(parents=True, exist_ok=True)
        LOG.write_text(f"$ python3 {NAME}.py"
                       + (" --recompute" if FORCE else "") + "\n"
                       + tee.buf.getvalue())


if __name__ == "__main__":
    main()
