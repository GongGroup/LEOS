#!/usr/bin/env python3
"""Slide-scale panels of the three r3-T-B figures (16:9 talk, 13.33 x 7.5 in).

Run from this directory:   python3 r3tb_slide_panels.py        (about 10 s)

Writes _cache/r3tb/v2_panels/<name>.png (300 dpi, white background, tight
bounding box, no panel letter) and _cache/r3tb/v2_panels/INDEX.md (natural
size, aspect ratio, content and numbers of every image).

Every panel is drawn by the panel function of its journal figure script
(r3tb_fig1_method.py, r3tb_fig2_benchmark.py, r3tb_fig3_case_studies.py in
this directory), which take their numbers from the frozen plugin, results/
and the caches under _cache/r3tb/ (see the caption files).  Nothing is typed
in here; the numbers listed in INDEX.md are read from the objects that the
panel functions plot.  The caches are only read: if one of them is stale the
script stops and names the command that rebuilds it.

Scale.  A panel is laid out at journal scale (style.py: 7.5 pt text, 7 pt
ticks and legends, 8 pt axis labels) and the finished page is magnified M
times when it is written: savefig(dpi = 300 M) with the resolution tag of
the PNG set to 300 dpi.  Fonts, markers, line weights and dashes therefore
scale together, and the image has its natural size, M x the layout, at
300 dpi.  M = 2 gives 14 / 15 / 16 pt.  save() refuses a panel whose smallest
text is below 12 pt at natural size.

Slide-specific choices (data and drawing code are those of the figure scripts):
  legends         a panel that relied on a legend of its journal figure gets
                  its own, with the colours and markers of style.py
  dividing_point  drawn at its design width (text at journal size), then the
                  axes are widened; the labels outside the frame keep their
                  design distance from it
  accuracy        bars keep the journal ratio of thickness to row pitch
  accuracy_atoms  draw_accuracy with its second and third block removed
  margins         31 rows at 14 pt do not fit the height of a slide in one
                  column: the two blocks of draw_margins are shown side by
                  side (the function is drawn twice, each axes shows one
                  block); margins_column is the journal layout at the
                  largest scale that fits 7.5 in
  margins_summary the ratios of draw_margins in two rows, one per deciding
                  rule: one marker per bond pair, pairs with the same ratio
                  stacked about the row (in order of the ratio, each at the
                  free place nearest to the centre line; the stacks of a row
                  are compressed by one factor to fit it), the median of
                  each row marked; drawn here, not by a panel function

python3 r3tb_slide_panels.py NAME [NAME ...] writes only the named images and
leaves INDEX.md as it is.
"""
import json
import math
import struct
import sys
from collections import Counter
from pathlib import Path

sys.dont_write_bytecode = True      # leave no __pycache__ in the repository

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import numpy as np                                    # noqa: E402
import style                                          # noqa: E402
from style import C, CM, MARKER                       # noqa: E402
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.axes import Axes                      # noqa: E402
from matplotlib.lines import Line2D                   # noqa: E402
from matplotlib.patches import Patch, Rectangle       # noqa: E402
from matplotlib.ticker import FixedLocator, NullFormatter  # noqa: E402
from matplotlib.transforms import Bbox                # noqa: E402

import r3tb_fig1_method as f1                         # noqa: E402
import r3tb_fig2_benchmark as f2                      # noqa: E402
import r3tb_fig3_case_studies as f3                   # noqa: E402

apps = f3.apps
OUT = HERE / "_cache" / "r3tb" / "v2_panels"
DPI = 300
PAD = 0.03              # inches around the tight box, at journal scale
PT_MIN = 12.0           # smallest text allowed at natural size
SLIDE = (13.33, 7.5)    # inches

THIS = C["aux1"]        # in-situ electronegativity (this work)


# ------------------------------------------------------------------ caches
def check_caches():
    """The panel functions read three caches and rebuild a stale one in
    place.  This script only reads: stop if a cache is stale."""
    stale = []
    sha = f2.sha256(f2.FROZEN)
    an = json.loads(f2.ANALYSIS_JSON.read_text()) if f2.ANALYSIS_JSON.is_file() else {}
    sc = json.loads(f2.SCAN_JSON.read_text()) if f2.SCAN_JSON.is_file() else {}
    grid = all([p["value"] for p in sc.get("scan", {}).get(k, {}).get("points", [])]
               == f2.SCAN[k]["values"] for k in f2.SCAN)
    if an.get("plugin_sha256") != sha or sc.get("plugin_sha256") != sha or not grid:
        stale.append("python3 r3tb_fig2_benchmark.py --recompute")
    inputs = []
    for _study, case, run, _role, _sid, _stat in apps.case_list():
        inputs += apps.case_files(apps.CASES / case / run)
    inputs += [apps.CASES / "INDEX.tsv", apps.CASES / "STRUCTURES.tsv"]
    sig = {apps.rel(p): apps.md5(p) for p in inputs + apps.SOURCES}
    sig["cache_version"] = apps.CACHE_VERSION
    blob = json.loads(apps.CACHE.read_text()) if apps.CACHE.is_file() else {}
    if blob.get("inputs") != sig:
        stale.append("python3 r3tb_applications.py")
    if stale:
        sys.exit("stale cache under _cache/r3tb/; rebuild first with: "
                 + "; ".join(stale))


# ----------------------------------------------------------------- helpers
def figure(w_cm, h_cm):
    return plt.figure(figsize=(w_cm * CM, h_cm * CM))


def axes(fig, left, top, width, height):
    """Axes placed in cm from the left / top edges of the figure."""
    fw, fh = fig.get_size_inches() / CM
    return fig.add_axes([left / fw, 1 - (top + height) / fh,
                         width / fw, height / fh])


def legend(fig, handles, x, y, loc, **kw):
    """Figure legend anchored x cm from the left and y cm from the top."""
    fw, fh = fig.get_size_inches() / CM
    opts = dict(borderaxespad=0.0, borderpad=0.0, handletextpad=0.4,
                handlelength=1.3, columnspacing=1.8, labelspacing=0.4)
    opts.update(kw)
    return fig.legend(handles=handles, loc=loc,
                      bbox_to_anchor=(x / fw, 1 - y / fh), **opts)


def drawn_texts(fig):
    """Every text that is drawn: free texts, legends, labels of visible axes."""
    fig.canvas.draw()
    out, legends, stack = list(fig.texts), list(fig.legends), list(fig.axes)
    while stack:
        ax = stack.pop()
        stack += ax.child_axes
        out += ax.texts
        if ax.get_legend() is not None:
            legends.append(ax.get_legend())
        if ax.axison:
            for axis in (ax.xaxis, ax.yaxis):
                out.append(axis.label)
                for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                    out += [tick.label1, tick.label2]
    for leg in legends:
        out += leg.get_texts() + [leg.get_title()]
    return [t for t in out if t.get_visible() and t.get_text().strip()]


def content_box(fig):
    """Bounding box (inches) of everything that is drawn.  Unlike
    bbox_inches="tight" it leaves out the rectangle of an axes whose axis is
    off (the frames of draw_accuracy, draw_rule, draw_molecule, ...)."""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()

    def boxes(ax):
        if ax.axison:
            return [ax.get_tightbbox(rend)]
        hidden = [ax.patch, ax.xaxis, ax.yaxis, *ax.spines.values()]
        out = []
        for a in ax.get_children():
            if a.get_visible() and not any(a is h for h in hidden):
                out += (boxes(a) if isinstance(a, Axes)
                        else [a.get_tightbbox(rend)])
        return out

    found = [a.get_window_extent(rend) for a in fig.texts + fig.legends]
    for ax in fig.axes:
        found += boxes(ax)
    found = [b for b in found if b is not None and 0 < b.width < float("inf")
             and 0 < b.height < float("inf")]
    return Bbox.union(found).transformed(fig.dpi_scale_trans.inverted())


def save(fig, name, mag):
    """Write <name>.png magnified `mag` times; return its natural size."""
    OUT.mkdir(parents=True, exist_ok=True)
    pt = sorted({round(t.get_fontsize() * mag, 1) for t in drawn_texts(fig)})
    assert pt[0] >= PT_MIN, (name, pt)
    path = OUT / f"{name}.png"
    fig.savefig(path, dpi=DPI * mag, bbox_inches=content_box(fig).padded(PAD),
                facecolor="white", pil_kwargs={"dpi": (DPI, DPI)})
    plt.close(fig)
    with path.open("rb") as fh:
        w, h = struct.unpack(">II", fh.read(24)[16:24])
    assert w / DPI <= SLIDE[0] and h / DPI <= SLIDE[1], (name, w, h)
    return dict(file=path.name, px=(w, h), w=w / DPI, h=h / DPI, mag=mag,
                pt=(pt[0], pt[-1]))


def plain(s):
    """Label text without mathtext markup, for INDEX.md."""
    return (s.replace(r"$\beta$", "β").replace(r"$_{", "").replace("}$", "")
             .replace("$", ""))


def name_of(case):
    return f2.NAME.get(case.lower(), case)


def signed(v):
    return f3.signed(v)


def marker(key, label, ms=4.0):
    """Legend handle of a method: colour and marker of style.py."""
    return Line2D([], [], ls="none", marker=MARKER[key], ms=ms, color=C[key],
                  mec=C[key], label=label)


def case_study_handles(*keys):
    """Legend handles of the case-study panels (as f3.method_legend)."""
    made = {
        "this": Line2D([], [], label=f3.LABEL_THIS, **f3.this_kw()),
        "greedy": Line2D([], [], label=f3.LABEL_GREEDY, **f3.greedy_kw()),
        "flow": Line2D([], [], color=f3.FLOW, marker=MARKER["flow"],
                       ms=f3.ms(1.0), mew=0, lw=0, label=f3.LABEL_FLOW)}
    return [made[k] for k in keys]


# ------------------------------------------------------- method (figure 1)
def panel_molecule():
    fig = figure(f1.WA, f1.H)
    f1.draw_molecule(axes(fig, 0, 0, f1.WA, f1.H))
    mol = f1.molecule()
    n = Counter(mol["kinds"])
    return fig, 2.0, dict(
        shows="Atoms and Wannier-function centres of one H2S2 molecule: bond "
              "pairs (contested, filled) and lone pairs (own, open); own legend.",
        printed=[f"atom labels {', '.join(mol['labels'])}; formula "
                 f"{mol['formula']}; scale bar 1 Å"],
        drawn=[f"{len(mol['kinds'])} WF centres = {n['bond']} bond pairs "
               f"(2 S–H, 1 S–S) + {n['lone']} lone pairs"],
        source="f1.draw_molecule; cases/h2s2/v3")


def draw_dividing(fig, left, top, width, height):
    """f1.draw_dividing_point on axes `width` cm wide with the text at journal
    size: drawn at the design width WB (k = 1), then the axes are widened and
    the labels outside the frame are moved back to their design distance."""
    ax = axes(fig, left, top, f1.WB, height)
    out = f1.draw_dividing_point(ax)
    fw, fh = fig.get_size_inches() / CM
    ax.set_position([left / fw, 1 - (top + height) / fh, width / fw, height / fh])
    stretch = width / f1.WB
    for t in ax.texts:
        x = t.get_position()[0]
        if x < 0 or x > 1:
            assert t.get_transform() is not ax.transData
            t.set_x(x / stretch if x < 0 else 1 + (x - 1) / stretch)
    return ax, out


def dividing_handles():
    wf = [Line2D([], [], ls="none", ms=3.0, mew=0.0, label=f1.WF_LABEL["bond"],
                 **f1.WF_STYLE["bond"]),
          Line2D([], [], ls="none", ms=2.9, mew=0.5, label=f1.WF_LABEL["lone"],
                 **f1.WF_STYLE["lone"])]
    lines = [Line2D([], [], color=f1.MID, lw=1.0, ls=(0, (4, 2)),
                    label="Bond midpoint: nearest-atom rule"),
             Line2D([], [], color=THIS, lw=1.2,
                    label="Force-balance point $f\\,$*: this work")]
    return wf + lines


def dividing(width, height, mag):
    wb = f1.WB
    left, right, top = 0.625 * wb, 0.335 * wb, 0.46
    fig = figure(left + width + right, top + height + 2.0)
    _ax, (order, fstar, _table) = draw_dividing(fig, left, top, width, height)
    legend(fig, dividing_handles(), 0.5 * (left + width + right),
           top + height + 1.08, "upper center", ncol=2, handlelength=1.9,
           columnspacing=2.0)
    groups = f1.survey()["groups"]
    recs = [r for x in order for _, r in groups[x]]
    z = {x: next(r["zX"] for _, r in groups[x] if r["fstar"] is not None)
         for x in order}
    return fig, mag, dict(
        shows="Position f of the WF centre along every terminal X–H bond of "
              "clean38 (0 at X, 1 at H), one row per bond type, with the bond "
              "midpoint (dashed) and the force-balance point f* (solid, "
              "shaded right of it); own legend.",
        printed=["rows (Z_eff of X), n: " + "; ".join(
            f"{x}–H ({z[x]:.2f}) n = {len(groups[x])}" for x in order),
            "abscissa 0, 0.25, 0.5, 0.75, 1; pair to X: H(+1); pair to H: H(−1)"],
        drawn=["f* (position of the solid line): " + ", ".join(
            f"{x}–H {fstar[x]:.3f}" for x in order),
            f"{len(recs)} bonds; {sum(r['nearer_H'] for r in recs)} centres "
            f"nearer H than X; pair to X {sum(r['to'] == 'X' for r in recs)}, "
            f"pair to H {sum(r['to'] == 'H' for r in recs)}"],
        source="f1.draw_dividing_point; frozen plugin on the nine clean38 "
               "cases that contain H")


def panel_dividing_point():
    return dividing(4.75, f1.HB, 2.0)


def panel_dividing_point_wide():
    return dividing(9.10, 3.95, 2.4)


def panel_rule():
    fig = figure(f1.WC, f1.H)
    ex = f1.draw_rule(axes(fig, 0, 0, f1.WC, f1.H))
    rows = []
    for key, fmt, third in (("chi", "{:.2f}", "rho2"), ("force", "{:.3f}", "d")):
        e = ex[key]
        last = "chi" if key == "chi" else "force"
        rows.append(f"{name_of(e['case'])}, " + "; ".join(
            f"{end['el']} {end['zeff']:.2f}, {fmt.format(end[third])}, "
            f"{end[last]:.2f}" for end in e["ends"])
            + f" -> pair to {e['ends'][0]['el']}")
    return fig, 2.0, dict(
        shows="The decision rule: sketch of the quantities (own WFs, "
              "contested WF, ρ, d), Z_eff = Z_val − n_inner − 5/16, the two "
              "criteria χ = Z_eff/ρ² and Z_eff/d², and one worked example of each.",
        printed=["5/16",
                 "χ example (Z_eff, ρ² in Å², χ): " + rows[0],
                 "force example (Z_eff, d in Å, Z_eff/d²): " + rows[1]],
        drawn=[],
        source="f1.draw_rule; cases/as4s4/v3 and cases/h2s2/v3")


# ---------------------------------------------------- benchmark (figure 2)
BAR = 0.20 * 5.0 / f2.H_A      # bar thickness in rows, as in the journal figure


def accuracy(bar_w, height):
    """f2.draw_accuracy with bars `bar_w` cm long on a strip `height` cm high."""
    name_w, count_w, gap = 4.85, 1.28, 0.42          # of draw_accuracy, k = 1
    frame = name_w + 3 * (bar_w + count_w) + 2 * gap
    fig = figure(0.1 + frame + 0.2, 0.15 + height + 0.95)
    subs = f2.draw_accuracy(axes(fig, 0.1, 0.15, frame, height))
    for sub in subs:
        for p in sub.patches:
            mid = p.get_y() + 0.5 * p.get_height()
            p.set_height(BAR)
            p.set_y(mid - 0.5 * BAR)
    return fig, subs


def accuracy_numbers(metrics):
    bench = f2.data()["bench"]
    return ["; ".join(f"{plain(f2.MLABEL[key])} "
                      f"{bench[key]['summary'][num]}/{bench[key]['summary'][den]}"
                      for key, *_r in f2.METHODS)
            + f"  [{label}; axis 0, 50, 100]" for num, den, label in metrics]


def panel_accuracy():
    fig, _subs = accuracy(1.78, 2.30)
    return fig, 2.0, dict(
        shows="Accuracy of the five methods on clean38: atoms correct, cases "
              "fully correct, element-mean pass (bars in %, counts printed); "
              "method names with colour and marker on the left.",
        printed=accuracy_numbers(f2.METRICS), drawn=[],
        source="f2.draw_accuracy; results/clean38_*.tsv")


def panel_accuracy_atoms():
    fig, subs = accuracy(4.90, 2.70)
    subs[1].remove()
    subs[2].remove()
    return fig, 2.6, dict(
        shows="Atoms correct on clean38, five methods: the first block of the "
              "accuracy panel alone, as a large bar chart.",
        printed=accuracy_numbers(f2.METRICS[:1]), drawn=[],
        source="f2.draw_accuracy (blocks 2 and 3 removed); results/clean38_*.tsv")


def panel_cases():
    w, h = 3.60, f2.H_B
    ml, mr, mt = 1.95, 2.45, 0.50
    fig = figure(ml + w + mr, mt + h + 2.2)
    clean, rest = f2.draw_cases(axes(fig, ml, mt, w, h))
    handles = [marker(f2.STYLE[m], f2.MLABEL[m]) for m in f2.TRIO]
    handles.append(Line2D([], [], ls="none", label="of: labelled atoms"))
    legend(fig, handles, 0.15, mt + h + 1.10, "upper left", ncol=2)
    tot = {m: sum(r["wrong"][m] for r in clean + rest) for m in f2.TRIO}
    row = lambda r: (f"{name_of(r['key'][1])} "
                     + "/".join(str(r["wrong"][m]) for m in f2.TRIO)
                     + f" of {r['n']}")
    return fig, 2.0, dict(
        shows="Incorrect atoms (% of the labelled atoms) of this work, the "
              "min-cost flow and the per-WF greedy rule in the 14 cases that "
              "at least one of the three does not get fully right; counts in "
              "the columns; own legend.",
        printed=["this work / min-cost flow / greedy of labelled atoms, top to "
                 "bottom: " + "; ".join(row(r) for r in clean),
                 "below the dashed line: " + "; ".join(row(r) for r in rest),
                 "abscissa 0, 50, 100"],
        drawn=["column sums (not printed): " + ", ".join(
            f"{plain(f2.MLABEL[m])} {tot[m]}" for m in f2.TRIO)],
        source="f2.draw_cases; results/clean38_r3-T-B.tsv, clean38_flow.tsv, "
               "clean38_greedy.tsv")


def panel_elements():
    w, h = 11.4, (f2.F_ROW_Y + 1.0) * 0.76
    fig = figure(w + 0.2, h + 0.2)
    full, part = f2.draw_elements(axes(fig, 0.1, 0.1, w, h))
    els = f2.data()["analysis"]["elements"]
    by_n = sorted(els, key=lambda s: (-els[s]["labelled"], s))
    return fig, 2.1, dict(
        shows="Periodic-table map of clean38 for this work: each tile is "
              "filled in proportion to the labelled atoms of the element that "
              "are correct; the number is the count of labelled atoms; own key.",
        printed=[f"every labelled atom correct ({len(full)} elements); "
                 + ", ".join(f"{s} {els[s]['ok']}/{els[s]['labelled']}"
                             for s in part) + " correct",
                 "labelled atoms: " + ", ".join(
                     f"{s} {els[s]['labelled']}" for s in by_n)],
        drawn=[f"{len(els)} elements in clean38"],
        source="f2.draw_elements; _cache/r3tb/fig2_analysis.json")


def panel_chi():
    _fn, w, h, (ml, mr, mt, mb) = f2.PANELS["chi"]
    fig = figure(ml + w + mr, mt + h + mb)
    rows, rs_p, rs_all = f2.draw_chi(axes(fig, ml, mt, w, h))
    n_p = sum(r["block"] == "p" for r in rows)
    return fig, 2.0, dict(
        shows="In-situ χ = Z_eff/ρ² of the eight elements that carry lone "
              "pairs at a contested WF (median over atoms, bar: interquartile "
              "range) against the Allen electronegativity; own legend.",
        printed=[f"p-block elements: r_s = {rs_p:.2f} (n = {n_p}); d-block "
                 f"element; all elements: r_s = {rs_all:.2f} (n = {len(rows)})",
                 "ordinate (log) 2, 5, 10, 20; abscissa 1.5 to 4.0"],
        drawn=["median χ (e Å⁻²): " + ", ".join(
            f"{r['el']} {r['chi']:.2f}" for r in sorted(rows, key=lambda r: r["chi"]))],
        source="f2.draw_chi; _cache/r3tb/fig2_analysis.json; abscissa "
               "pyleos.capacity.ENEG_ALLEN")


def panel_constants():
    w, h = 13.0, 3.0
    ml, mt = 1.25, 0.15
    fig = figure(ml + w + 0.2, mt + h + 2.1)
    _subs, info = f2.draw_constants(axes(fig, ml, mt, w, h))
    handles = [
        Line2D([], [], color=C["muted"], lw=0.5, ls=(0, (4, 3)),
               marker=MARKER[f2.STYLE[f2.THIS]], ms=4.6, mfc=THIS, mec="white",
               mew=0.4, label="Default"),
        Patch(facecolor=C["light"], alpha=0.55, lw=0,
              label="All 38 cases as at the default")]
    legend(fig, handles, ml + 0.5 * w, mt + h + 1.42, "upper center", ncol=2,
           handlelength=2.2, columnspacing=2.6)
    sc = f2.data()["scan"]["scan"]
    n_all = sc["PAIR_SCREEN"]["points"][0]["atoms_labelled"]
    label = lambda k: f2.SCAN[k]["label"].replace("\n", " ")
    value = lambda k, v: (f"10^{round(math.log10(v))}" if f2.SCAN[k].get("log")
                          else f"{v:g}")
    text = lambda k: plain(f2.SCAN[k]["default_text"].replace("$^{", "^"))
    return fig, 2.0, dict(
        shows="Atoms correct on clean38 against each of the four constants of "
              "the method (the other three at their defaults); triangle and "
              "dashed line: default; grey band: every case unchanged; own key.",
        printed=["defaults: " + "; ".join(f"{label(k)} {text(k)}" for k in info),
                 f"ordinate: atoms correct (of {n_all}), ticks 700, 720, 740, 760"],
        drawn=["grey band (every case unchanged): " + "; ".join(
            f"{label(k)} {value(k, v['plateau'][0])} to "
            f"{value(k, v['plateau'][1])}" for k, v in info.items()),
            "atoms correct at the default: " + ", ".join(
                sorted({str(v["at_default"]) for v in info.values()}))
            + f"; {sum(len(s['points']) for s in sc.values())} runs"],
        source="f2.draw_constants; _cache/r3tb/fig2_scan.json")


def margins_numbers():
    single, shared = f2.margin_classes()
    cell = lambda n: str(n) if n else "–"
    lab = lambda c: (f"{c['winner']}–{c['loser']}" if c["single"] else
                     f"{c['winner']}{c['hosts'][0]}"
                     + (f"–{c['hosts'][-1]}" if len(c["hosts"]) > 1 else "")
                     + f"–{c['loser']}")
    out = []
    for title, rows in (("One atom takes the pair", single),
                        ("Like atoms share the pair", shared)):
        out.append(f"{title} {sum(len(c['force']) for c in rows)}/"
                   f"{sum(len(c['chi']) for c in rows)}: " + "; ".join(
                       f"{lab(c)} {cell(len(c['force']))}/{cell(len(c['chi']))}"
                       for c in rows))
    both = single + shared
    out.append(f"All {sum(len(c['force']) for c in both)}/"
               f"{sum(len(c['chi']) for c in both)}; abscissa 1, 10, 100, ∞")
    return ["columns: WFs decided by force / by χ"] + out


MARGINS_SHOWS = ("Winner/loser ratio of the deciding quantity (force Z_eff/d² "
                 "or χ = Z_eff/ρ²) of every contested WF with a heteronuclear "
                 "partner; rows are bond classes, winner first; counts by rule "
                 "in the columns; own legend.")
MARGINS_SOURCE = "f2.draw_margins; _cache/r3tb/fig2_analysis.json"


def panel_margins():
    """The two blocks of draw_margins side by side at the journal row pitch."""
    single, shared = f2.margin_classes()
    n1, n2 = len(single), len(shared)
    top, bottom = -0.85, n1 + n2 + 2 * 1.45 - 0.25 + 0.75    # of draw_margins
    pitch = f2.H_F / (bottom - top)                           # cm per row
    windows = ((top, n1 + 0.75), (n1 + 1.45 + top, bottom))
    w, ml, mr, mt, sep = 4.30, 1.75, 1.35, 0.55, 0.45
    col = ml + w + mr
    h1 = (windows[0][1] - windows[0][0]) * pitch
    fig = figure(2 * col + sep, mt + h1 + 0.95)
    for i, (lo, hi) in enumerate(windows):
        h = (hi - lo) * pitch
        ax = axes(fig, i * (col + sep) + ml, mt, w, h)
        f2.draw_margins(ax)
        assert ax.get_ylim() == (bottom, top), ax.get_ylim()
        ax.set_ylim(hi, lo)
        for t in ax.texts:
            t.set_visible(lo < t.get_position()[1] < hi)
        if i:
            ax.get_legend().remove()
            fw, fh = fig.get_size_inches() / CM
            fig.text((i * (col + sep) + 0.30) / fw, 1 - (mt + h + 1.55) / fh,
                     "Rows: bond classes, winner–loser\n"
                     "Columns: WFs decided by force, by $\\chi$\n"
                     "Grey slot (∞): $Z_{\\mathrm{eff}}$ of the loser ≤ 0\n"
                     "Dashed line: ratio 1",
                     ha="left", va="top", linespacing=1.55)
    return fig, 1.85, dict(shows=MARGINS_SHOWS + " Two columns: left, one "
                           "atom takes the pair; right, like atoms share the "
                           "pair, and the totals.",
                           printed=margins_numbers(), drawn=[],
                           source=MARGINS_SOURCE)


def panel_margins_column():
    _fn, w, h, (ml, mr, mt, mb) = f2.PANELS["margins"]
    fig = figure(ml + w + mr, mt + h + mb)
    f2.draw_margins(axes(fig, ml, mt, w, h))
    return fig, 1.72, dict(shows=MARGINS_SHOWS + " Journal layout in one "
                           "column; it fills the height of the slide (no "
                           "room for a title band).",
                           printed=margins_numbers(), drawn=[],
                           source=MARGINS_SOURCE)


def swarm(x, pitch):
    """Offsets from the centre line that keep points with the abscissae x at
    least `pitch` apart (x, pitch and the offsets in one unit): in order of
    x, each point takes the free place nearest to the centre line."""
    x = np.asarray(x, float)
    y = np.zeros(len(x))
    placed = []
    for i in np.argsort(x, kind="stable"):
        near = [j for j in placed if abs(x[j] - x[i]) < pitch]
        places = [0.0]
        for j in near:
            dy = math.sqrt(pitch ** 2 - (x[j] - x[i]) ** 2)
            places += [y[j] + dy, y[j] - dy]
        y[i] = next(p for p in sorted(places, key=lambda p: (abs(p), p))
                    if all((x[j] - x[i]) ** 2 + (y[j] - p) ** 2
                           >= pitch ** 2 * (1 - 1e-9) for j in near))
        placed.append(i)
    return y


def panel_margins_summary():
    """The ratios of draw_margins in two rows, one per deciding rule; the
    axis, the dashed line and the grey slot are those of draw_margins (an
    infinite ratio is drawn at INF_X).  Every pair is one marker; markers
    with the same ratio form a stack about the row, and the stacks of a row
    are compressed by one factor so that the highest fills the row: the
    height of a stack is proportional to the number of pairs in it.
    Ordinate in cm at journal scale."""
    wfs = f2.data()["analysis"]["wfs"]
    mk = MARKER[f2.STYLE[f2.THIS]]
    # rule, label, marker, row height.  Markers as in draw_margins, with
    # thinner edges (there 0.25 and 0.7 pt) so that a stack shows its markers.
    rows = (
        ("force", r"Force $Z_{\mathrm{eff}}$/$d^{2}$",
         dict(ls="none", marker="o", ms=2.7, mfc=THIS, mec="white", mew=0.15),
         1.36),
        ("chi", r"$\chi$ = $Z_{\mathrm{eff}}$/$\rho^{2}$",
         dict(ls="none", marker=mk, ms=4.3, mfc="white", mec=C["ink"], mew=0.5),
         1.54))
    w, ml, mt = 8.74, 0.20, 0.10
    head, foot, gap = 0.44, 0.36, 0.04      # label line, median line, between
    h = 0.16 + sum(head + r[3] + foot for r in rows) + gap * (len(rows) - 1)
    lo, hi = 0.8, f2.INF_HI                                 # as draw_margins
    cm = lambda v: w * math.log(v / lo) / math.log(hi / lo)       # ratio -> cm
    ratio = lambda u: lo * (hi / lo) ** (u / w)                   # cm -> ratio
    fig = figure(ml + w + 0.20, mt + h + 0.95)
    ax = axes(fig, ml, mt, w, h)
    ax.set_xscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(h, 0)
    ink, top, out = C["ink"], 0.06, {}
    for rule, label, kw, height in rows:
        vals = sorted(math.inf if x["margin"] == "inf" else x["margin"]
                      for x in wfs if x["rule"] == rule)
        fin = [v for v in vals if not math.isinf(v)]
        x = np.array([f2.INF_X if math.isinf(v) else v for v in vals])
        y0, y1 = top + head, top + head + height
        size = (kw["ms"] + kw["mew"]) / 72 * 2.54             # marker, cm
        ax.add_patch(Rectangle((f2.INF_LO, y0), f2.INF_HI - f2.INF_LO, height,
                               fc=C["light"], ec="none", alpha=0.55, zorder=0))
        ax.plot([1.0, 1.0], [y0, y1], color=ink, lw=0.5, ls=(0, (4, 2.5)),
                zorder=0)
        off = swarm([cm(v) for v in x], size)
        squeeze = min(1.0, 0.5 * (height - size) / np.abs(off).max())
        order = np.argsort(off, kind="stable")      # drawn from the top down
        ax.plot(x[order], 0.5 * (y0 + y1) + squeeze * off[order], zorder=3, **kw)
        med = float(np.median(vals))
        ax.plot([med, med], [y1 + 0.05, y1 + foot - 0.03], color=ink, lw=1.0,
                solid_capstyle="butt", zorder=5)
        ax.text(ratio(cm(med) + 0.10), y1 + 0.5 * foot + 0.01,
                f"median {med:.2f}", ha="left", va="center")
        ax.text(ratio(0.10), top + 0.5 * head - 0.03,
                f"{label} — {len(vals)} pairs", ha="left", va="center")
        out[rule] = dict(n=len(vals), n_inf=len(vals) - len(fin), med=med,
                         lo=fin[0], hi=fin[-1], squeeze=squeeze)
        top = y1 + foot + gap
    ax.set_xticks([1, 10, 100, f2.INF_X])
    ax.set_xticklabels(["1", "10", "100", "∞"])
    ax.xaxis.set_minor_locator(FixedLocator(
        [m * 10 ** e for e in (0, 1) for m in range(2, 10)]))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_yticks([])
    ax.set_xlabel("Ratio winner / loser")
    f, c = out["force"], out["chi"]
    return fig, 2.0, dict(
        shows="Winner/loser ratio of the deciding quantity of every bond pair "
              "with ends of two elements, in two rows: decided by the force "
              "Z_eff/d² (filled circles) and by χ = Z_eff/ρ² (open "
              "triangles); one marker per pair, pairs with the same ratio "
              "stacked about the row (stack height proportional to the number "
              "of pairs); tick: median of the row; dashed line: ratio 1; grey "
              "slot: ratio ∞.",
        printed=[f"Force Z_eff/d² — {f['n']} pairs; median {f['med']:.2f}",
                 f"χ = Z_eff/ρ² — {c['n']} pairs; median {c['med']:.2f}",
                 "abscissa 1, 10, 100, ∞"],
        drawn=[f"{f['n'] + c['n']} pairs in all; {f['n_inf'] + c['n_inf']} in "
               f"the grey slot (ratio ∞: {f['n_inf']} force, {c['n_inf']} χ)",
               f"finite ratios: force {f['lo']:.3f} to {f['hi']:.3f}, χ "
               f"{c['lo']:.3f} to {c['hi']:.3f}; the median of the force row "
               "counts the pairs at ∞",
               "markers of a stack overlap: the step is "
               f"{100 * f['squeeze']:.0f} % (force) and "
               f"{100 * c['squeeze']:.0f} % (χ) of the marker size"],
        source="drawn in r3tb_slide_panels.py from "
               "_cache/r3tb/fig2_analysis.json (block wfs, the records of "
               "f2.draw_margins)")


# -------------------------------------------------- case studies (figure 3)
def panel_delithiation():
    w, h, ml, mt = 6.60, 5.80, 1.80, 0.60
    fig = figure(ml + w + 0.15, mt + h + 0.95)
    f3.draw_delithiation(axes(fig, ml, mt, w, h))
    legend(fig, case_study_handles("this", "greedy"), ml, mt - 0.13,
           "lower left", ncol=2, columnspacing=2.0)
    rows = f3.delithiation_rows()
    count = lambda d: ", ".join(f"{signed(v)}: {n}" for v, n in d.items())
    return fig, 2.0, dict(
        shows="Delithiation of Li2O2 (Li16O16 cell, 0, 1, 2 Li removed): mean "
              "oxidation state of O and of Li against charge counting, and "
              "the oxidation state of every O atom (number of atoms at each "
              "value); this work and the per-WF greedy rule; own legend.",
        printed=["mean O (ordinate ticks): " + ", ".join(
            f3.minus(f"{r['ref']['o']:.4f}".rstrip("0").rstrip(".")) for r in rows)
            + "; mean Li +1 (ticks 0, +1); Li removed 0, 1, 2",
            "O atoms, this work (value: atoms): " + " | ".join(
                count(r["this"]["o_atoms"]) for r in rows),
            "O atoms, greedy (value: atoms): " + " | ".join(
                count(r["greedy"]["o_atoms"]) for r in rows)],
        drawn=["mean O, this work = charge counting = greedy: " + ", ".join(
            f3.minus(f"{r['this']['o']:.4f}") for r in rows)],
        source="f3.draw_delithiation; cases/{li16o16,li15o16,li14o16}/v31 via "
               "_cache/r3tb/applications.json")


def panel_vacancy():
    w, h, ml, mt = 7.80, 3.20, 1.05, 0.60
    fig = figure(ml + w + 0.1, mt + h + 0.95)
    f3.draw_vacancy(axes(fig, ml, mt, w, h))
    legend(fig, case_study_handles("this"), ml, mt - 0.13, "lower left")
    d = f3.vacancy_data()
    ce = lambda name: " ".join(signed(v) for v in d[name]["this"])
    n3 = lambda name: d[name]["this"].count(3) or "no"
    return fig, 2.0, dict(
        shows="O vacancy in CeO2: oxidation state of each of the four Ce atoms "
              "in Ce4O8 and Ce4O7 (this work), and the four Ce atoms around "
              "the vacancy V_O; own legend.",
        printed=["; ".join(f"{lab}: Ce1–Ce4 {ce(name)} ({n3(name)} Ce3+)"
                           for name, lab in f3.CEO2),
                 f"4 × {d['d_ce'][0]:.2f} Å; ordinate +3, +4"],
        drawn=[],
        source="f3.draw_vacancy; cases/{ce4o8,ce4o7}/v3 via "
               "_cache/r3tb/applications.json")


def panel_heldout():
    w, h, ml, mt = 11.0, 4.0, 1.05, 1.30
    fig = figure(ml + w + 0.2, mt + h + 1.05)
    f3.draw_heldout(axes(fig, ml, mt, w, h))
    legend(fig, case_study_handles("this", "flow", "greedy"), ml, mt - 0.09,
           "lower left", labelspacing=0.35)
    rows, tot = f3.heldout_rows()
    return fig, 2.0, dict(
        shows="The 13 held-out compounds with homonuclear bonds: atoms that "
              "receive the reference oxidation state (%), per compound, for "
              "this work, the min-cost flow and the per-WF greedy rule; own "
              "legends.",
        printed=["Atoms correct, 13 compounds: " + ", ".join(
            f"{lab} {tot[k][0]}/{tot[k][1]}" for k, lab in
            (("this", "this work"), ("flow", "min-cost flow"),
             ("greedy", "greedy"))), "ordinate 0, 50, 100"],
        drawn=["atoms correct of labelled, this work / flow / greedy: "
               + "; ".join(f"{r['label']} {r['this']}/{r['flow']}/{r['greedy']} "
                           f"of {r['n']}" for r in rows),
               "compounds fully correct: " + ", ".join(
                   f"{lab} {tot[k][2]}/{len(rows)}" for k, lab in
                   (("this", "this work"), ("flow", "min-cost flow"),
                    ("greedy", "greedy")))],
        source="f3.draw_heldout; _cache/r3tb/applications.json "
               "(results/testset_results_cases.tsv)")


# -------------------------------------------------------------------- main
PANELS = (
    ("Method (r3tb_fig1_method.py)",
     (("molecule", panel_molecule), ("dividing_point", panel_dividing_point),
      ("rule", panel_rule))),
    ("Benchmark on clean38 (r3tb_fig2_benchmark.py)",
     (("accuracy", panel_accuracy), ("cases", panel_cases),
      ("elements", panel_elements), ("chi", panel_chi),
      ("constants", panel_constants), ("margins", panel_margins))),
    ("Case studies (r3tb_fig3_case_studies.py)",
     (("delithiation", panel_delithiation), ("vacancy", panel_vacancy),
      ("heldout", panel_heldout))),
    ("Composites and variants",
     (("accuracy_atoms", panel_accuracy_atoms),
      ("dividing_point_wide", panel_dividing_point_wide),
      ("margins_column", panel_margins_column),
      ("margins_summary", panel_margins_summary))),
)


def write_index(done):
    key = [(f"{k}", C[k], MARKER[k], lab) for k, lab in (
        ("aux1", "In-situ electronegativity (this work)"),
        ("flow", "Min-cost flow"), ("nearest", "Per-WF greedy rule"),
        ("geo", "Min-cost flow, β = 0"), ("aux2", "Lewis-first with Allen table"))]
    L = ["# Slide panels of the r3-T-B figures",
         "",
         "Written by `python3 r3tb_slide_panels.py` (run in `paper/`; it "
         "regenerates every image and this file).",
         "",
         f"PNG, {DPI} dpi, white background, tight bounding box, no panel "
         f"letter. Natural size = pixels / {DPI}: place at 100 % on a "
         f"{SLIDE[0]} x {SLIDE[1]} in slide. Each panel is the journal panel "
         "magnified M times (fonts, markers and line weights together); "
         "\"text\" is the smallest and largest font size at natural size.",
         "",
         "Method colours and markers (paper/style.py): " + "; ".join(
             f"{lab} {col} `{mk}`" for _k, col, mk, lab in key) + ".",
         "",
         "Added for the slides, not in the journal panels: the legends of "
         "`dividing_point`, `cases`, `constants`, `delithiation`, `vacancy` "
         "and `heldout` (method names), and the four lines of notes in "
         "`margins`. `margins` shows the two blocks of the journal panel side "
         "by side; `margins_column` is the journal layout; `margins_summary` "
         "shows the same ratios in two rows, one per deciding rule (not a "
         "journal panel).",
         "",
         "| file | natural size (in) | aspect (w : h) | pixels | M | text (pt) |",
         "|---|---|---|---|---|---|"]
    for _group, items in done:
        for info, _meta in items:
            L.append(f"| `{info['file']}` | {info['w']:.2f} x {info['h']:.2f} | "
                     f"{info['w'] / info['h']:.2f} : 1 | {info['px'][0]} x "
                     f"{info['px'][1]} | {info['mag']:g} | {info['pt'][0]:g}–"
                     f"{info['pt'][1]:g} |")
    for group, items in done:
        L += ["", f"## {group}"]
        for info, meta in items:
            L += ["", f"### {info['file']}", "",
                  f"- Size: {info['w']:.2f} x {info['h']:.2f} in, aspect "
                  f"{info['w'] / info['h']:.2f} : 1.",
                  f"- Shows: {meta['shows']}",
                  "- Numbers printed in the image:"]
            L += [f"  - {s}" for s in meta["printed"]]
            if meta["drawn"]:
                L.append("- Drawn but not printed as text:")
                L += [f"  - {s}" for s in meta["drawn"]]
            L.append(f"- Source: {meta['source']}.")
    (OUT / "INDEX.md").write_text("\n".join(L) + "\n")


def main():
    only = sys.argv[1:]
    known = [name for _g, items in PANELS for name, _fn in items]
    if set(only) - set(known):
        sys.exit(f"unknown panel; choose from: {', '.join(known)}")
    check_caches()
    done = []
    for group, items in PANELS:
        rows = []
        for name, fn in items:
            if only and name not in only:
                continue
            fig, mag, meta = fn()
            info = save(fig, name, mag)
            rows.append((info, meta))
            print(f"  {info['file']:26s} {info['w']:5.2f} x {info['h']:4.2f} in  "
                  f"M = {mag:<4g} text {info['pt'][0]:g}-{info['pt'][1]:g} pt")
        done.append((group, rows))
    if only:
        print(f"wrote {sum(len(r) for _g, r in done)} images in {OUT}; "
              f"INDEX.md not rewritten")
        return
    write_index(done)
    print(f"wrote {sum(len(r) for _g, r in done)} images and INDEX.md in {OUT}")


if __name__ == "__main__":
    main()
