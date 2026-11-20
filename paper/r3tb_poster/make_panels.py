#!/usr/bin/env python3
"""Data panels of the r3-T-B poster (900 x 1200 mm), one SVG each, at final size.

Run:  python3 make_panels.py          (from any directory; a few seconds)
Writes panels/dividing.svg, panels/chi.svg and panels/sizes.json (final width
and height, mm).

Both panels are drawn by the panel functions of the journal figure scripts in
the parent directory (r3tb_fig1_method.py, r3tb_fig2_benchmark.py), which read
their numbers from the frozen plugin and from results/ (see the caption
files).  Nothing is typed in here.

Scale.  A panel is laid out in journal centimetres with the journal fonts
(7.5 pt text) and written as an SVG whose size attributes are multiplied by
S = 3.2: 1 cm of design is 32 mm on the poster and 7.5 pt becomes 24 pt.
Text stays text (style.py: svg.fonttype none; Arial).  The tick and legend
sizes of the journal style (7 pt) are raised to 7.5 pt so that no label is
below 24 pt; the panel of fig. 1 scales with its axes (dividing: k = 1.07, so
7 pt x k = 7.5 pt).

Poster-specific choice: the second legend entry of chi reads "Mo (d block)".
The schematics and the accuracy chart are hand-laid: make_schematics.py.
"""
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
sys.path.insert(0, str(PAPER))

import matplotlib as mpl                              # noqa: E402
import style                                          # noqa: E402
from style import CM                                  # noqa: E402
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.transforms import Bbox                # noqa: E402

import r3tb_fig1_method as f1                         # noqa: E402
import r3tb_fig2_benchmark as f2                      # noqa: E402

S = 3.2                       # poster mm per design mm (7.5 pt -> 24 pt)
OUT = HERE / "panels"
SIZES = {}

# journal style with the 7 pt labels raised to the 7.5 pt of the body text
POSTER_RC = {"xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
             "legend.fontsize": 7.5, "svg.hashsalt": "r3tb-poster"}


def figure(w_cm, h_cm):
    return plt.figure(figsize=(w_cm * CM, h_cm * CM))


def axes(fig, left, top, width, height):
    """Axes placed in cm from the left / top edges of the figure."""
    fw, fh = fig.get_size_inches() / CM
    return fig.add_axes([left / fw, 1 - (top + height) / fh,
                         width / fw, height / fh])


def save(fig, name, crop=None):
    """Write panels/<name>.svg at final size.  crop = (x0, y0, x1, y1) in cm
    from the bottom-left corner of the figure."""
    OUT.mkdir(exist_ok=True)
    path = OUT / f"{name}.svg"
    kw = {}
    if crop is not None:
        kw["bbox_inches"] = Bbox([[crop[0] * CM, crop[1] * CM],
                                  [crop[2] * CM, crop[3] * CM]])
    fig.savefig(path, format="svg", metadata={"Date": None}, **kw)
    plt.close(fig)
    svg = path.read_text()
    m = re.search(r'<svg[^>]*?width="([\d.]+)pt" height="([\d.]+)pt"', svg)
    w_mm = float(m.group(1)) / 72 * 25.4 * S
    h_mm = float(m.group(2)) / 72 * 25.4 * S
    svg = svg.replace(f'width="{m.group(1)}pt" height="{m.group(2)}pt"',
                      f'width="{w_mm:.2f}mm" height="{h_mm:.2f}mm"', 1)
    path.write_text(svg)
    SIZES[name] = {"w_mm": round(w_mm, 1), "h_mm": round(h_mm, 1)}
    print(f"  {name:13s} {w_mm:6.1f} x {h_mm:6.1f} mm")


# -------------------------------------------------------------------- fig. 1
def panel_dividing(k=1.07, hb=3.25):
    """Strip plot of f; axes k * WB wide, labels outside (see the docstring
    of draw_dividing_point)."""
    wb = f1.WB * k
    left, right, top, bottom = 0.625 * wb, 0.335 * wb, 0.46 * k, 1.02 * k
    fig = figure(left + wb + right, top + hb + bottom)
    ax = axes(fig, left, top, wb, hb)
    f1.draw_dividing_point(ax)
    save(fig, "dividing")



# -------------------------------------------------------------------- fig. 2
def panel_chi(w=5.29, h=4.3):
    ml, mr, mt, mb = 1.45, 0.20, 0.12, 0.90
    fig = figure(ml + w + mr, mt + h + mb)
    ax = axes(fig, ml, mt, w, h)
    f2.draw_chi(ax)
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi * 1.38)          # head room: the legend clears F
    # poster legend: the p-block entry as drawn, the open symbol named only
    ax.get_legend().get_texts()[1].set_text("Mo (d block)")
    ax.yaxis.set_label_coords(-0.92 / w, 0.5)
    save(fig, "chi")



def main():
    print(f"panels at {S} x journal scale -> {OUT}")
    with mpl.rc_context(POSTER_RC):
        panel_dividing()
        panel_chi()
    (OUT / "sizes.json").write_text(json.dumps(SIZES, indent=1) + "\n")


if __name__ == "__main__":
    main()
