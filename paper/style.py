"""Shared style for the LEOS paper figures (journal style, English labels).

Every figure script in this directory imports this module so that the five
figures read as one set: same fonts, sizes, line weights, colours and
method-to-colour mapping.

Sizes follow the RSC template: single column 8.3 cm, double column 17.1 cm.
The categorical palette is Okabe-Ito; it was checked with the dataviz
validator (all checks pass; the purple/green pair is only 7.6 apart for
deuteranopes, so every series also carries its own marker shape).
"""
import re
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parent

CM = 1 / 2.54
SINGLE = 8.3 * CM      # single-column width, inches
DOUBLE = 17.1 * CM     # double-column width, inches

# One colour and one marker per method, fixed across all figures.
C = {
    "flow": "#0072B2",      # min-cost flow (this work)
    "nearest": "#D55E00",   # per-WF greedy baseline (EOSCAR OS(int) column); not a pure nearest-atom rule
    "geo": "#009E73",       # geometric flow, beta = 0
    "aux1": "#CC79A7",
    "aux2": "#E69F00",
    "aux3": "#56B4E9",
    "ref": "#4D4D4D",       # reference values, charge counting
    "ink": "#1A1A1A",       # text, axes
    "muted": "#7F7F7F",     # secondary text, guide lines
    "light": "#D9D9D9",     # fills that must recede
}
MARKER = {"flow": "o", "nearest": "s", "geo": "^", "ref": "D",
          "aux1": "v", "aux2": "P", "aux3": "X"}
LABEL = {
    "flow": "Min-cost flow (this work)",
    "nearest": "Per-WF greedy rule",
    "geo": r"Min-cost flow, $\beta$ = 0",
    "ref": "Reference",
}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "font.size": 7.5,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.6,
    "axes.edgecolor": C["ink"],
    "axes.labelcolor": C["ink"],
    "text.color": C["ink"],
    "xtick.color": C["ink"],
    "ytick.color": C["ink"],
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.size": 1.5,
    "ytick.minor.size": 1.5,
    "xtick.minor.width": 0.5,
    "ytick.minor.width": 0.5,
    "lines.linewidth": 1.0,
    "lines.markersize": 3.5,
    "patch.linewidth": 0.5,
    "hatch.linewidth": 0.5,
    "legend.frameon": False,
    "legend.handlelength": 1.4,
    "legend.handletextpad": 0.5,
    "legend.borderaxespad": 0.4,
    "axes.grid": False,
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})


def formula(s):
    """'As4S4' -> 'As$_4$S$_4$'; 'SrBe(OH)4' -> 'SrBe(OH)$_4$'."""
    return re.sub(r"(?<=[A-Za-z\)])(\d+)", r"$_{\1}$", s)


def panel_label(ax, letter, dx=-0.02, dy=1.03):
    """Bold '(a)' above the top-left corner of an axes."""
    ax.text(dx, dy, f"({letter})", transform=ax.transAxes, fontsize=9,
            fontweight="bold", ha="right", va="bottom")


def save(fig, name):
    """Write <name>.pdf (vector) and <name>.png (300 dpi preview) here."""
    for ext, kw in (("pdf", {}), ("png", {"dpi": 300})):
        fig.savefig(HERE / f"{name}.{ext}", bbox_inches="tight",
                    pad_inches=0.02, **kw)
    plt.close(fig)
