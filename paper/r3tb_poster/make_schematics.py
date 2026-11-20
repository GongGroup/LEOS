#!/usr/bin/env python3
"""Hand-laid schematics of the r3-T-B poster (1 SVG unit = 1 mm on the poster).

Run:  python3 make_schematics.py
Writes into panels/:
  basis.svg      atom-centred basis, plane waves, Wannier functions (360 mm wide)
  pipeline.svg   from DFT to oxidation states, four steps (360 mm)
  rules.svg      the assignment rule (846 mm)
  accuracy.svg   atoms assigned correctly by five methods (340 mm)
  gpuw90.svg     the Wannier step as batched tensor algebra (846 mm)
  qr.svg         QR code of the repository link (assets/qr_leos.txt, written by
                 make_qr.swift: swift make_qr.swift <url> Q | sed 1d > assets/qr_leos.txt)

Nothing is typed in: the molecule, the classes of its Wannier functions (WFs),
the owners and the two worked examples come from r3tb_fig1_method (frozen
plugin on cases/h2s2/v3 and cases/as4s4/v3); the accuracies are summed from
results/clean38_*.tsv (see ACCURACY for the plugins that produce them); the timings come from fig3_gpuw90.load_iteration_times
(results/w90_threads_scan.tsv and the V100 table of gpuw90/README.md) and,
for gpuw90 on a CPU with a known thread count, results/gpuw90_cpu_bench_local.tsv;
the agreement with wannier90 is checked against the lines of gpuw90/README.md that
record it.  Every number drawn is printed with its source.
"""
import csv
import math
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
sys.path.insert(0, str(PAPER))

import r3tb_fig1_method as f1                         # noqa: E402
import fig3_gpuw90 as g3                              # noqa: E402

OUT = HERE / "panels"
FONT = "Avenir Next, Helvetica Neue, Arial, sans-serif"
RED, INK, GREY = "#B5223A", "#1B1B1F", "#66666C"
SOFT, CARD, LINE = "#A5A29E", "#F4F2F0", "#DDD9D6"
TINT, BLOB, CHEV = "#F6E3E6", "#F1CDD3", "#D9A3AC"
THIS = "#CC79A7"                 # this work in the data panels (paper/style.py: aux1)
BODY, SMALL = 10.6, 8.6          # 30 pt and 24.4 pt


# ------------------------------------------------------------------ SVG text
def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def spans(markup, size):
    """'*Z*_{eff} / *d*^{2}' -> tspans.  *..* italic, _{..} sub, ^{..} sup."""
    out, shift, italic = [], 0.0, False
    for tok in re.findall(r"_\{[^}]*\}|\^\{[^}]*\}|\*|[^*_^]+|[_^]", markup):
        if tok == "*":
            italic = not italic
            continue
        level, body = 0.0, tok
        if tok.startswith("_{"):
            level, body = 0.26 * size, tok[2:-1]
        elif tok.startswith("^{"):
            level, body = -0.40 * size, tok[2:-1]
        fs = size if level == 0.0 else 0.68 * size
        parts, it = [], italic
        for piece in re.split(r"(\*)", body):
            if piece == "*":
                it = not it
            elif piece:
                parts.append((piece, it))
        for piece, it in parts:
            attr = f' font-size="{fs:.2f}"'
            if it:
                attr += ' font-style="italic"'
            if abs(level - shift) > 1e-9:
                attr += f' dy="{level - shift:.2f}"'
                shift = level
            out.append(f"<tspan{attr}>{esc(piece)}</tspan>")
    return "".join(out)


def text(x, y, markup, size=BODY, anchor="start", weight=400, fill=INK):
    return (f'<text x="{x:.2f}" y="{y:.2f}" font-size="{size:.2f}" '
            f'font-weight="{weight}" text-anchor="{anchor}" fill="{fill}">'
            f"{spans(markup, size)}</text>")


def duo(x, y, lead, rule, size=BODY):
    """Condition in ink, rule in red, one line."""
    return (f'<text x="{x:.2f}" y="{y:.2f}" font-size="{size:.2f}" fill="{INK}" '
            f'font-weight="600">{spans(lead, size)}'
            f'<tspan fill="{RED}" font-weight="700">{spans(rule, size)}</tspan></text>')


# ---------------------------------------------------------------- SVG shapes
def circle(x, y, r, fill="#fff", stroke=INK, sw=1.1, extra=""):
    return (f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r:.2f}" fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{sw}"{extra}/>')


def line(x1, y1, x2, y2, stroke=INK, sw=1.0, extra=""):
    return (f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
            f'stroke="{stroke}" stroke-width="{sw}"{extra}/>')


def block_arrow(x, y, length, fill, t=5.0, hw=12.5, hl=10.0):
    """Horizontal arrow from x over `length` mm (negative: to the left)."""
    s = 1.0 if length >= 0 else -1.0
    L = abs(length)
    hl = min(hl, 0.7 * L)
    pts = [(0, -t / 2), (L - hl, -t / 2), (L - hl, -hw / 2), (L, 0),
           (L - hl, hw / 2), (L - hl, t / 2), (0, t / 2)]
    d = " ".join(f"{x + s * px:.2f},{y + py:.2f}" for px, py in pts)
    return f'<polygon points="{d}" fill="{fill}"/>'


def dimension(x1, x2, y, label, size=SMALL):
    return (line(x1, y, x2, y, GREY, 0.7) + line(x1, y - 2.6, x1, y + 2.6, GREY, 0.7)
            + line(x2, y - 2.6, x2, y + 2.6, GREY, 0.7)
            + text(0.5 * (x1 + x2), y + 0.98 * size, label, size, "middle", 400, GREY))


def atom(x, y, r, label, winner=False, size=None):
    size = size or (0.92 * r if len(label) == 1 else 0.74 * r)
    return (circle(x, y, r, TINT if winner else "#fff", RED if winner else INK,
                   1.5 if winner else 1.1)
            + text(x, y + 0.36 * size, label, size, "middle", 700,
                   RED if winner else INK))


def chevron(x, y, h=22.0):
    return (f'<polyline points="{x:.2f},{y - h:.2f} {x + 0.5 * h:.2f},{y:.2f} '
            f'{x:.2f},{y + h:.2f}" fill="none" stroke="{CHEV}" stroke-width="3.4" '
            f'stroke-linecap="round" stroke-linejoin="round"/>')


def svg(width, height, body):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}mm" '
            f'height="{height}mm" viewBox="0 0 {width} {height}" '
            f'font-family="{FONT}">\n' + "\n".join(body) + "\n</svg>\n")


PROV = []


def prov(what, value, source):
    PROV.append((what, value, source))


# ------------------------------------------------------------------ pipeline
def chip(x, y, w, h, label, kind="plain", size=9.2):
    """Rounded box with a centred label: plain (ink outline), soft (grey
    outline) or solid (red, white text)."""
    fill, stroke, ink, weight = {"plain": ("#fff", INK, INK, 600),
                                 "soft": ("#fff", SOFT, GREY, 600),
                                 "solid": (RED, RED, "#fff", 700)}[kind]
    return (f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" rx="{h / 2:.2f}" '
            f'fill="{fill}" stroke="{stroke}" stroke-width="1.2"/>'
            + text(x + w / 2, y + h / 2 + 0.35 * size, label, size, "middle", weight, ink))


def thin_arrow(x1, y1, x2, y2, colour=GREY, sw=1.3, head=4.2):
    """Line with a filled head at (x2, y2); horizontal or vertical."""
    dx, dy = x2 - x1, y2 - y1
    L = math.hypot(dx, dy)
    ux, uy = dx / L, dy / L
    bx, by = x2 - head * 1.5 * ux, y2 - head * 1.5 * uy
    pts = [(x2, y2), (bx - head * 0.62 * uy, by + head * 0.62 * ux),
           (bx + head * 0.62 * uy, by - head * 0.62 * ux)]
    d = " ".join(f"{px:.2f},{py:.2f}" for px, py in pts)
    return line(x1, y1, bx, by, colour, sw) + f'<polygon points="{d}" fill="{colour}"/>'


def badge(x, y, n, r=5.6):
    return (circle(x, y, r, RED, RED, 0)
            + text(x, y + 2.75, str(n), 7.8, "middle", 700, "#fff"))


def draw_pipeline(W=360.0, H=51.0):
    """DFT -> Wannier functions -> owner of each pair -> oxidation states."""
    b = []
    y, h = 1.0, 24.0
    nodes = [("DFT", 40.0, "soft"), ("Wannier functions", 93.0, "plain"),
             ("owner of each pair", 100.0, "solid"), ("oxidation states", 88.0, "plain")]
    gap = (W - 1.6 - sum(n[1] for n in nodes)) / 3     # 0.8 mm for the outlines
    x, mids = 0.8, []
    for i, (lab, w, kind) in enumerate(nodes):
        b.append(chip(x, y, w, h, lab, kind))
        mids.append(x + w / 2)
        if i < 3:
            b.append(thin_arrow(x + w + 2.2, y + h / 2, x + w + gap - 2.2, y + h / 2))
        x += w + gap
    yy = y + h + 15.5
    b.append(text(mids[0], yy, "VASP", SMALL, "middle", 400, GREY))
    b.append(badge(mids[1] - 34.5, yy - 3.0, 5))
    b.append(text(mids[1] - 26.5, yy, "gpuw90, on GPU", SMALL, "start", 400, INK))
    b.append(badge(mids[2] - 28.5, yy - 3.0, 1))
    b.append(text(mids[2] - 20.5, yy, "the new rule", SMALL, "start", 700, RED))
    b.append(text(mids[3], yy, "one integer per atom", SMALL, "middle", 400, GREY))
    return svg(W, H, b)


def wave(x0, x1, y, amp, lam, phase, colour, sw=1.3):
    pts = []
    n = 80
    for i in range(n + 1):
        x = x0 + (x1 - x0) * i / n
        pts.append(f"{x:.2f},{y + amp * math.sin(2 * math.pi * (x - x0) / lam + phase):.2f}")
    return (f'<polyline points="{" ".join(pts)}" fill="none" stroke="{colour}" '
            f'stroke-width="{sw}" stroke-linecap="round"/>')


def draw_basis(W=360.0, H=73.0):
    """Atom-centred basis functions, plane waves, Wannier functions."""
    b = []
    tw, gap, left = 111.5, 12.05, 0.7              # the outlines stay inside
    tiles = [("GTO / LCAO basis", INK, "each function has an atom",
              "electrons can be projected"),
             ("Plane-wave basis", INK, "no function has an atom",
              "nothing to project on"),
             ("Wannier functions", RED, "pairs with a position",
              "electrons can be counted")]
    for i, (title, colour, c1, c2) in enumerate(tiles):
        x0 = left + i * (tw + gap)
        ya = 31.0
        xa, xb = x0 + 30.0, x0 + 82.0
        b.append(text(x0, 9.6, title, 9.4, "start", 700, colour))
        b.append(f'<rect x="{x0:.2f}" y="15" width="{tw}" height="32" rx="4" fill="#fff" '
                 f'stroke="{RED if i == 2 else LINE}" stroke-width="{1.2 if i == 2 else 0.8}"/>')
        if i == 0:
            for x in (xa, xb):
                b.append(circle(x, ya, 13.5, BLOB, "none", 0, ' opacity="0.6"'))
                b.append(circle(x, ya, 6.8, "#E3A6B0", "none", 0, ' opacity="0.85"'))
        elif i == 1:
            b.append(wave(x0 + 5, x0 + tw - 5, ya - 7.5, 3.6, 20.0, 0.4, "#E3A6B0"))
            b.append(wave(x0 + 5, x0 + tw - 5, ya, 4.4, 34.0, 1.9, SOFT))
            b.append(wave(x0 + 5, x0 + tw - 5, ya + 7.5, 3.6, 51.0, 0.0, "#E3A6B0"))
        else:
            b.append(line(xa, ya, xb, ya, SOFT, 1.0))
            for x, y in ((xa - 13.5, ya - 4.5), (xb + 12.0, ya - 7.2), (xb + 12.0, ya + 7.2)):
                b.append(circle(x, y, 6.6, BLOB, "none", 0, ' opacity="0.75"'))
                b.append(circle(x, y, 3.0, "#fff", INK, 1.0))
            xm = 0.5 * (xa + xb) + 5.0
            b.append(circle(xm, ya, 8.0, BLOB, "none", 0, ' opacity="0.75"'))
            b.append(circle(xm, ya, 3.4, INK, INK, 0))
        for x in (xa, xb):
            b.append(circle(x, ya, 2.8, INK, INK, 0))
        b.append(text(x0, 58.2, c1, SMALL, "start", 600, INK))
        b.append(text(x0, 70.0, c2, SMALL, "start", 400, GREY))
    b.append(block_arrow(left + tw + gap + tw + 1.2, 31.0, gap - 2.4, RED, t=3.0, hw=8.0, hl=5.5))
    return svg(W, H, b)


# --------------------------------------------------------------------- rules
def molecule_xy(mol, cx, cy, scale):
    pts = list(mol["qa"]) + list(mol["qw"])
    x0 = 0.5 * (min(p[0] for p in pts) + max(p[0] for p in pts))
    y0 = 0.5 * (min(p[1] for p in pts) + max(p[1] for p in pts))

    def to(p):
        return cx + scale * (p[0] - x0), cy - scale * (p[1] - y0)

    return [to(p) for p in mol["qa"]], [to(p) for p in mol["qw"]]


def draw_rules(W=846.0, H=349.0):
    mol = f1.molecule()
    ex = f1.examples()
    els, geo = mol["els"], mol["an"]["geo"]
    X0, X1 = 196.0, 652.0                      # the card of the three cases
    TX, TABX = X0 + 16.0, 500.0                # text edge; left edge of the tables
    HEAD = 11.3                                # 32 pt
    b = [f'<rect x="{X0}" y="0" width="{X1 - X0}" height="{H}" rx="6" fill="{CARD}"/>']
    scale, cy = 66.0, 136.0
    r_atom = {"S": 14.0, "H": 10.0}

    # ---- zone 1: the WF centres as computed --------------------------------
    A, Wc = molecule_xy(mol, 86.0, cy, scale)
    b.append(text(2, 20, "H_{2}S_{2}", BODY, "start", 600))
    for i, j in mol["bonds"]:
        b.append(line(*A[i], *A[j], INK, 1.0))
    for (x, y), e in zip(A, els):
        b.append(atom(x, y, r_atom[e], e))
    for (x, y), kind in zip(Wc, mol["kinds"]):
        b.append(circle(x, y, 4.6, INK if kind == "bond" else "#fff", INK, 1.2))
    n_bond = sum(k == "bond" for k in mol["kinds"])
    n_lone = sum(k == "lone" for k in mol["kinds"])
    assert n_bond + n_lone == len(geo)
    prov("zone 1: WFs of H2S2", f"{n_lone} one-centre + {n_bond} bond pairs",
         f"{mol['dir'].relative_to(f1.REPO)} (frozen plugin)")
    ly = 288.0
    b.append(circle(7, ly - 3.5, 4.6, "#fff", INK, 1.2))
    b.append(text(17, ly, "pair on one atom: it stays", 9.6))
    b.append(circle(7, ly + 21 - 3.5, 4.6, INK, INK, 1.2))
    b.append(text(17, ly + 21, "pair between two atoms:", 9.6))
    b.append(text(17, ly + 34.5, "the atoms compete", 9.6))
    b.append(chevron(177, cy))

    # ---- row A: two elements, force ----------------------------------------
    win, lose = ex["force"]["ends"]
    assert (win["el"], lose["el"]) == ("S", "H") and ex["force"]["rule"] == "force"
    b.append(duo(TX, 25, "Different elements: ",
                 "the larger *Z*_{eff} / *d*^{2} takes the pair", HEAD))
    ya, xs, span = 78.0, 268.0, 198.0
    xh = xs + span
    xd = xs + span * win["d"] / (win["d"] + lose["d"])
    b.append(line(xs, ya, xh, ya, SOFT, 1.1))
    k = 9.6                                        # mm per e/A^2
    b.append(block_arrow(xd - 10, ya, -k * win["force"], RED, 6.4, 15.5, 12.5))
    b.append(block_arrow(xd + 10, ya, k * lose["force"], SOFT, 6.4, 15.5, 12.5))
    b.append(atom(xs, ya, 24, "S", winner=True))
    b.append(atom(xh, ya, 15, "H"))
    b.append(circle(xd, ya, 7.2, INK, INK, 0))
    b.append(dimension(xs, xd, ya + 31, "*d*_{S}"))
    b.append(dimension(xd, xh, ya + 31, "*d*_{H}"))
    b.append(table(TABX, 52, ("*Z*_{eff}", "*d* (Å)", "*Z*_{eff}/*d*^{2}"),
                   [("S", f"{win['zeff']:.2f}", f"{win['d']:.3f}", f"{win['force']:.2f}"),
                    ("H", f"{lose['zeff']:.2f}", f"{lose['d']:.3f}", f"{lose['force']:.2f}")]))
    b.append(text(TABX + 140, 116, "S–H in H_{2}S_{2}", SMALL, "end", 400, GREY))
    prov("row A: S-H pair of H2S2, Z_eff / d / force",
         f"S {win['zeff']:.4f} / {win['d']:.4f} / {win['force']:.4f};  "
         f"H {lose['zeff']:.4f} / {lose['d']:.4f} / {lose['force']:.4f}",
         f"{ex['force']['dir']}, WF {ex['force']['wf']}")

    # ---- row B: lone pairs on both ends, chi -------------------------------
    b.append(line(TX, 127, X1 - 12, 127, LINE, 0.7))
    win, lose = ex["chi"]["ends"]
    assert (win["el"], lose["el"]) == ("S", "As") and ex["chi"]["rule"] == "chi"
    b.append(duo(TX, 151, "Both carry lone pairs: ",
                 "the larger *χ* = *Z*_{eff} / *ρ*^{2} takes it instead", HEAD))
    yb, xl, span = 205.5, 266.0, 168.0
    xr = xl + span
    xd = xl + span * lose["d"] / (win["d"] + lose["d"])
    r_s = 43.0
    r_as = r_s * math.sqrt(lose["rho2"] / win["rho2"])
    dash = ' stroke-dasharray="3.4 2.8"'
    for x, r in ((xl, r_as), (xr, r_s)):
        b.append(circle(x, yb, r, "none", GREY, 0.9, dash))
    lobes = [(xl - 32.0, yb), (xr + 20.7, yb - 23.1), (xr + 20.7, yb + 23.1)]
    for x, y in lobes:
        b.append(circle(x, y, 11.0, BLOB, "none", 0, ' opacity="0.75"'))
    b.append(line(xl, yb, xr, yb, SOFT, 1.1))
    kc = 10.5                                      # mm per unit of chi
    b.append(block_arrow(xd + 10, yb, kc * win["chi"], RED, 6.4, 15.5, 12.5))
    b.append(block_arrow(xd - 10, yb, -kc * lose["chi"], SOFT, 6.4, 15.5, 12.5))
    b.append(atom(xl, yb, 22, "As"))
    b.append(atom(xr, yb, 22, "S", winner=True))
    for x, y in lobes:
        b.append(circle(x, y, 4.6, "#fff", INK, 1.2))
    b.append(circle(xd, yb, 7.2, INK, INK, 0))
    # radius of the lone-pair cloud
    for x, r, ang, lab in ((xl, r_as, 60.0, "*ρ*_{As}"), (xr, r_s, 120.0, "*ρ*_{S}")):
        ca, sa = math.cos(math.radians(ang)), math.sin(math.radians(ang))
        b.append(line(x + 22 * ca, yb - 22 * sa, x + r * ca, yb - r * sa, GREY, 0.9))
        b.append(text(x + (r + 3.5) * ca + (1.0 if ca > 0 else -1.0),
                      yb - (r + 3.5) * sa + 1.0, lab, SMALL,
                      "start" if ca > 0 else "end", 400, GREY))
    b.append(table(TABX, 180, ("*Z*_{eff}", "*ρ*^{2} (Å^{2})", "*χ*"),
                   [("S", f"{win['zeff']:.2f}", f"{win['rho2']:.2f}", f"{win['chi']:.2f}"),
                    ("As", f"{lose['zeff']:.2f}", f"{lose['rho2']:.2f}", f"{lose['chi']:.2f}")]))
    b.append(text(TABX + 140, 244, "As–S in As_{4}S_{4}", SMALL, "end", 400, GREY))
    prov("row B: As-S pair of As4S4, Z_eff / rho^2 / chi",
         f"S {win['zeff']:.4f} / {win['rho2']:.4f} / {win['chi']:.4f};  "
         f"As {lose['zeff']:.4f} / {lose['rho2']:.4f} / {lose['chi']:.4f}",
         f"{ex['chi']['dir']}, WF {ex['chi']['wf']}")

    # ---- row C: one element, shared ----------------------------------------
    b.append(line(TX, 257, X1 - 12, 257, LINE, 0.7))
    ss = [g for g in geo if g["contested"] and len(g["hosts"]) == 2]
    assert len(ss) == 1 and all(abs(v - 0.5) < 2e-3 for v in ss[0]["hosts"].values())
    assert {els[a] for a in ss[0]["hosts"]} == {"S"}
    b.append(duo(TX, 289.5, "Same element: ", "the pair is shared", HEAD))
    yc, x1, x2 = 285.5, 464.0, 624.0
    xm = 0.5 * (x1 + x2)
    b.append(line(x1, yc, x2, yc, SOFT, 1.1))
    b.append(block_arrow(xm - 10, yc, -30, RED, 6.4, 15.5, 12.5))
    b.append(block_arrow(xm + 10, yc, 30, RED, 6.4, 15.5, 12.5))
    b.append(atom(x1, yc, 16, "S", winner=True))
    b.append(atom(x2, yc, 16, "S", winner=True))
    b.append(circle(xm, yc, 8.6, CARD, "none", 0))
    b.append(split_dot(xm, yc, 7.2, 90.0))
    n_pair = mol["an"]["n_pair"]
    share = round(n_pair * 0.5, 2)
    b.append(text(xm - 42, yc - 11, f"{share:g} e^{{−}}", SMALL, "middle", 400, GREY))
    b.append(text(xm + 42, yc - 11, f"{share:g} e^{{−}}", SMALL, "middle", 400, GREY))
    prov("row C: S-S pair of H2S2",
         f"shares {[round(v, 4) for v in sorted(ss[0]['hosts'].values())]}, "
         f"{n_pair} electrons per WF", str(mol["dir"].relative_to(f1.REPO)))

    # ---- Z_eff --------------------------------------------------------------
    b.append(f'<rect x="{X0 + 8}" y="312" width="{X1 - X0 - 16}" height="30" rx="4" fill="#fff"/>')
    screen = f1.plugin().PAIR_SCREEN
    assert abs(screen - 5 / 16) < 1e-12
    b.append(text(TX, 331.8, "*Z*_{eff} = *Z*_{val} − *n*_{inner} − 5/16", 11.6,
                  "start", 700, RED))
    b.append(text(404, 324.6, "*n*_{inner}: the atom’s own electrons inside *d*",
                  SMALL, "start", 400, INK))
    b.append(text(404, 336.6, "5/16: the other electron of the pair",
                  SMALL, "start", 400, INK))
    b.append(chevron(660, cy))

    # ---- zone 3: owners and oxidation states -------------------------------
    A, Wc = molecule_xy(mol, 761.0, cy, scale)
    s_atoms = [i for i, e in enumerate(els) if e == "S"]
    (xa, ya_), (xb, yb_) = A[s_atoms[0]], A[s_atoms[1]]
    r_disc = 0.5 * math.hypot(xa - xb, ya_ - yb_) + 2.0
    for i in s_atoms:
        b.append(circle(*A[i], r_disc, TINT, "none", 0))
    for i, j in mol["bonds"]:
        b.append(line(*A[i], *A[j], INK, 1.0))
    held = [0.0] * len(els)
    ang = math.degrees(math.atan2(-(yb_ - ya_), xb - xa)) + 90.0
    for (x, y), g in zip(Wc, geo):
        for a, f in g["hosts"].items():
            held[a] += n_pair * f
        assert all(els[a] == "S" for a in g["hosts"])
        if len(g["hosts"]) == 2:
            b.append(circle(x, y, 6.6, TINT, "none", 0))
            b.append(split_dot(x, y, 5.4, ang, 1.1))
        else:
            b.append(circle(x, y, 4.6, RED, RED, 0))
    for (x, y), e in zip(A, els):
        b.append(atom(x, y, r_atom[e], e, winner=(e == "S")))
    zval = mol["an"]["zval"]
    os_ = mol["os"]
    for a, e in enumerate(els):
        assert abs(held[a] - round(held[a])) < 5e-3
        held[a] = float(round(held[a]))
        assert int(round(zval[e] - held[a])) == os_[a]
    # oxidation-state labels, placed in the free corner next to each atom
    spot = {0: (14, 5, "start"), 1: (14, 5, "start"), 2: (21, -25, "start"),
            3: (-22, 36, "end")}
    for a, (dx, dy, anc) in spot.items():
        x, y = A[a]
        b.append(text(x + dx, y + dy, f"{os_[a]:+d}".replace("-", "−"), 13.5, anc, 700,
                      RED if els[a] == "S" else INK))
    for e in ("S", "H"):
        a = els.index(e)
        prov(f"zone 3: {e} of H2S2", f"Z_val {zval[e]:g} - {held[a]:g} e held = "
             f"{os_[a]:+d}", str(mol["dir"].relative_to(f1.REPO)))
    ty = 288.0
    b.append(text(674, ty, "oxidation state =", 9.6))
    b.append(text(674, ty + 13.8, "*Z*_{val} − electrons held", 9.6))
    row = []
    for e in ("S", "H"):
        a = els.index(e)
        row.append(f"{e}: {zval[e]:g} − {held[a]:g} = "
                   + f"{os_[a]:+d}".replace("-", "−"))
    b.append(text(674, ty + 34.5, row[0], 10.0, "start", 700, RED))
    b.append(text(674, ty + 48.5, row[1], 10.0, "start", 700, INK))
    return svg(W, H, b)


def table(x, y, head, rows, size=SMALL):
    """Element | two inputs | result (bold); the first row is the winner."""
    cols = (0.0, 52.0, 97.0, 140.0)          # right edges of the number columns
    big = BODY
    out = [text(x + c, y, h, size, "end", 400, GREY) for c, h in zip(cols[1:], head)]
    out.append(line(x, y + 4.4, x + cols[-1], y + 4.4, SOFT, 0.6))
    for n, (el, *vals) in enumerate(rows):
        yy = y + 18.5 + 16.5 * n
        top = n == 0
        if top:
            out.append(f'<rect x="{x - 3.5:.2f}" y="{yy - 11.4:.2f}" '
                       f'width="{cols[-1] + 7:.2f}" height="15.6" rx="2.2" fill="{TINT}"/>')
        out.append(text(x, yy, el, big, "start", 700, RED if top else INK))
        for c, v in zip(cols[1:-1], vals[:-1]):
            out.append(text(x + c, yy, v, big, "end", 400, INK))
        out.append(text(x + cols[-1], yy, vals[-1], big, "end",
                        700 if top else 400, RED if top else INK))
    return "".join(out)


def split_dot(x, y, r, angle, gap=1.3):
    """A pair shared by two atoms: two half discs `2 * gap` apart; `angle`
    (degrees, y up) is the direction of the cut."""
    ux, uy = math.cos(math.radians(angle)), -math.sin(math.radians(angle))
    nx, ny = -uy, ux
    out = []
    for s in (1.0, -1.0):
        pts = []
        for i in range(25):
            t = math.pi * i / 24
            pts.append((x + s * gap * nx + r * (math.cos(t) * ux + s * math.sin(t) * nx),
                        y + s * gap * ny + r * (math.cos(t) * uy + s * math.sin(t) * ny)))
        d = " ".join(f"{px:.2f},{py:.2f}" for px, py in pts)
        out.append(f'<polygon points="{d}" fill="{RED}"/>')
    return "".join(out)


# ------------------------------------------------------------------ accuracy
# Each row: algorithm, where its bond polarity comes from, results file
# (python3 -m pyleos bench clean --split dev --method <plugin> --out results/<file>).
# Plugins: pyleos/methods/ (this work, the two flows; greedy: --greedy) and
# baselines/ next to this script (midpoint.py, shifted_midpoint.py, voronoi.py,
# density_basins.py; the last two need numpy and scipy).
ACCURACY = [
    ("This work", "no table", "clean38_r3-T-B.tsv"),
    ("Shifted-midpoint rule", "no table", "clean38_shifted_midpoint.tsv"),
    ("Minimum-cost flow", "Allen table", "clean38_flow.tsv"),
    ("Model-density basins", "no table", "clean38_density_basins.tsv"),
    ("Greedy assignment", "Allen table", "clean38_greedy.tsv"),
    ("Midpoint rule", "no table", "clean38_midpoint.tsv"),
    ("Voronoi cells", "no table", "clean38_voronoi.tsv"),
    ("Minimum-cost flow", "distance only", "clean38_geo.tsv"),
]


def score(name):
    rows = list(csv.DictReader((f1.REPO / "results" / name).open(), delimiter="\t"))
    ok = sum(int(r["atoms_ok"]) for r in rows)
    n = sum(int(r["atoms_labelled"]) for r in rows)
    full = sum(r["atoms_ok"] == r["atoms_labelled"] for r in rows)
    return ok, n, full, len(rows)


def draw_accuracy(W=340.0, H=171.0):
    """Atoms assigned the reference oxidation state by the assignment
    algorithms, best first."""
    b = []
    x0, full_len, pitch, bh = 152.0, 118.0, 19.3, 11.0
    scores = [(lab, tag, name) + score(name) for lab, tag, name in ACCURACY]
    assert [r[3] for r in scores] == sorted((r[3] for r in scores), reverse=True)
    assert scores[0][0] == "This work" and scores[0][3] > scores[1][3]
    n_all = scores[0][4]
    b.append(text(0, 8.4, f"atoms assigned correctly, of {n_all} "
                  "(a small set of difficult cases)", SMALL, "start", 400, GREY))
    y = 15.0
    for lab, tag, name, ok, n, full, cases in scores:
        assert n == n_all and cases == scores[0][6]
        this = lab == "This work"
        b.append(f'<text x="0" y="{y + bh - 1.7:.2f}" font-size="9.6" '
                 f'font-weight="{700 if this else 400}" fill="{RED if this else INK}">'
                 f'{esc(lab)}<tspan font-size="{SMALL}" font-weight="400" fill="{GREY}">'
                 f'  {esc(tag)}</tspan></text>')
        b.append(f'<rect x="{x0}" y="{y:.2f}" width="{full_len}" height="{bh}" fill="{CARD}"/>')
        b.append(f'<rect x="{x0}" y="{y:.2f}" width="{full_len * ok / n:.2f}" '
                 f'height="{bh}" fill="{THIS if this else SOFT}"/>')
        b.append(text(x0 + full_len + 4.5, y + bh - 1.7, f"{ok}", 9.8, "start",
                      700 if this else 400, RED if this else INK))
        b.append(text(x0 + full_len + 29.0, y + bh - 1.7, f"{100 * ok / n:.1f}%", 9.8, "start",
                      700 if this else 400, RED if this else INK))
        prov(f"accuracy: {lab} ({tag})", f"{ok}/{n} atoms ({100 * ok / n:.1f}%), "
             f"{full}/{cases} cases fully correct", f"results/{name}")
        y += pitch
    return svg(W, H, b)


# -------------------------------------------------------------------- gpuw90
def box(x, y, w, h, title, sub):
    return (f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" rx="3.5" '
            f'fill="#fff" stroke="{RED}" stroke-width="1.2"/>'
            + text(x + 9, y + 14.2, title, 9.8, "start", 700, INK)
            + text(x + 9, y + 27.2, sub, SMALL, "start", 400, INK))


def draw_gpuw90(W=846.0, H=183.0):
    rows = {cfg: t for _tag, cfg, _kind, t, _src in g3.load_iteration_times()}
    ref, gpu = rows["wannier90, 1 core"], rows["gpuw90, GPU"]
    for p in g3.PROV:
        if "1 core" in p[1] or "gpuw90, GPU" in p[1] or "hand-copied" in p[1]:
            prov("gpuw90 bars: " + p[1], p[2], p[3])
    # gpuw90 on a CPU: the run whose thread count is on record
    local = f1.REPO / "results" / "gpuw90_cpu_bench_local.tsv"
    lines = [ln for ln in local.read_text().splitlines() if not ln.startswith("#")]
    recs = list(csv.DictReader(lines, delimiter="\t"))
    threads = {r["threads"] for r in recs}
    assert len(threads) == 1
    n_thr = int(threads.pop())
    mean = lambda st: sum(float(r["s_per_iter"]) for r in recs if r["stage"] == st) / \
        sum(r["stage"] == st for r in recs)
    cpu = {"dis": mean("disentangle"), "wan": mean("wannierise")}
    prov(f"gpuw90 bars: gpuw90, {n_thr}-core CPU: s/iteration dis | loc (speed-up)",
         f"{cpu['dis']:.3f} | {cpu['wan']:.3f}  ({ref['dis'] / cpu['dis']:.1f}x | "
         f"{ref['wan'] / cpu['wan']:.1f}x)", "results/gpuw90_cpu_bench_local.tsv")
    b = []

    # ---- column 1: the Wannier step ----------------------------------------
    w1 = 238.0
    b.append(chip(34, 0, 170, 22, "DFT wavefunctions", "soft"))
    b.append(thin_arrow(w1 / 2, 24, w1 / 2, 39))
    b.append(f'<rect x="0" y="41" width="{w1}" height="103" rx="5" fill="{TINT}"/>')
    b.append(text(9, 55.0, "gpuw90", 10.6, "start", 700, RED))
    b.append(text(w1 - 9, 55.0, "PyTorch, CPU or GPU", SMALL, "end", 400, INK))
    b.append(box(8, 61, w1 - 16, 33, "Disentanglement",
                 "keep *N*_{WF} states at every *k*-point"))
    b.append(thin_arrow(w1 / 2, 95.2, w1 / 2, 102.8, RED))
    b.append(box(8, 104, w1 - 16, 33, "Localisation",
                 "rotate them to the smallest spread"))
    b.append(thin_arrow(w1 / 2, 146, w1 / 2, 159))
    b.append(chip(14, 161, 210, 22, "WF centres and spreads", "plain"))

    # ---- column 2: one batched operation; the loop -------------------------
    c2 = 268.0
    s, y0 = 50.0, 21.0
    yc = y0 + s / 2
    b.append(text(c2, 10, "wannier90", 9.6, "start", 600, GREY))
    b.append(f'<rect x="{c2 + 2}" y="{y0}" width="{s}" height="{s}" fill="#fff" '
             f'stroke="{GREY}" stroke-width="1.2"/>')
    b.append(text(c2 + 2 + s / 2, yc + 3.6, "*M*^{(*k*,*b*)}", BODY, "middle", 400, GREY))
    xa, r = c2 + 2 + s + 15.0, 12.5
    b.append(f'<path d="M{xa - 7.5:.2f},{yc - 10:.2f} A{r},{r} 0 1 1 '
             f'{xa - 7.5:.2f},{yc + 10:.2f}" fill="none" stroke="{GREY}" '
             f'stroke-width="1.2"/>')
    b.append(f'<polygon points="{xa - 12.6:.2f},{yc + 10.2:.2f} {xa - 5.2:.2f},'
             f'{yc + 5.4:.2f} {xa - 5.0:.2f},{yc + 14.2:.2f}" fill="{GREY}"/>')
    b.append(text(c2 + 2, y0 + s + 13, "one (*k*, *b*) at a time", SMALL, "start", 400, GREY))
    b.append(block_arrow(c2 + 104, yc, 26, RED, t=5.2, hw=13, hl=10.5))
    x1, n, off = c2 + 144.0, 5, 4.6
    b.append(text(x1, 10, "gpuw90", 9.6, "start", 700, RED))
    for i in reversed(range(n)):
        b.append(f'<rect x="{x1 + i * off:.2f}" y="{y0 + (n - 1 - i) * off:.2f}" '
                 f'width="{s}" height="{s}" fill="{TINT if i else "#fff"}" stroke="{RED}" '
                 f'stroke-width="1.2"/>')
    yf = y0 + (n - 1) * off
    b.append(text(x1 + s / 2, yf + s / 2 + 3.6, "*M*", BODY, "middle", 600, RED))
    xt = x1 + s + (n - 1) * off + 8
    b.append(text(xt, yc + 2.0, "all (*k*, *b*)", SMALL, "start", 600, INK))
    b.append(text(xt, yc + 14.2, "in one call", SMALL, "start", 600, INK))
    b.append(text(x1, yf + s + 12.5, "[*N*_{*k*}, *N*_{*b*}, *N*_{WF}, *N*_{WF}]", SMALL,
                  "start", 400, INK))
    # the four tensor operations of one localisation iteration
    ty = 118.0
    b.append(text(c2, ty, "One localisation iteration, every *k*-point at once", 9.6,
                  "start", 600, INK))
    steps = [("rotate *M*", "matmul"), ("gradient", "elementwise"),
             ("exp(Δ*W*)", "eigh"), ("update *U*", "matmul")]
    cw, cg, ch = 64.0, 12.0, 24.0
    for i, (lab, op) in enumerate(steps):
        x = c2 + i * (cw + cg)
        b.append(chip(x, ty + 8, cw, ch, lab, "plain", 9.0))
        b.append(text(x + cw / 2, ty + 8 + ch + 11.5, op, SMALL, "middle", 400, GREY))
        if i < len(steps) - 1:
            b.append(thin_arrow(x + cw + 1.6, ty + 8 + ch / 2, x + cw + cg - 1.6,
                                ty + 8 + ch / 2, GREY, 1.2, 3.4))
    xr, xl_, yl = c2 + 3 * (cw + cg) + cw / 2, c2 + cw / 2, ty + 8 + ch + 21.0
    b.append(f'<path d="M{xr:.2f},{yl - 5:.2f} V{yl:.2f} H{xl_:.2f}" fill="none" '
             f'stroke="{GREY}" stroke-width="1.2"/>')
    b.append(thin_arrow(xl_, yl, xl_, yl - 6.5, GREY, 1.2, 3.4))
    c2_right = c2 + 4 * cw + 3 * cg

    # ---- column 3: timings and agreement with wannier90 --------------------
    c3 = 596.0
    x0, mm_per_s, bh = c3 + 86.0, 112.0, 8.4
    y = 10.0
    b.append(text(W, y + 0.2, "seconds per iteration", SMALL, "end", 400, GREY))
    for title, key in (("Disentanglement", "dis"), ("Localisation", "wan")):
        b.append(text(c3, y, title, 9.6, "start", 600, INK))
        for n_row, (lab, t, col) in enumerate((("wannier90, 1 core", ref, SOFT),
                                               (f"gpuw90, {n_thr} cores", cpu, "#E3A6B0"),
                                               ("gpuw90, 1 GPU", gpu, RED))):
            yy = y + 6.5 + 11.4 * n_row
            w = mm_per_s * t[key]
            b.append(text(c3, yy + bh - 1.5, lab, SMALL, "start", 400, INK))
            b.append(f'<rect x="{x0}" y="{yy:.2f}" width="{w:.2f}" height="{bh}" '
                     f'fill="{col}"/>')
            b.append(text(x0 + w + 3, yy + bh - 1.5, f"{t[key]:.3f}", SMALL, "start",
                          400, INK))
            if col == RED:
                b.append(text(x0 + w + 30, yy + bh - 1.5,
                              f"{ref[key] / t[key]:.1f}× faster", SMALL + 0.8, "start",
                              700, RED))
            elif col != SOFT:
                b.append(text(x0 + w + 30, yy + bh - 1.5,
                              f"{ref[key] / t[key]:.1f}× faster", SMALL, "start", 400, INK))
        y += 51.0
    b.append(text(c3, y - 1.5, "96-atom Li_{2}MnO_{3} cell, 288 WFs, one V100 GPU", SMALL,
                  "start", 400, GREY))
    # agreement with wannier90 (gpuw90/README.md; reproduced with the check scripts)
    readme = (f1.REPO / "gpuw90" / "README.md").read_text()
    m = re.search(r"GaAs (8\.3)e-7 Å, Diamond (3\.1)e-9 Å", readme)
    m2 = re.search(r"\|ΔΩ_I\| = (4\.4)e-9", readme)
    assert m and m2, "gpuw90/README.md: the acceptance lines changed"
    y = 130.0
    b.append(f'<rect x="{c3 - 6}" y="{y - 12}" width="{W - c3 + 6}" height="65" rx="4" '
             f'fill="{CARD}"/>')
    b.append(text(c3 + 2, y, "Same result as wannier90", 9.6, "start", 700, INK))
    checks = [("GaAs", f"centres within {m.group(1)} × 10^{{−7}} Å"),
              ("Diamond", f"centres within {m.group(2)} × 10^{{−9}} Å"),
              ("Li_{2}MnO_{3}", f"spread Ω_{{I}} within {m2.group(1)} × 10^{{−9}} Å^{{2}}")]
    for i, (name, val) in enumerate(checks):
        yy = y + 15.0 + 13.0 * i
        b.append(text(c3 + 2, yy, name, SMALL, "start", 600, INK))
        b.append(text(W - 6, yy, val, SMALL, "end", 400, INK))
    prov("gpuw90: agreement with wannier90",
         f"GaAs {m.group(1)}e-7 A, diamond {m.group(2)}e-9 A (max centre deviation); "
         f"Li2MnO3 |dOmega_I| {m2.group(1)}e-9 A^2",
         "gpuw90/README.md; check_against_w90.py testdata/{gaas,diamond}, "
         "check_dis_against_w90.py testdata/li2mno3 --seed wannier90.1")
    assert c2_right < c3 - 20
    return svg(W, H, b)


def draw_qr(size=52.0, quiet=4):
    """QR code of the repository link, dark modules on white, `quiet` modules
    of margin (the matrix file carries one already)."""
    rows = [ln.strip() for ln in (HERE / "assets" / "qr_leos.txt").read_text().split()]
    n = len(rows)
    assert n and all(len(r) == n and set(r) <= {"0", "1"} for r in rows)
    extra = quiet - 1
    m = size / (n + 2 * extra)
    b = [f'<rect x="0" y="0" width="{size}" height="{size}" rx="2.2" fill="#fff"/>']
    d = []
    for j, r in enumerate(rows):
        i = 0
        while i < n:
            if r[i] == "1":
                k = i
                while k < n and r[k] == "1":
                    k += 1
                d.append(f"M{(i + extra) * m:.3f},{(j + extra) * m:.3f}h{(k - i) * m:.3f}"
                         f"v{m:.3f}h{-(k - i) * m:.3f}z")
                i = k
            else:
                i += 1
    b.append(f'<path d="{"".join(d)}" fill="#111114" shape-rendering="crispEdges"/>')
    prov("QR code", f"{n - 2} x {n - 2} modules, {m:.2f} mm each, quiet zone {quiet} modules",
         "assets/qr_leos.txt (make_qr.swift, level Q)")
    return svg(size, size, b)


def main():
    OUT.mkdir(exist_ok=True)
    for name, draw in (("basis", draw_basis), ("pipeline", draw_pipeline), ("rules", draw_rules),
                       ("accuracy", draw_accuracy), ("gpuw90", draw_gpuw90), ("qr", draw_qr)):
        (OUT / f"{name}.svg").write_text(draw())
    w = max(len(p[0]) for p in PROV)
    for what, value, source in PROV:
        print(f"  {what:{w}s}  {value}   [{source}]")
    print(f"schematics -> {OUT}/basis.svg, pipeline.svg, rules.svg, accuracy.svg, gpuw90.svg")


if __name__ == "__main__":
    main()
