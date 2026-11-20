#!/usr/bin/env python3
"""Figures of the repository README, cut from the poster PDF.

Run:  python3 make_readme_figs.py        (after sh build.sh; needs PyMuPDF)
Reads ../LEOS_poster_RSC_FIRST_2026.pdf (900 x 1200 mm, vector) and writes
readme/*.png: one image per poster block, plus a preview of the whole poster.

The README shows rasters and not panels/*.svg because the schematics are set
in Avenir Next and have no background: on GitHub they would fall back to
another font and vanish in the dark theme.  These PNGs are the only rendered
figures kept in git; everything else under paper/ is regenerated.

BOXES are in poster millimetres, taken from the layout report of
`sh build.sh --measure` (image boxes of #intro, #rule, #s2-#s5).  ANCHORS
guards them: each heading must still sit where the boxes assume, otherwise
the script stops and the boxes have to be read from the report again.
"""
import sys
from pathlib import Path

import fitz                                           # PyMuPDF

HERE = Path(__file__).resolve().parent
PDF = HERE.parent / "LEOS_poster_RSC_FIRST_2026.pdf"
OUT = HERE / "readme"
MM = 72 / 25.4                                        # PDF points per mm

# name: (x, y, width, height in mm, width of the PNG in pixels)
BOXES = {
    "basis":    (497.0, 172.0, 368.0, 138.0, 1300),   # basis sets and pipeline
    "rule":     (27.0, 350.7, 846.0, 349.0, 2100),    # 1  the assignment rule
    "hydrogen": (27.0, 733.5, 251.7, 154.7, 1000),    # 2  X-H bonds
    "accuracy": (294.8, 733.5, 340.0, 171.0, 1360),   # 3  eight algorithms
    "chi":      (650.9, 733.5, 222.1, 170.2, 900),    # 4  in-situ chi vs Allen
    "gpuw90":   (27.0, 966.8, 846.0, 183.0, 2100),    # 5  gpuw90
    "poster":   (0.0, 0.0, 900.0, 1200.0, 1500),      # whole poster
}
PAD = {"basis": 0.0, "poster": 0.0}                   # others: 2 mm of white

# heading text: (left, top) of its first glyph box in mm, tolerance 1.5 mm
ANCHORS = {
    "Why Wannier functions?": (39.0, 184.4),
    "Solving the hydrogen problem": (47.6, 711.6),
    "Comparison of assignment algorithms": (315.4, 711.6),
    "Emergent electronegativity": (671.5, 711.6),
    "gpuw90: computing Wannier functions on a GPU": (47.6, 927.5),
}


def main():
    doc = fitz.open(PDF)
    page = doc[0]
    w, h = page.rect.width / MM, page.rect.height / MM
    assert len(doc) == 1 and abs(w - 900) < 0.5 and abs(h - 1200) < 0.5, (w, h)
    bad = []
    for text, (x, y) in ANCHORS.items():
        hits = page.search_for(text)
        got = (hits[0].x0 / MM, hits[0].y0 / MM) if hits else None
        if got is None or abs(got[0] - x) > 1.5 or abs(got[1] - y) > 1.5:
            bad.append(f"  {text!r}: expected ({x}, {y}) mm, found {got}")
    if bad:
        sys.exit("poster layout changed; update BOXES and ANCHORS from "
                 "`sh build.sh --measure`:\n" + "\n".join(bad))
    OUT.mkdir(exist_ok=True)
    for name, (x, y, bw, bh, px) in BOXES.items():
        p = PAD.get(name, 2.0)
        clip = fitz.Rect((x - p) * MM, (y - p) * MM,
                         (x + bw + p) * MM, (y + bh + p) * MM) & page.rect
        zoom = px / clip.width
        pix = page.get_pixmap(matrix=fitz.Matrix(zoom, zoom), clip=clip,
                              alpha=False)
        path = OUT / f"{name}.png"
        pix.save(path)
        print(f"  {path.relative_to(HERE.parent.parent)}  "
              f"{pix.width} x {pix.height} px  {path.stat().st_size / 1024:.0f} kB")


if __name__ == "__main__":
    main()
