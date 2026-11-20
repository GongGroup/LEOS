#!/bin/sh
# Poster build: panels and schematics (SVG) -> PDF (900 x 1200 mm, one page, vector)
# -> PNG preview.
#   sh build.sh            everything
#   sh build.sh --no-panels   reuse panels/*.svg
#   sh build.sh --measure     also print the layout report (boxes in mm, overflow flags)
# The figures of the repository README are cut from the PDF afterwards:
#   python3 make_readme_figs.py   (writes readme/*.png, which are tracked in git)
# Needs: python3 with numpy, matplotlib, PyMuPDF; the headless Chromium of Playwright
# (override with CHROME=/path/to/chrome-headless-shell).
set -e
cd "$(dirname "$0")"
HERE=$(pwd)
OUT_PDF="$HERE/../LEOS_poster_RSC_FIRST_2026.pdf"
OUT_PNG="$HERE/../LEOS_poster_RSC_FIRST_2026.png"
CHROME=${CHROME:-"$HOME/Library/Caches/ms-playwright/chromium_headless_shell-1228/chrome-headless-shell-mac-arm64/chrome-headless-shell"}
FLAGS="--no-sandbox --disable-gpu --allow-file-access-from-files --virtual-time-budget=8000"

case " $* " in *" --no-panels "*) ;; *) python3 make_panels.py; python3 make_schematics.py ;; esac

"$CHROME" $FLAGS --no-pdf-header-footer --print-to-pdf="$OUT_PDF" "file://$HERE/poster.html" 2>&1 \
  | grep -v -E "^$|bytes written" || true

case " $* " in *" --measure "*)
  "$CHROME" $FLAGS --dump-dom "file://$HERE/poster.html#measure" 2>/dev/null \
    | sed -n '/<pre id="measure">/,/<\/pre>/p' | sed -e 's/<[^>]*>//g' -e 's/&amp;/\&/g' ;;
esac

python3 - "$OUT_PDF" "$OUT_PNG" <<'PY'
import sys, fitz
doc = fitz.open(sys.argv[1]); page = doc[0]
w, h = page.rect.width / 72 * 25.4, page.rect.height / 72 * 25.4
assert len(doc) == 1, f"{len(doc)} pages: content overflows the page"
assert abs(w - 900) < 0.5 and abs(h - 1200) < 0.5, (w, h)
zoom = 3000 / page.rect.height                  # preview 3000 px tall
pix = page.get_pixmap(matrix=fitz.Matrix(zoom, zoom), alpha=False)
pix.save(sys.argv[2])
print(f"{sys.argv[1]}: 1 page, {w:.0f} x {h:.0f} mm")
print(f"{sys.argv[2]}: {pix.width} x {pix.height} px")
PY
