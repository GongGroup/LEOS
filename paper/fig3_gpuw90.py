#!/usr/bin/env python3
"""Figure 3 - Accelerating the Wannierisation step.

(a) Where each stage of the per-frame pipeline runs: the in-situ reference
    and the six modifications M1-M6; on the right, the internal change made
    by gpuw90 (M4).
(b) Wall time per iteration of disentanglement and of localisation: serial
    wannier90, wannier90 with MPI (M1) and with MPI + threads (M1 + M2),
    gpuw90 (M4) on CPU and on one V100 GPU.
(c) Wall time per frame of the Wannier step and of the GPU job for the
    in-situ reference and for M1, M1 + M2, M1 + M3.
(d) wannier90 run time with and without the SCDM-only mode (M5), and time
    per analysed molecular-dynamics step with and without warm start (M6).

Run from this directory:   python3 fig3_gpuw90.py
Needs numpy + matplotlib + the standard library and the repository itself
(pyleos, pure standard library, is called once to report - not to plot - the
oxidation states of the two SCDM-only cases of panel (d)).

Every plotted number is read from a file under results/ or cases/, with one
exception: the gpuw90 timings of panel (b) exist only in gpuw90/README.md and
are kept in the hand-copied block below, which the script checks against the
README.
(The same block holds one documented SCDM-only timing that no data file
supports; it is drawn only when SHOW_DOCUMENTED_SCDM is switched on.)
No estimate or extrapolation is plotted.
Writes fig3_gpuw90.pdf / .png (page 17.1 cm x 11.4 cm) and prints a
provenance table.

COLOUR AND MARK MAPPING (fixed for the whole figure)
    grey  style.C["ref"],  square marker - reference configuration: serial
          wannier90 library inside the GPU job, full iterative
          Wannierisation, warm start off
    blue  style.C["flow"], round marker  - with the modification(s) named by
          the tag(s) M1-M6 next to the mark
    solid / open bar (panel c)           - Wannier step / GPU job
    dashed outline (panel a)             - stage whose iterations are skipped
    open marker, dashed line (panel d)   - documented value without a data
          file (drawn only with SHOW_DOCUMENTED_SCDM = True)
    All text is ink; colour is never the only cue (every mark is labelled
    with its tag or sits in a labelled row).
"""
import csv
import re
import statistics
import sys
from pathlib import Path

import numpy as np

import style
from style import C, CM, formula

import matplotlib.pyplot as plt                              # noqa: E402
from matplotlib.lines import Line2D                          # noqa: E402
from matplotlib.patches import FancyArrowPatch, Rectangle    # noqa: E402
from matplotlib.ticker import (FixedLocator, FuncFormatter,  # noqa: E402
                               MultipleLocator, NullFormatter, NullLocator)

REPO = style.REPO
sys.path.insert(0, str(REPO))
from pyleos.parsers import load_case                         # noqa: E402
from pyleos.transport import assign_case                     # noqa: E402

DATA = REPO / "cases"                   # <case>/<run>/ directories
RESULTS = REPO / "results"
E2E = RESULTS / "w90_e2e_pipeline.tsv"
W90_SCAN = RESULTS / "w90_threads_scan.tsv"
W90_REF_WOUT = REPO / "gpuw90" / "testdata" / "li2mno3" / "wannier90.1_ref.wout"
GPUW90_README = REPO / "gpuw90" / "README.md"
RESULTS_MD = REPO / "docs" / "RESULTS.md"
WS_DIR = RESULTS / "cathode_md"         # warmstart_{INCAR,loop}.{on,off}
# SCDM-only runs and the full Wannierisations of the same cells.  Only the
# pairs with identical POSCAR, k grid, number of bands, number of WFs and
# spin treatment are used (CSe2 and GaSe); see fig3_caption.md for the four
# v3-scdm-only runs that have no like-for-like partner.
SCDM_PAIRS = [
    ("CSe2", DATA / "cse2" / "v31-wp", DATA / "cse2" / "v3-scdm-only"),
    ("GaSe", DATA / "sega" / "v31-wp", DATA / "sega" / "v3-scdm-only"),
]

# Panel (b): rows of W90_SCAN that are shown (arm, np, nt); spin channel 1.
SCAN_REF = ("seq", 1, 1)
SCAN_MPI = ("thr", 8, 1)
SCAN_MPI_THREADS = ("thr", 8, 4)
# Panel (d): also draw the SCDM-only timing that docs/RESULTS.md quotes for
# the 96-atom Li2MnO3 frame (no data file; open markers, dashed connector).
SHOW_DOCUMENTED_SCDM = False

# ==========================================================================
# HAND-COPIED VALUES - the only numbers in this script that are not read
# from a data file.
#
# (1) gpuw90 timings, plotted in panel (b).
# Source: gpuw90/README.md, section "性能实测", subsection "V100（SAI 集群）":
#   | 阶段 | V100 | 计算节点 torch-CPU | 集群 wannier90.x 串行 | V100 相对串行 |
#   | disentangle | 0.079 | 0.423 | ~1.10 | ~14× |
#   | wannierise | 0.041 | 0.303 | ~0.92 | ~22× |
# Measured 2026-07-11, SAI job 644120 (gpuw90/run_bench_v100.sbatch);
# li2mno3 frame, spin channel 1, 288 WFs, 356 bands, 8 k-points; seconds per
# iteration = stage wall time / fixed number of iterations (bench_gpu.py).
README_TIMINGS = {
    "gpuw90_cpu": {"dis": 0.423, "wan": 0.303},
    "gpuw90_v100": {"dis": 0.079, "wan": 0.041},
}
README_ROWS = [r"\|\s*disentangle\s*\|\s*0\.079\s*\|\s*0\.423\s*\|",
               r"\|\s*wannierise\s*\|\s*0\.041\s*\|\s*0\.303\s*\|"]
#
# (2) SCDM-only timing of the Li2MnO3 frame, docs/RESULTS.md section 9
# ("SCDM-only 快速模式", table row "96 原子 Li2MnO3 脱锂帧，12 对"):
#   "降低 67 倍（554 → 8.3 秒/帧）"
# No file under results/ or cases/ holds these two numbers.  They are drawn
# only when SHOW_DOCUMENTED_SCDM is True, and then with open markers.
DOCUMENTED_SCDM = {"label": "Li2MnO3", "full": 554.0, "scdm": 8.3}
DOCUMENTED_SCDM_ROW = r"554\s*→\s*8\.3\s*秒/帧"
# ==========================================================================

# Statements of docs/RESULTS.md section 9 that the caption quotes or that
# the provenance table compares with the data (pattern -> what it is).
RESULTS_CHECKS = [
    (r"元素平均 32/38（84%），逐原子 19/38",
     "SCDM-only vs full, 38 cases: element means 32/38, per atom 19/38"),
    (r"元素平均 0/12", "SCDM-only vs full, Li2MnO3 frames: 0/12"),
    (DOCUMENTED_SCDM_ROW, "SCDM-only, Li2MnO3 frame: 554 -> 8.3 s per frame"),
    (r"CO2：2 秒 → 0\.023 秒", "SCDM-only, CO2: 2 s -> 0.023 s"),
    (r"38 案例对比的原始数据不在仓库中",
     "raw data of the 38-case comparison are not in the repository"),
    (r"\| 开启热启动 \| 572 秒 \| 6742\.6 秒 \|", "warm start on: 572 s, 6742.6 s"),
    (r"\| 关闭热启动 \| 546 秒 \| 6532\.4 秒 \|", "warm start off: 546 s, 6532.4 s"),
    (r"快约 14 倍和 22 倍", "gpuw90 on V100: about 14x and 22x"),
    (r"快约 4\.6 倍", "gpuw90 on CPU: about 4.6x"),
    (r"两个自旋的 Wannier 部分可降到约半分钟",
     "gpuw90 extrapolation to about half a minute (estimate, not plotted)"),
]

COL_REF, COL_MOD = C["ref"], C["flow"]
MK_REF, MK_MOD = "s", "o"
NWF = r"$N_{\mathrm{WF}}$"

PROV = []       # (panel, quantity, value, source)


def prov(panel, quantity, value, source):
    PROV.append((panel, quantity, str(value), str(source)))


def rel(p):
    return str(Path(p).relative_to(REPO))


def read_tsv(path):
    """Tab-separated table; lines starting with '#' are comments."""
    with open(path, newline="") as fh:
        return list(csv.DictReader((ln for ln in fh if not ln.startswith("#")),
                                   delimiter="\t"))


def doc_line(path, pattern):
    """1-based number of the first line of `path` matching `pattern`."""
    for i, line in enumerate(Path(path).read_text().splitlines(), 1):
        if re.search(pattern, line):
            return i
    return None


def tag(*names):
    """Bold modification tag(s): tag('M1', 'M2') -> 'M1 + M2' in bold."""
    return r"$\mathbf{" + r"\,+\,".join(names) + "}$"


# ------------------------------------------------------------------ data --
def load_iteration_times():
    """Panel (b): seconds per iteration, spin channel 1 of the Li2MnO3 frame."""
    scan = read_tsv(W90_SCAN)

    def pick(key):
        arm, np_, nt = key
        r, = [r for r in scan if (r["arm"], r["np"], r["nt"], r["spin"]) ==
              (arm, str(np_), str(nt), "1")]
        src = f"{rel(W90_SCAN)} (arm {arm}, np {np_}, nt {nt}, spin 1)"
        return {"dis": float(r["dis_s_per_it"]),
                "wan": float(r["wan_s_per_it"])}, src

    readme = rel(GPUW90_README) + " (V100 table, job 644120) [hand-copied]"
    rows = [
        ("Reference", "wannier90, 1 core", "ref") + pick(SCAN_REF),
        (tag("M1"), "wannier90, 8 cores", "mod") + pick(SCAN_MPI),
        (tag("M1", "M2"), "wannier90, 32 cores", "mod")
        + pick(SCAN_MPI_THREADS),
        (tag("M4"), "gpuw90, CPU", "mod", README_TIMINGS["gpuw90_cpu"], readme),
        (tag("M4"), "gpuw90, GPU", "mod", README_TIMINGS["gpuw90_v100"], readme),
    ]
    base = rows[0][3]
    plain = ["reference", "M1", "M1 + M2", "M4", "M4"]
    for (_, cfg, _, t, src), name in zip(rows, plain):
        prov("b", f"{name}: {cfg}: s/iteration dis | loc (speed-up)",
             f"{t['dis']:.3f} | {t['wan']:.3f}  ({base['dis'] / t['dis']:.1f}x"
             f" | {base['wan'] / t['wan']:.1f}x)", src)
    txt = GPUW90_README.read_text()
    ok = all(re.search(p, txt) for p in README_ROWS)
    prov("b", "hand-copied gpuw90 rows still present in the README",
         "yes" if ok else "NO - UPDATE README_TIMINGS",
         f"{rel(GPUW90_README)}:{doc_line(GPUW90_README, README_ROWS[0])}-"
         f"{doc_line(GPUW90_README, README_ROWS[1])}")
    if not ok:
        raise SystemExit("gpuw90/README.md no longer contains the rows copied "
                         "into README_TIMINGS - update the hand-copied block")
    # spin channel 2 of the same scan rows (not plotted)
    for key, name in ((SCAN_REF, "reference"), (SCAN_MPI_THREADS, "M1 + M2")):
        r, = [r for r in scan if (r["arm"], r["np"], r["nt"], r["spin"]) ==
              (key[0], str(key[1]), str(key[2]), "2")]
        prov("b", f"{name}, spin channel 2: s/iteration dis | loc; dis "
             "iterations (not plotted)",
             f"{r['dis_s_per_it']} | {r['wan_s_per_it']}; {r['dis_iter']}",
             rel(W90_SCAN))
    # serial baseline of the README (same workload, another run; not plotted)
    if W90_REF_WOUT.is_file():
        w = W90_REF_WOUT.read_text(errors="replace")
        n_dis = len(re.findall(r"^\s+\d+\s+\S+\s+\S+\s+\S+\s+\S+\s+<-- DIS$",
                               w, flags=re.M))
        t_dis = float(re.search(r"Time to disentangle bands\s+([\d.]+)", w)[1])
        t_wan = float(re.search(r"Time for wannierise\s+([\d.]+)", w)[1])
        conv = re.findall(r"^\s+(\d+)\s+\S+\s+\S+\s+\S+\s+([\d.]+)\s+<-- CONV$",
                          w, flags=re.M)
        n_wan = int(conv[-1][0])
        prov("b", "serial wannier90 run used by the README: stage time / "
             "iterations, dis | loc (not plotted)",
             f"{t_dis:.3f} s / {n_dis} = {t_dis / n_dis:.3f} | "
             f"{t_wan:.3f} s / {n_wan} = {t_wan / n_wan:.3f}",
             rel(W90_REF_WOUT))
    return rows


def load_pipeline():
    """Panel (c): rows A-D of the per-frame pipeline table."""
    rows = {r["row"]: r for r in read_tsv(E2E)}
    assert list(rows) == ["A", "B", "C", "D"], "unexpected rows in " + rel(E2E)
    for r, must in (("A", "in-situ"), ("B", "offline w90"),
                    ("C", "row B + threaded MKL"),
                    ("D", "row B + per-spin num_wann")):
        assert must in rows[r]["config"], f"row {r} of {rel(E2E)} has changed"
    labels = {"A": "Reference", "B": tag("M1"), "C": tag("M1", "M2"),
              "D": tag("M1", "M3")}
    kind = {"A": "ref", "B": "mod", "C": "mod", "D": "mod"}
    out = []
    for key, r in rows.items():
        w, g = float(r["w90_wall_s"]), float(r["gpu_leg_s"])
        out.append((labels[key], kind[key], w, g))
        prov("c", f"row {key} ({r['config']}): wannier90 wall | GPU job (s)",
             f"{w:.1f} | {g:.1f}", rel(E2E))
    a, b = out[0], out[1]
    prov("c", "derived: GPU-job time saved by M1; whole frame A vs B "
         "(GPU job + wannier90)",
         f"{100 * (1 - b[3] / a[3]):.1f} %; {a[3]:.1f} vs {b[3] + b[2]:.1f} s",
         rel(E2E))
    return out


def parse_wout(path):
    """Problem size and timings of a wannier90 .wout (last occurrence)."""
    w = Path(path).read_text(errors="replace")

    def last(pattern, cast):
        return cast(re.findall(pattern, w)[-1])

    return {
        "nwf": last(r"Number of Wannier Functions\s*:\s*(\d+)", int),
        "nb": last(r"Number of input Bloch states\s*:\s*(\d+)", int),
        "nk": last(r"Total points\s*=\s*(\d+)", int),
        "n_iter": [int(v) for v in re.findall(
            r"Total number of iterations\s*:\s*(\d+)", w)[-2:]],
        "t_dis": last(r"Time to disentangle\s+([\d.]+)", float),
        "t_wan": last(r"Time for wannierise\s+([\d.]+)", float),
        "t_tot": last(r"Total Execution Time\s+([\d.]+)", float),
    }


def os_summary(case_dir):
    """Element-mean oxidation states (min-cost flow, defaults) and strays."""
    case = load_case(case_dir)
    res = assign_case(case)
    el = case.structure.elements
    means = {}
    for e in dict.fromkeys(el):
        v = [res.os_int[i] for i, x in enumerate(el) if x == e]
        means[e] = sum(v) / len(v)
    return ("  ".join(f"{e} {m:+.2f}" for e, m in means.items())
            + f"  (stray WFs: {len(res.stray)})")


def load_scdm():
    """Panel (d), top: wannier90 run time, full vs SCDM-only."""
    out = []
    for name, d_full, d_scdm in SCDM_PAIRS:
        full = parse_wout(d_full / "wannier90.wout")
        scdm = parse_wout(d_scdm / "wannier90.wout")
        same = (d_full / "POSCAR").read_bytes() == (d_scdm / "POSCAR").read_bytes()
        assert same, f"{name}: the two runs use different cells"
        assert all(full[k] == scdm[k] for k in ("nwf", "nb", "nk")), \
            f"{name}: the two runs differ in problem size"
        assert scdm["n_iter"] == [0, 0], f"{name}: not an SCDM-only run"
        assert "EOS_SCDM_ONLY = .TRUE." in (d_scdm / "INCAR").read_text()
        assert "EOS_SCDM_ONLY" not in (d_full / "INCAR").read_text()
        out.append({"label": formula(name), "full": full["t_tot"],
                    "scdm": scdm["t_tot"], "documented": False})
        prov("d", f"{name}: wannier90 Total Execution Time, full | "
             "SCDM-only (ratio)",
             f"{full['t_tot']:.3f} | {scdm['t_tot']:.3f} s  "
             f"({full['t_tot'] / scdm['t_tot']:.0f}x)",
             f"{rel(d_full)}/wannier90.wout | {rel(d_scdm)}/wannier90.wout")
        prov("d", f"{name}: WFs / bands / k-points; iterations set in the "
             "SCDM-only run (loc, dis)",
             f"{full['nwf']} / {full['nb']} / {full['nk']}; "
             f"{scdm['n_iter'][0]}, {scdm['n_iter'][1]}",
             f"{rel(d_scdm)}/wannier90.wout")
        prov("d", f"{name}: element means, full | SCDM-only (not plotted)",
             f"{os_summary(d_full)} | {os_summary(d_scdm)}",
             "pyleos assign_case, defaults, on the two EOSCAR files")
    line = doc_line(RESULTS_MD, DOCUMENTED_SCDM_ROW)
    prov("d", "Li2MnO3 frame, full | SCDM-only: documented only "
         + ("(drawn with open markers)" if SHOW_DOCUMENTED_SCDM
            else "(not plotted)"),
         f"{DOCUMENTED_SCDM['full']:g} | {DOCUMENTED_SCDM['scdm']:g} s",
         f"{rel(RESULTS_MD)}:{line} [hand-copied, no data file]")
    if SHOW_DOCUMENTED_SCDM:
        if line is None:
            raise SystemExit("docs/RESULTS.md no longer quotes 554 -> 8.3 s; "
                             "update DOCUMENTED_SCDM")
        out.append({"label": formula(DOCUMENTED_SCDM["label"]) + "*",
                    "full": DOCUMENTED_SCDM["full"],
                    "scdm": DOCUMENTED_SCDM["scdm"], "documented": True})
    return out


def load_warm_start():
    """Panel (d), bottom: loop time of the analysed MD steps, on vs off."""
    incar = {k: (WS_DIR / f"warmstart_INCAR.{k}").read_text().splitlines()
             for k in ("on", "off")}
    diff = [(a, b) for a, b in zip(incar["on"], incar["off"]) if a != b]
    assert diff == [("EOS_WARM_START = .TRUE.", "EOS_WARM_START = .FALSE.")], \
        "the two warm-start runs differ in more than EOS_WARM_START"
    text = "\n".join(incar["on"])
    every = int(re.search(r"EOS_EVERY_N\s*=\s*(\d+)", text)[1])
    system = re.search(r"SYSTEM\s*=\s*(.+)", text)[1].strip()
    magmom = re.search(r"MAGMOM\s*=\s*(.+)", text)[1]
    n_atoms = sum(int(n) for n in re.findall(r"(\d+)\*", magmom))
    prov("d", "warm-start runs: system; atoms in the cell (sum of the "
         "MAGMOM multiplicities)", f"{system}; {n_atoms}",
         rel(WS_DIR / "warmstart_INCAR.on"))
    out = {}
    for k in ("off", "on"):
        f = WS_DIR / f"warmstart_loop_{k}.txt"
        t = [float(x) for x in f.read_text().split()]
        analysed = [v for i, v in enumerate(t, 1) if i % every == 0]
        plain = [v for i, v in enumerate(t, 1) if i % every and i > 1]
        out[k] = {"analysed": analysed, "mean": statistics.mean(analysed)}
        prov("d", f"warm start {k}: analysed steps (every {every}th of "
             f"{len(t)}), loop time (s)",
             " ".join(f"{v:.1f}" for v in analysed), rel(f))
        prov("d", f"warm start {k}: mean +- s.d. of the {len(analysed)} "
             "analysed steps; sum of all steps",
             f"{out[k]['mean']:.1f} +- {statistics.stdev(analysed):.1f} s; "
             f"{sum(t):.1f} s", rel(f))
        prov("d", f"warm start {k}: other steps (first step excluded), "
             "mean / median (not plotted)",
             f"{statistics.mean(plain):.1f} / {statistics.median(plain):.1f} s",
             rel(f))
    return out


def check_results_md():
    for pattern, what in RESULTS_CHECKS:
        line = doc_line(RESULTS_MD, pattern)
        prov("doc", what, "found" if line else "NOT FOUND",
             f"{rel(RESULTS_MD)}:{line}" if line else rel(RESULTS_MD))


# ================================================================ drawing ==
PAD_CM = 0.02 / CM                     # pad_inches of style.save, in cm
PAGE_W = style.DOUBLE / CM             # 17.1 cm, saved page
PAGE_H = 11.4                          # cm, saved page
FS = 7.0                               # text inside the panels
# All positions below are in cm on the saved page, origin top-left.


def fx(x):
    return (x - PAD_CM) / (PAGE_W - 2 * PAD_CM)


def fy(y):
    return 1 - (y - PAD_CM) / (PAGE_H - 2 * PAD_CM)


def ax_cm(fig, left, top, width, height, **kw):
    return fig.add_axes([fx(left), fy(top + height),
                         width / (PAGE_W - 2 * PAD_CM),
                         height / (PAGE_H - 2 * PAD_CM)], **kw)


def label_at(page, letter, x, y):
    """style.panel_label with its right-bottom corner at page position x, y."""
    style.panel_label(page, letter, dx=fx(x), dy=fy(y))


# ------------------------------------------------------------- panel (a) --
HDR_Y, HDR_H = 0.45, 0.75
STAGES = [            # key, left edge, width, label
    ("scf", 3.05, 0.90, "SCF"),
    ("ovl", 4.25, 1.30, "Overlaps\n" + r"$M^{(k,b)}$"),
    ("scdm", 5.85, 1.55, "SCDM start\n" + r"$A^{(k)}$"),
    ("dis", 7.70, 2.05, "Disentanglement"),
    ("loc", 10.05, 1.95, "Localisation"),
    ("asg", 12.45, 1.45, "Assignment"),
]
SX = {k: (x, w) for k, x, w, _ in STAGES}
ROW_H = 0.38
ROW_Y = {"ref1": 1.52, "ref2": 1.97, "M1": 2.60, "M2": 3.05, "M3": 3.50,
         "M4": 3.95, "M5": 4.40, "M6": 4.85}
TAG_X, NAME_X = 0.30, 0.87
INSET_X = 14.40


def bar(page, row, first, last, text, kind, dashed=False):
    """Outlined rectangle under the stage columns first..last."""
    x0, x1 = SX[first][0], SX[last][0] + SX[last][1]
    y = ROW_Y[row]
    colour = COL_REF if kind == "ref" else COL_MOD
    page.add_patch(Rectangle(
        (x0, y), x1 - x0, ROW_H, fc="none", ec=colour,
        lw=0.6 if kind == "ref" else 0.8,
        ls=(0, (3, 2)) if dashed else "-", joinstyle="miter"))
    page.text((x0 + x1) / 2, y + ROW_H / 2 - 0.01, text, ha="center",
              va="center_baseline", fontsize=FS)


def arrow(page, p0, p1, colour=C["ink"], lw=0.6, **kw):
    page.add_patch(FancyArrowPatch(
        p0, p1, arrowstyle="-|>,head_length=3,head_width=1.5", color=colour,
        lw=lw, shrinkA=0, shrinkB=0, mutation_scale=1, **kw))


def draw_a(page):
    mid = HDR_Y + HDR_H / 2
    # stage chain
    for i, (key, x, w, text) in enumerate(STAGES):
        page.add_patch(Rectangle((x, HDR_Y), w, HDR_H, fc="none",
                                 ec=C["ink"], lw=0.6, joinstyle="miter"))
        page.text(x + w / 2, mid, text, ha="center", va="center",
                  fontsize=FS, linespacing=1.15)
        if i:
            px, pw = STAGES[i - 1][1:3]
            arrow(page, (px + pw, mid), (x, mid))
    lx = SX["loc"][0] + SX["loc"][1]
    page.text((lx + SX["asg"][0]) / 2 - 0.03, mid - 0.10, r"$r_n$",
              ha="center", va="bottom", fontsize=FS)
    # row labels
    for row, t, name in (("ref1", None, "Reference"),
                         ("M1", "M1", "Offline wannier90"),
                         ("M2", "M2", "Threaded MKL"),
                         ("M3", "M3", "Spin-resolved " + NWF),
                         ("M4", "M4", "gpuw90"),
                         ("M5", "M5", "SCDM-only"),
                         ("M6", "M6", "Warm start")):
        y = ROW_Y[row] + ROW_H / 2
        if t:
            page.text(TAG_X, y, t, ha="left", va="center_baseline",
                      fontsize=FS, fontweight="bold")
        page.text(NAME_X if t else TAG_X, y, name, ha="left",
                  va="center_baseline", fontsize=FS)
    # where each stage runs
    bar(page, "ref1", "scf", "asg", "GPU job (4 GPUs)", "ref")
    bar(page, "ref2", "dis", "loc", "wannier90 library, 1 process", "ref")
    bar(page, "M1", "scf", "scdm", r"GPU job: export of $M$, $A$", "mod")
    bar(page, "M1", "dis", "loc", r"CPU node, MPI over $k$-points", "mod")
    bar(page, "M1", "asg", "asg", "CPU", "mod")
    bar(page, "M2", "dis", "loc", "MPI processes × threads", "mod")
    bar(page, "M3", "scdm", "loc", "fewer WFs in the minority-spin channel",
        "mod")
    bar(page, "M4", "dis", "loc", r"batched over $k$ and $b$; CPU or GPU",
        "mod")
    bar(page, "M5", "dis", "loc", "0 iterations (SCDM centres)", "mod",
        dashed=True)
    bar(page, "M6", "loc", "loc", r"$r_n$ of frame $t-1$", "mod")
    draw_inset(page)
    label_at(page, "a", 0.50, HDR_Y + 0.33)


def draw_inset(page):
    """gpuw90 (M4): a loop over (k, b) becomes one batched operation."""
    x0, s = INSET_X, 0.74                 # left edge, side of a matrix
    xm = x0 + 0.08                        # left edge of the matrices
    # wannier90: one matrix per operation, inside a loop
    y_top = 0.92
    page.text(x0, 0.60, "wannier90", ha="left", va="center_baseline",
              fontsize=FS)
    page.add_patch(Rectangle((xm, y_top), s, s, fc="white", ec=COL_REF,
                             lw=0.6, joinstyle="miter"))
    page.text(xm + s / 2, y_top + s / 2, r"$M^{(k,b)}$", ha="center",
              va="center", fontsize=FS)
    page.add_patch(FancyArrowPatch(
        (xm + s + 0.10, y_top + 0.10), (xm + s + 0.10, y_top + s - 0.10),
        connectionstyle="arc3,rad=-1.1",
        arrowstyle="-|>,head_length=2.6,head_width=1.3", color=COL_REF, lw=0.6,
        shrinkA=0, shrinkB=0, mutation_scale=1))
    page.text(xm + s + 0.56, y_top + s / 2, r"$N_k N_b$" + "\noperations",
              ha="left", va="center", fontsize=FS, linespacing=1.2)
    # arrow between the two, with the tag
    ya = y_top + s + 0.22
    arrow(page, (xm + s / 2, ya), (xm + s / 2, ya + 0.48), colour=COL_MOD,
          lw=0.8)
    page.text(xm + s / 2 + 0.16, ya + 0.22, "M4", ha="left",
              va="center_baseline", fontsize=FS, fontweight="bold")
    # gpuw90: the same matrices stacked into one tensor
    yl = ya + 0.48 + 0.30
    page.text(x0, yl, "gpuw90", ha="left", va="center_baseline", fontsize=FS)
    n, off = 4, 0.085
    yb = yl + 0.30 + (n - 1) * off
    for i in reversed(range(n)):
        page.add_patch(Rectangle((xm + i * off, yb - i * off), s, s,
                                 fc="white", ec=COL_MOD, lw=0.8,
                                 joinstyle="miter"))
    page.text(xm + s / 2, yb + s / 2, r"$M$", ha="center", va="center",
              fontsize=FS)
    page.text(xm + s + (n - 1) * off + 0.22, yb + s / 2 - (n - 1) * off / 2,
              "1 batched\noperation", ha="left", va="center", fontsize=FS,
              linespacing=1.2)
    page.text(x0, yb + s + 0.24,
              r"$[N_k,\ N_b,\ N_{\mathrm{WF}},\ N_{\mathrm{WF}}]$",
              ha="left", va="center_baseline", fontsize=FS)


# ---------------------------------------------------------- panels (b-d) --
P_TOP, P_H = 6.25, 4.30                # common top edge and height of b-d
B_TAG_X, B_CFG_X = 0.30, 1.58
B_AX = (3.95, 2.15)                    # axes left edge, width
B_VAL_X, B_FAC_X = 7.05, 7.75          # right edges of the two text columns
C_AX = (9.55, 2.70)
D_AX = (13.90, 2.75)
LEG_Y = P_TOP - 0.10                   # baseline of the legends above c, d


def mark(kind, **kw):
    if kind == "ref":
        style_ = dict(marker=MK_REF, ms=4.0, mfc=COL_REF, mec=COL_REF)
    else:
        style_ = dict(marker=MK_MOD, ms=4.4, mfc=COL_MOD, mec=COL_MOD)
    style_.update(mew=0.8, ls="none", zorder=4)
    style_.update(kw)
    return style_


def plain_y(ax):
    ax.tick_params(axis="y", which="both", left=False, right=False)
    ax.yaxis.set_minor_locator(NullLocator())


def log_x(ax, lo, hi):
    ax.set_xscale("log")
    ax.set_xlim(lo, hi)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())


def legend_row(page, handles, labels, x, **kw):
    """One-line legend whose lower-left corner is at page position x, LEG_Y."""
    opts = dict(loc="lower left", bbox_to_anchor=(x, LEG_Y),
                bbox_transform=page.transData, ncol=len(handles),
                columnspacing=1.2, borderaxespad=0, borderpad=0,
                handletextpad=0.5)
    opts.update(kw)
    page.add_artist(page.legend(handles, labels, **opts))


def draw_b(fig, page, rows):
    ax = ax_cm(fig, B_AX[0], P_TOP, B_AX[1], P_H)
    n = len(rows)
    y_lo, y_hi = -0.75, 2 * n + 2.25
    ax.set_ylim(y_hi, y_lo)

    def y_cm(y):
        return P_TOP + (y - y_lo) / (y_hi - y_lo) * P_H

    base = rows[0][3]
    for g, (stage, name) in enumerate((("dis", "Disentanglement"),
                                       ("wan", "Localisation"))):
        y0 = g * (n + 1.5)
        page.text(B_TAG_X, y_cm(y0), name, ha="left", va="center_baseline",
                  fontsize=FS)
        ax.plot([base[stage]] * 2, [y0 + 0.55, y0 + n + 0.45],
                color=C["muted"], lw=0.4, zorder=1)
        for i, (label, cfg, kind, t, _src) in enumerate(rows):
            y = y0 + 1 + i
            if kind == "mod":
                ax.plot([t[stage], base[stage]], [y, y], color=COL_MOD,
                        lw=0.8, zorder=2, solid_capstyle="butt")
            ax.plot([t[stage]], [y], **mark(kind))
            page.text(B_TAG_X, y_cm(y), label, ha="left",
                      va="center_baseline", fontsize=FS)
            page.text(B_CFG_X, y_cm(y), cfg, ha="left", va="center_baseline",
                      fontsize=FS)
            page.text(B_VAL_X, y_cm(y), f"{t[stage]:.3f} s", ha="right",
                      va="center_baseline", fontsize=FS)
            page.text(B_FAC_X, y_cm(y),
                      "1×" if kind == "ref" else f"{base[stage] / t[stage]:.1f}×",
                      ha="right", va="center_baseline", fontsize=FS)
    log_x(ax, 0.025, 1.6)
    ax.set_yticks([])
    plain_y(ax)
    ax.set_xlabel("Wall time per iteration (s)")
    label_at(page, "b", B_TAG_X + 0.20, P_TOP - 0.12)
    return ax


def draw_c(fig, page, rows):
    ax = ax_cm(fig, C_AX[0], P_TOP, C_AX[1], P_H)
    dy, h = 0.19, 0.32
    for i, (label, kind, t_w90, t_gpu) in enumerate(rows):
        colour = COL_REF if kind == "ref" else COL_MOD
        ax.barh(i - dy, t_w90, height=h, color=colour, edgecolor=colour,
                linewidth=0.8, zorder=3)
        ax.barh(i + dy, t_gpu, height=h, color="white", edgecolor=colour,
                linewidth=0.8, zorder=3)
        for y, v in ((i - dy, t_w90), (i + dy, t_gpu)):
            ax.annotate(f"{v:.1f}", (v, y), xytext=(2.5, 0),
                        textcoords="offset points", ha="left",
                        va="center_baseline", fontsize=FS)
    ax.set_ylim(len(rows) - 0.42, -0.58)
    ax.set_xlim(0, 1150)
    ax.xaxis.set_major_locator(MultipleLocator(300))
    ax.xaxis.set_minor_locator(MultipleLocator(100))
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows])
    ax.tick_params(axis="y", pad=4)
    plain_y(ax)
    ax.set_xlabel("Wall time per frame (s)")
    handles = [Rectangle((0, 0), 1, 1, fc=COL_REF, ec=COL_REF, lw=0.8),
               Rectangle((0, 0), 1, 1, fc="white", ec=COL_REF, lw=0.8)]
    legend_row(page, handles, ["Wannier step", "GPU job"], C_AX[0] - 0.85,
               handlelength=1.1, handleheight=0.75)
    label_at(page, "c", C_AX[0] - 1.00, P_TOP - 0.12)
    return ax


def draw_d(fig, page, scdm, ws):
    n = len(scdm)
    h_top = 0.44 + 0.44 * n
    gap = 1.12                           # tick labels + x label of the top axes
    # --- SCDM-only (M5)
    ax = ax_cm(fig, D_AX[0], P_TOP, D_AX[1], h_top)
    for i, r in enumerate(scdm):
        open_ = dict(mfc="white") if r["documented"] else {}
        ax.plot([r["scdm"], r["full"]], [i, i], color=COL_MOD, lw=0.8,
                ls=(0, (3, 2)) if r["documented"] else "-", zorder=2,
                solid_capstyle="butt")
        ax.plot([r["full"]], [i], **mark("ref", **open_))
        ax.plot([r["scdm"]], [i], **mark("mod", **open_))
        ax.annotate(f"{r['full']:.3g}", (r["full"], i), xytext=(4, 0),
                    textcoords="offset points", ha="left",
                    va="center_baseline", fontsize=FS)
        ax.annotate(f"{r['scdm']:.3f}" if r["scdm"] < 1 else f"{r['scdm']:.3g}",
                    (r["scdm"], i), xytext=(-4, 0),
                    textcoords="offset points", ha="right",
                    va="center_baseline", fontsize=FS)
        ax.annotate(f"{r['full'] / r['scdm']:.0f}×",
                    (np.sqrt(r["full"] * r["scdm"]), i), xytext=(0, 2.5),
                    textcoords="offset points", ha="center", va="bottom",
                    fontsize=FS)
    if any(r["documented"] for r in scdm):
        log_x(ax, 5e-4, 5e4)
    else:
        log_x(ax, 1e-3, 1e3)
    ax.xaxis.set_major_locator(FixedLocator([0.01, 0.1, 1, 10, 100, 1000]))
    ax.xaxis.set_major_formatter(FuncFormatter(
        lambda v, _: f"{v:g}" if round(np.log10(v)) % 2 == 0 else ""))
    ax.set_ylim(n - 0.42, -0.80)
    ax.set_yticks(range(n))
    ax.set_yticklabels([r["label"] for r in scdm])
    ax.tick_params(axis="y", pad=4)
    plain_y(ax)
    ax.set_xlabel("wannier90 run time (s)")
    handles = [Line2D([], [], **mark("ref")), Line2D([], [], **mark("mod"))]
    legend_row(page, handles, ["Reference", tag("M5")], D_AX[0] - 0.85,
               handlelength=0.9, handletextpad=0.2)
    label_at(page, "d", D_AX[0] - 1.00, P_TOP - 0.12)
    # --- warm start (M6)
    top2 = P_TOP + h_top + gap
    h2 = P_TOP + P_H - top2
    ax2 = ax_cm(fig, D_AX[0], top2, D_AX[1], h2)
    y_lo, y_hi = -0.62, 1.62
    x_lo, x_hi, ms = 510, 610, 3.4
    # one row per setting; dots closer than one marker width are displaced
    # sideways (0, up, down, ...) so that all eight steps stay visible
    size = 1.2 * (ms + 0.8) / 72 / CM          # marker + clearance, cm
    min_gap = size / D_AX[1] * (x_hi - x_lo)
    dy = max(0.23, size / h2 * (y_hi - y_lo))
    for i, (key, kind) in enumerate((("off", "ref"), ("on", "mod"))):
        last, steps = {}, sorted(ws[key]["analysed"])
        for v in steps:
            k = next(k for k in (0, -1, 1, -2, 2)
                     if k not in last or v - last[k] >= min_gap)
            last[k] = v
            ax2.plot([v], [i + dy * k], **mark(kind, ms=ms, mfc="white"))
        m = ws[key]["mean"]
        ax2.plot([m, m], [i - 0.36, i + 0.36], color=C["ink"], lw=1.0,
                 zorder=5, solid_capstyle="butt")
        # the mean is labelled on the side with the wider gap to the dots
        right = (min(v for v in steps if v > m) - m
                 > m - max(v for v in steps if v < m))
        ax2.annotate(f"{m:.0f}", (m, i), xytext=(2.5 if right else -2.5, 0),
                     textcoords="offset points", fontsize=FS,
                     ha="left" if right else "right", va="center_baseline")
    ax2.set_ylim(y_hi, y_lo)
    ax2.set_xlim(x_lo, x_hi)
    ax2.xaxis.set_major_locator(MultipleLocator(25))
    ax2.xaxis.set_minor_locator(MultipleLocator(5))
    ax2.set_yticks([0, 1])
    ax2.set_yticklabels(["Reference", tag("M6")])
    ax2.tick_params(axis="y", pad=4)
    plain_y(ax2)
    ax2.set_xlabel("Time per analysed step (s)")
    return ax, ax2


def print_provenance():
    w = [max(len(r[i]) for r in PROV) for i in range(3)]
    w[1], w[2] = min(w[1], 84), min(w[2], 58)
    print("\nPROVENANCE (panel | quantity | value | source)")
    print("-" * 170)
    for p, q, v, s in PROV:
        print(f"({p + ')':<4} {q:<{w[1]}}  {v:<{w[2]}}  {s}")
    print("-" * 170)


def main():
    iters = load_iteration_times()
    pipeline = load_pipeline()
    scdm = load_scdm()
    ws = load_warm_start()
    check_results_md()

    pad = PAD_CM * CM
    fig = plt.figure(figsize=(style.DOUBLE - 2 * pad, PAGE_H * CM - 2 * pad))
    # invisible frame: pins the tight bounding box to the full canvas, so the
    # saved page is exactly one double column (17.1 cm) wide
    fig.add_artist(Rectangle((0, 0), 1, 1, transform=fig.transFigure,
                             fill=False, edgecolor="none", lw=0))
    # full-canvas axes in page centimetres: schematic and free text
    page = fig.add_axes([0, 0, 1, 1])
    page.set_xlim(PAD_CM, PAGE_W - PAD_CM)
    page.set_ylim(PAGE_H - PAD_CM, PAD_CM)
    page.axis("off")

    draw_a(page)
    draw_b(fig, page, iters)
    draw_c(fig, page, pipeline)
    draw_d(fig, page, scdm, ws)
    style.save(fig, "fig3_gpuw90")
    print_provenance()
    print(f"wrote {style.HERE / 'fig3_gpuw90.pdf'} and .png "
          f"(page {PAGE_W:.1f} x {PAGE_H:.1f} cm)")


if __name__ == "__main__":
    main()
