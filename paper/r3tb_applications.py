#!/usr/bin/env python3
"""r3-T-B on the application cases: the test table behind r3tb_fig3.

    python3 paper/r3tb_applications.py [--recompute]

(from the repo root; `python3 r3tb_applications.py` inside paper/ works too).

The frozen plugin pyleos/methods/frozen/r3-T-B-cap-repair.py is run on every
application case whose EOSCAR is in the repo, exactly as
`python3 -m pyleos bench clean` calls a plugin:

    bench.static_check(source); bench.install_trackers()
    plugin = bench.load_plugin(METHOD)
    plugin.assign(bench.scrub(load_case(run_dir)))

(the 13 held-out homonuclear cases go through bench.score_case, as
`python3 -m pyleos bench test --set homonuclear` does).  Next to it:

    flow    pyleos.transport.assign_case(load_case(run_dir)), beta = 10, Allen
            (methane: Co4 total from assign_case(fragments=[Co atoms]), as in
            paper/fig5_case_studies.py)
    greedy  OS(int) column of the EOSCAR (per-WF greedy rule)
    stray   len(assign_case(load_case(run_dir)).stray)

Cases (cases/INDEX.tsv): the Li2O2 delithiation series and every other
registered run of li16o16 / li15o16 / li14o16 (role app-li2o2), methane on
Co4/CeO2 (role app-methane, plus the v31 reruns that paper/fig5 draws),
CeO2 with and without an O vacancy (app-ceo2), the 30 runs of the
Comas-Vila set (comas; eight of them are the MH4 hydrides), the 13
held-out homonuclear cases (test-homonuclear) and Fe3O4 (app-fe3o4).  The
cathode MD runs (app-cathode) have no EOSCAR in the repo and are listed as
not tested.

Output
    results/r3tb_applications.tsv            one row per case and quantity
    paper/_cache/r3tb/applications.json      per-atom results (read by
                                             r3tb_fig3_case_studies.py)
    paper/_cache/r3tb/applications.run.log   what this script printed

Nothing is typed by hand except the literature values LIT_MH4 (same table as
paper/fig5_case_studies.py) and LIT_COMAS, both transcribed from
docs/comas_vila_2025_notes.md (N<k> = line k of that file).
"""
import argparse
import csv
import hashlib
import io
import itertools
import json
import sys
import tempfile
import time
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from pyleos import bench                                        # noqa: E402
from pyleos.parsers import (CaseData, Structure, WannierChannel,  # noqa: E402
                            load_case, parse_eoscar)
from pyleos.transport import assign_case                        # noqa: E402

METHOD = REPO / "pyleos" / "methods" / "frozen" / "r3-T-B-cap-repair.py"
CASES = REPO / "cases"
OUT_TSV = REPO / "results" / "r3tb_applications.tsv"
CACHE_DIR = HERE / "_cache" / "r3tb"
CACHE = CACHE_DIR / "applications.json"
LOG = CACHE_DIR / "applications.run.log"
CACHE_VERSION = 4           # bump when compute() or the cached layout changes
COMMAND = "python3 paper/r3tb_applications.py"
SOURCES = [METHOD] + [
    REPO / "pyleos" / f for f in ("parsers.py", "transport.py", "capacity.py",
                                  "bench.py", "methods/baseline.py",
                                  "sets/test_set.json")]
STUDY_ORDER = ["delithiation", "polaron", "delithiation-other", "li15o16-md",
               "methane", "ceo2", "mh4", "comas-other", "homonuclear",
               "fe3o4", "cathode-md"]

# ---------------------------------------------------------------- inputs --
MH4 = [("mgh4", "Mg"), ("cah4", "Ca"), ("srh4", "Sr"), ("sch4", "Sc"),
       ("yh4", "Y"), ("zrh4", "Zr"), ("ceh4", "Ce"), ("thh4", "Th")]
H2_CUT = 1.0                # two H atoms closer than this (A) are an H2 unit
# Reruns of the methane images that paper/fig5_case_studies.py draws as
# filled symbols (INDEX role "diag"); tested next to the app-methane runs.
METHANE_EXTRA = [("ce18o36co4ch4_is", "v31-mag6"),
                 ("ce18o36co4ch4_is", "v31-mag10"),
                 ("ce18o36co4ch4_ts", "v31-mag6"),
                 ("ce18o36co4ch4_ts", "v31-mag10"),
                 ("ce18o36co4ch4_fs", "v31-mag10"),
                 ("ce18o36co4ch4_fs", "v31-4gpu")]
METHANE_SETS = {          # image -> run, as CH4_SETS of fig5_case_studies.py
    "legacy (v3)": {"is": "v3", "ts": "v3", "fs": "v3"},
    "reruns (v31)": {"is": "v31-mag6", "ts": "v31-mag6", "fs": "v31-mag10"},
}
PAIR_CUTOFF = 2.4           # A: longest O-O distance still paired (fig. 4)
MD_RUN = ("li15o16", "v3-md2")      # trajectory EOSCAR, 20 analysed steps

# =============================================================================
# LITERATURE VALUES (hand-typed).  G. Comas-Vila et al., J. Chem. Theory
# Comput. 2025, 21, 7075-7086, as transcribed in
# docs/comas_vila_2025_notes.md (N<k> = line k of that file).
# LIT_MH4 is the table of paper/fig5_case_studies.py:
#   case: (metal, H_b, H_a per atom, notes line, page)
# =============================================================================
LIT_MH4 = {
    "mgh4": (2, -1, 0.0, 352, "7081"),         # M(+2), Hb(-1), neutral Ha-Ha
    "cah4": (2, -1, 0.0, 352, "7081"),
    "srh4": (2, -1, 0.0, 352, "7081"),
    "sch4": (3, -1, -0.5, 353, "7081-7082"),   # M(+3), Hb(-1), Ha(-1/2)
    "yh4": (3, -1, -0.5, 353, "7081-7082"),
    "zrh4": (4, -1, -1.0, 354, "7082"),        # M4+ (H-)4
    "ceh4": (4, -1, -1.0, 354, "7082"),
    "thh4": (4, -1, -1.0, 354, "7082"),
}
# Other cases of the Comas-Vila set: element-mean oxidation states that the
# notes state (None where the paper gives no element value).
#   case: ({element: mean}, notes line, remark)
LIT_COMAS = {
    "ba3si4": ({"Ba": 2, "Si": -1.5}, 347, "Si4 fragment -6"),
    "cs2pt": ({"Pt": -2}, 348, ""),
    "ybafe2o5_a": ({"Fe": 2.5}, 351, "one phase Fe +3/+2, the other +2.5"),
    "ybafe2o5_b": ({"Fe": 2.5}, 351, "one phase Fe +3/+2, the other +2.5"),
    "xeo_a": ({"Xe": 2}, 355, "100 GPa: +2; 200 GPa: Xe1 +4, Xe2 0"),
    "xeo_b": ({"Xe": 2}, 355, "100 GPa: +2; 200 GPa: Xe1 +4, Xe2 0"),
    "xeo2_a": ({"Xe": 4}, 356, ""),
    "xeo2_b": ({"Xe": 4}, 356, ""),
    "xeo3_a": ({"Xe": 4}, 357, "by atoms Xe +4 (O1 -2, O2 0); "
                               "O2 as fragment: +5"),
    "xeo3_b": ({"Xe": 5}, 358, "Xe1 +6, Xe2 +4; called ambiguous"),
    "xef2": ({"Xe": 2}, 359, ""),
    "xef3": ({"Xe": 3}, 359, "Xe1 +2, Xe2 +4"),
    "xef4": ({"Xe": 4}, 359, ""),
    "mgxe": ({"Mg": 2, "Xe": -2}, 360, ""),
    "hgf4_a": ({"Hg": 4}, 361, ""),
    "hgf4_b": ({"Hg": 4}, 361, ""),
    "hgf4_c": ({"Hg": 4}, 361, ""),
    "csh3": ({"H": -1 / 3}, 362, "H3 fragment -1"),
    "csf3": ({"F": -1 / 3}, 363, "F3 fragment -1"),
    "csf5": ({"F": -0.2}, 364, "F5 fragment -1"),
    "lin5": ({}, 365, "not decided (R = 50)"),
    "k_hp": ({"K": 1}, 345, "electrons on non-nuclear attractors"),
}


# ---------------------------------------------------------------- helpers --
class Tee(io.TextIOBase):
    """stdout that also keeps a copy for the run log."""

    def __init__(self, stream):
        self.stream, self.buf = stream, io.StringIO()

    def write(self, s):
        self.stream.write(s)
        self.buf.write(s)
        return len(s)

    def flush(self):
        self.stream.flush()


def rel(p):
    return str(Path(p).relative_to(REPO))


def md5(path):
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def read_tsv(path):
    with Path(path).open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def signed(v):
    v = int(v)
    return "0" if v == 0 else f"{v:+d}"


def counts(vals):
    """'-1 x15, 0 x1' for a list of integers."""
    return ", ".join(f"{signed(v)} x{n}"
                     for v, n in sorted(Counter(vals).items())) or "-"


def listed(vals):
    return ",".join(signed(v) for v in vals)


def num(v):
    """'-0.5', '0', '-1' for a literature value per atom."""
    return "0" if v == 0 else f"{v:+g}"


def mmm(vals):
    """min/max/mean of a list of integers."""
    return (f"{signed(min(vals))}/{signed(max(vals))}/"
            f"{sum(vals) / len(vals):+.4f}")


def idx(elements, el):
    return [i for i, e in enumerate(elements) if e == el]


def image_dists(st, i, j, cut):
    """All i-j distances below cut, periodic images included (i = j
    allowed)."""
    df = [st.frac[j][k] - st.frac[i][k] for k in range(3)]
    df = [x - round(x) for x in df]
    out = []
    for n in itertools.product((-1, 0, 1), repeat=3):
        f = [df[k] + n[k] for k in range(3)]
        c = [sum(f[k] * st.lattice[k][x] for k in range(3)) for x in range(3)]
        d = sum(x * x for x in c) ** 0.5
        if 1e-6 < d < cut:
            out.append(d)
    return out


def o2_units(st, ox, cutoff=PAIR_CUTOFF):
    """Pair O atoms by repeatedly taking the shortest remaining O-O distance
    (up to cutoff), as paper/fig4_delithiation.py does.
    Returns (pairs [(i, j, d)], unpaired atoms)."""
    dist = sorted((bench.distance(st, i, j), i, j)
                  for a, i in enumerate(ox) for j in ox[a + 1:])
    free, pairs = set(ox), []
    for d, i, j in dist:
        if d > cutoff:
            break
        if i in free and j in free:
            free -= {i, j}
            pairs.append((i, j, d))
    return pairs, sorted(free)


def h_sites(st, hs, ms):
    """H atoms grouped into classes by their H-H and M-H neighbour distances
    (as paper/fig5_case_studies.py): per class the atoms, the nearest H-H
    distance and the pairs of mutually nearest atoms."""
    sig, near = {}, {}
    for i in hs:
        hh = sorted(d for j in hs for d in image_dists(st, i, j, 2.6))
        mh = sorted(d for j in ms for d in image_dists(st, i, j, 3.0))
        sig[i] = (tuple(round(x, 2) for x in hh),
                  tuple(round(x, 2) for x in mh))
        near[i] = min((min(image_dists(st, i, j, 9.0)), j) for j in hs)
    classes = {}
    for i in hs:
        classes.setdefault(sig[i], []).append(i)
    out = []
    for atoms in classes.values():
        pairs = sorted({tuple(sorted((i, near[i][1]))) for i in atoms
                        if near[near[i][1]][1] == i and near[i][1] != i})
        out.append({"atoms": atoms, "d_hh": near[atoms[0]][0],
                    "pairs": pairs})
    return sorted(out, key=lambda c: c["d_hh"])


def eoscar_header(d, key):
    for ln in (Path(d) / "EOSCAR").read_text().splitlines():
        if ln.startswith("# " + key) and "=" in ln:
            return float(ln.split("=", 1)[1])
    return None


def case_files(d):
    return [p for p in (Path(d) / f for f in
                        ("EOSCAR", "POSCAR", "wannier90.wout",
                         "wannier90.1.wout", "wannier90.2.wout", "XDATCAR"))
            if p.is_file()]


# ------------------------------------------------------------ case list ---
def study_of(row):
    case, run, role = row["case"], row["run"], row["role"]
    if role == "app-li2o2":
        if (case, run) == MD_RUN:
            return "li15o16-md"
        if run.startswith("polaron-"):
            return "polaron"
        return "delithiation" if run == "v31" else "delithiation-other"
    if role == "app-methane" or (case, run) in METHANE_EXTRA:
        return "methane"
    if role == "app-ceo2":
        return "ceo2"
    if role == "comas":
        return "mh4" if case in dict(MH4) else "comas-other"
    if role == "test-homonuclear":
        return "homonuclear"
    if role == "app-fe3o4":
        return "fe3o4"
    if role == "app-cathode":
        return "cathode-md"
    return None


def case_list():
    """[(study, case, run, role, structure id, structure status)]."""
    status = {r["id"]: r["status"] for r in read_tsv(CASES / "STRUCTURES.tsv")}
    out = []
    for r in read_tsv(CASES / "INDEX.tsv"):
        s = study_of(r)
        if s:
            out.append((s, r["case"], r["run"], r["role"], r["structure"],
                        status.get(r["structure"], "?")))
    return out


# ------------------------------------------------------------- compute ----
def diag_summary(plugin, case):
    dg = plugin.diagnose(case)
    return {"mode": dg["mode"], "excess": dg.get("excess", 0),
            "undetermined": bool(dg["undetermined"]), "flags": dg["flags"],
            "n_ties": len(dg.get("ties", [])),
            "electrons": [round(x, 6) for x in dg["populations"]]}


def run_r3tb(plugin, d):
    """The plugin on one run directory, as the clean benchmark calls it."""
    os_vals = [int(v) for v in plugin.assign(bench.scrub(load_case(d)))]
    return os_vals, diag_summary(plugin, bench.scrub(load_case(d)))


def run_others(d, elements):
    """Min-cost flow (defaults), stray WFs and the greedy column."""
    case = load_case(d)
    res = assign_case(case)
    eo = case.eoscar
    out = {"flow": list(res.os_int), "stray": len(res.stray),
           "n_wf": [len(ch.centres) for ch in case.channels],
           "ghost": res.ghost,
           "greedy": list(eo.greedy_os),
           "greedy_ne": [round(x, 6) for x in eo.greedy_ne],
           "zval": dict(eo.zval), "dne": eo.dne,
           "nelect_wann": eo.nelect_wann,
           "dne_warn": eoscar_header(d, "EOS_NELECT_WARN")}
    co = idx(elements, "Co")
    if co:                               # Co4 as one fragment (fig. 5a)
        frag = assign_case(load_case(d), fragments=[co])
        out["flow_frag_os"] = frag.frag_os[0]
    return out


def compute_md(plugin, d):
    """r3-T-B and the flow per analysed step of a trajectory EOSCAR.  The
    step-to-geometry map is not recorded (analysis every second MD step), so
    both possible XDATCAR configurations are tried, as fig4_delithiation.py
    does for the flow."""
    lines = (d / "EOSCAR").read_text().splitlines()
    starts = [i for i, ln in enumerate(lines) if ln.startswith("# STEP")]
    x = (d / "XDATCAR").read_text().splitlines()
    scale = float(x[1])
    lat = [[float(v) * scale for v in x[i].split()] for i in (2, 3, 4)]
    els = [e for e, n in zip(x[5].split(), x[6].split())
           for _ in range(int(n))]
    conf, i = [], 7
    while i < len(x):
        if x[i].startswith("Direct"):
            conf.append([[float(v) for v in x[i + 1 + k].split()[:3]]
                         for k in range(len(els))])
            i += 1 + len(els)
        else:
            i += 1
    ox, li = idx(els, "O"), idx(els, "Li")
    steps = []
    with tempfile.TemporaryDirectory() as tmp:
        for s, (a, b) in enumerate(zip(starts, starts[1:] + [len(lines)])):
            p = Path(tmp) / "EOSCAR"
            p.write_text("\n".join(lines[a:b]) + "\n")
            row = {"step": s + 1}
            for tag, k in (("odd", 2 * s), ("even", 2 * s + 1)):
                def make():
                    eo = parse_eoscar(p)
                    chans = [WannierChannel(
                        [list(w.centre) for w in eo.wf_detail if w.spin == sp],
                        [w.spread for w in eo.wf_detail if w.spin == sp])
                        for sp in (1, 2)]
                    return CaseData(d, Structure(lat, list(els), conf[k]),
                                    chans, eo, 2)
                o = [int(v) for v in plugin.assign(bench.scrub(make()))]
                dg = plugin.diagnose(bench.scrub(make()))
                res = assign_case(make())
                case = make()
                pairs, _ = o2_units(case.structure, ox)
                row[tag] = {
                    "r3tb_o": [o[i] for i in ox],
                    "r3tb_li": [o[i] for i in li],
                    "r3tb_units": [o[i] + o[j] for i, j, _ in pairs],
                    "undetermined": bool(dg["undetermined"]),
                    "flow_o": [res.os_int[i] for i in ox],
                    "flow_li": [res.os_int[i] for i in li],
                    "flow_units": [res.os_int[i] + res.os_int[j]
                                   for i, j, _ in pairs],
                    "greedy_o": [case.eoscar.greedy_os[i] for i in ox],
                    "greedy_li": [case.eoscar.greedy_os[i] for i in li],
                    "stray": len(res.stray),
                    "n_wf": sum(len(c.centres) for c in case.channels),
                    "ghost": res.ghost}
            steps.append(row)
    return {"n_conf": len(conf), "steps": steps}


def compute():
    """Everything that needs pyleos; returns the cached blob."""
    src = METHOD.read_text()
    violations, literals = bench.static_check(src)
    if violations:
        sys.exit(f"{METHOD}: forbidden constructs {violations}")
    bench.install_trackers()
    plugin = bench.load_plugin(METHOD)
    todo = case_list()
    runs, t_r3 = {}, 0.0
    # ---- pass 1: the plugin alone (so that the table trackers see only it)
    for study, case, run, role, sid, sstat in todo:
        d = CASES / case / run
        key = f"{case}/{run}"
        rec = {"study": study, "case": case, "run": run, "role": role,
               "structure": sid, "structure_status": sstat, "dir": rel(d)}
        runs[key] = rec
        if not (d / "EOSCAR").is_file():
            rec["tested"] = False
            rec["why"] = "no EOSCAR in the repo"
            continue
        if (case, run) == MD_RUN:
            rec["tested"] = "md"
            continue
        rec["tested"] = True
        rec["elements"] = list(load_case(d).structure.elements)
        t0 = time.time()
        rec["r3tb"], rec["diag"] = run_r3tb(plugin, d)
        t_r3 += time.time() - t0
        print(f"  r3-T-B {key}", flush=True)
    homo = {}
    for tset, batch, name, ref, spec, reg_dir, advisory in \
            bench.test_cases("homonuclear"):
        row, os_vals = bench.score_case(plugin, load_case(REPO / reg_dir),
                                        name, ref, spec, False)
        homo[name] = {"dir": reg_dir, "ref": ref, "site": spec,
                      "r3tb": {k: row[k] for k in
                               ("status", "dev", "atoms_ok", "atoms_labelled",
                                "wrong", "error")}, "r3tb_os": os_vals}
        key = reg_dir.split("cases/", 1)[1]
        assert runs[key]["r3tb"] == os_vals, key   # same call, same answer
    tables_r3tb = sorted(bench.TABLES_USED)
    # ---- pass 2: flow, stray, greedy
    for key, rec in runs.items():
        if rec["tested"] is True:
            rec.update(run_others(REPO / rec["dir"], rec["elements"]))
            print(f"  flow   {key}", flush=True)
        elif rec["tested"] == "md":
            rec["md"] = compute_md(plugin, REPO / rec["dir"])
            print(f"  md     {key}", flush=True)
    base = bench.load_plugin(bench.DEFAULT_METHOD)
    for tset, batch, name, ref, spec, reg_dir, advisory in \
            bench.test_cases("homonuclear"):
        d = REPO / reg_dir
        row, _ = bench.score_case(base, load_case(d), name, ref, spec, False)
        homo[name]["flow"] = {k: row[k] for k in
                              ("status", "dev", "atoms_ok", "atoms_labelled",
                               "wrong", "error")}
        case = load_case(d)
        labels = bench.site_labels(name, case.structure, ref, spec)
        g = list(case.eoscar.greedy_os)
        wrong = Counter((el, lab, o) for el, lab, o in
                        zip(case.structure.elements, labels, g)
                        if lab is not None and o != lab)
        n_lab = sum(lab is not None for lab in labels)
        homo[name]["labels"] = labels
        homo[name]["greedy"] = {
            "atoms_ok": n_lab - sum(wrong.values()), "atoms_labelled": n_lab,
            "wrong": ",".join(f"{el}{lab:+d}->{o:+d}x{k}" for (el, lab, o), k
                              in sorted(wrong.items()))}
    return {"runs": runs, "homonuclear": homo, "tables_used_r3tb": tables_r3tb,
            "element_literals": literals, "seconds_r3tb": round(t_r3, 1)}


def load(force=False):
    """The cached results; recomputed when an input, the plugin or pyleos
    changed (md5) or with force."""
    inputs = []
    for study, case, run, role, sid, sstat in case_list():
        inputs += case_files(CASES / case / run)
    inputs += [CASES / "INDEX.tsv", CASES / "STRUCTURES.tsv"]
    sig = {rel(p): md5(p) for p in inputs + SOURCES}
    sig["cache_version"] = CACHE_VERSION
    if CACHE.is_file() and not force:
        blob = json.loads(CACHE.read_text())
        if blob.get("inputs") == sig:
            return blob["data"]
    t0 = time.time()
    print("[r3tb_applications] running r3-T-B, flow and greedy ...",
          flush=True)
    data = compute()
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(
        {"generated_by": "paper/r3tb_applications.py", "inputs": sig,
         "data": data}, indent=1) + "\n")
    print(f"[r3tb_applications] {time.time() - t0:.0f} s -> {rel(CACHE)}")
    return data


# ----------------------------------------------------- derived quantities --
def structure_of(rec):
    return load_case(REPO / rec["dir"]).structure


def li2o2_view(rec):
    """O2 units, per-atom counts and means of one Li-O run."""
    st = structure_of(rec)
    el = rec["elements"]
    ox, li = idx(el, "O"), idx(el, "Li")
    pairs, free = o2_units(st, ox)
    out = {"n_li": len(li), "n_o": len(ox), "pairs": pairs, "free": free,
           "ref_o_mean": -len(li) / len(ox)}
    for m in ("r3tb", "flow", "greedy"):
        v = rec[m]
        out[m] = {"o": [v[i] for i in ox], "li": [v[i] for i in li],
                  "units": [v[i] + v[j] for i, j, _ in pairs],
                  "free": [v[i] for i in free]}
    z = rec["zval"]["O"]
    ne = rec["greedy_ne"]
    out["greedy"]["units_ne"] = [2 * z - ne[i] - ne[j] for i, j, _ in pairs]
    return out


def mh4_view(rec, metal):
    """Metal and H values of one MH4 run, H atoms ordered as in fig. 5(b)."""
    st = structure_of(rec)
    el = rec["elements"]
    hs, ms = idx(el, "H"), idx(el, metal)
    sites = h_sites(st, hs, ms)
    z = len(ms)
    d_min = min(s["d_hh"] for s in sites)
    out = {"metal": metal, "z": z, "two_sites": len(sites) == 2,
           "d_min": d_min, "h2": d_min < H2_CUT,
           "n_wf": sum(rec["n_wf"]), "stray": rec["stray"],
           "lost": rec["dne"] > rec["dne_warn"], "dne": rec["dne"],
           "lit": LIT_MH4[rec["case"]]}
    for m in ("r3tb", "flow", "greedy"):
        v = rec[m]
        o = {"m": [v[i] for i in ms], "h_all": [v[i] for i in hs],
             "h_mean": sum(v[i] for i in hs) / len(hs)}
        if out["two_sites"]:
            ha, hb = sites              # short H-H contact = H_a (pairs)
            assert len(ha["pairs"]) == z, (rec["case"], ha)
            o["hb"] = [v[i] for i in hb["atoms"]]
            o["ha_pairs"] = [sorted(v[i] for i in p) for p in ha["pairs"]]
            out["d_ha"], out["d_hb"] = ha["d_hh"], hb["d_hh"]
        else:
            (site,) = sites
            assert z == 1 and len(site["pairs"]) == 2, (rec["case"], site)
            # the two equivalent pairs; lower total first (drawing convention
            # of fig. 5(b), no site assignment)
            o["pairs"] = sorted((sorted(v[i] for i in p)
                                 for p in site["pairs"]),
                                key=lambda q: (sum(q), q))
            out["d_h"] = site["d_hh"]
        out[m] = o
    return out


def methane_view(rec):
    st = structure_of(rec)
    el = rec["elements"]
    co, hs, cs, ce, ox = (idx(el, x) for x in ("Co", "H", "C", "Ce", "O"))
    out = {"d_ch": [bench.distance(st, h, cs[0]) for h in hs],
           "h_atoms": [i + 1 for i in hs], "c_atom": cs[0] + 1}
    for m in ("r3tb", "flow", "greedy"):
        v = rec[m]
        out[m] = {"h": [v[i] for i in hs], "c": v[cs[0]],
                  "co": [v[i] for i in co], "co4": sum(v[i] for i in co),
                  "ce": [v[i] for i in ce], "o": [v[i] for i in ox]}
    out["flow"]["co4"] = rec["flow_frag_os"]       # fragment value (fig. 5a)
    return out


# --------------------------------------------------------------- the table --
COLUMNS = ["study", "case", "run", "role", "n_atoms", "n_wf", "n_stray",
           "quantity", "r3tb", "flow", "greedy", "reference", "r3tb_vs_ref",
           "note"]


class Table:
    def __init__(self):
        self.rows = []

    def add(self, rec, quantity, r3tb, flow="-", greedy="-", reference="-",
            match=None, note=""):
        tested = rec.get("tested") is True
        self.rows.append({
            "study": rec["study"], "case": rec["case"], "run": rec["run"],
            "role": rec["role"],
            "n_atoms": len(rec["elements"]) if tested else "-",
            "n_wf": sum(rec["n_wf"]) if tested else "-",
            "n_stray": rec["stray"] if tested else "-",
            "quantity": quantity, "r3tb": r3tb, "flow": flow,
            "greedy": greedy, "reference": reference,
            "r3tb_vs_ref": "-" if match is None else
            ("equal" if match else "differs"), "note": note})

    def write(self, path, header):
        rows = sorted(self.rows,                    # stable: keeps the order
                      key=lambda r: STUDY_ORDER.index(r["study"]))
        with Path(path).open("w") as fh:
            for ln in header:
                fh.write("# " + ln + "\n")
            fh.write("\t".join(COLUMNS) + "\n")
            for r in rows:
                fh.write("\t".join(str(r[c]) for c in COLUMNS) + "\n")


def element_rows(tab, rec, ref_means, ref_src=""):
    """One row per element: min/max/mean of the per-atom oxidation states."""
    by = defaultdict(list)
    for i, e in enumerate(rec["elements"]):
        by[e].append(i)
    for e, ii in by.items():
        ref = ref_means.get(e)
        mean = sum(rec["r3tb"][i] for i in ii) / len(ii)
        tab.add(rec, f"OS({e}) min/max/mean, {len(ii)} atoms",
                mmm([rec["r3tb"][i] for i in ii]),
                mmm([rec["flow"][i] for i in ii]),
                mmm([rec["greedy"][i] for i in ii]),
                "-" if ref is None else f"mean {ref:+.4f}",
                None if ref is None else abs(mean - ref) < 1e-3,
                ref_src if ref is not None else "")


def general_rows(tab, rec):
    dg = rec["diag"]
    note = (f"structure {rec['structure']}: {rec['structure_status']} "
            f"(cases/STRUCTURES.tsv); WFs beyond NELECT(Wannier): "
            f"{rec['ghost']} electrons; |dNe| {rec['dne']:.1e}")
    tab.add(rec, "r3-T-B diagnose(): mode; undetermined",
            f"{dg['mode']}; {'yes' if dg['undetermined'] else 'no'}",
            "-", "-", "-", None,
            note + (f"; flags: {' | '.join(dg['flags'])}" if dg["flags"]
                    else ""))


def build_table(data):
    runs, tab = data["runs"], Table()
    verdict = {}                        # study -> list of (label, ok, text)

    def judge(study, label, ok, text):
        verdict.setdefault(study, []).append((label, ok, text))

    # ---- Li2O2 delithiation (the cells of fig. 4a) and the other Li-O runs
    for key, rec in runs.items():
        if rec["study"] not in ("delithiation", "delithiation-other",
                                "polaron"):
            continue
        v = li2o2_view(rec)
        general_rows(tab, rec)
        element_rows(tab, rec, {"Li": 1.0, "O": v["ref_o_mean"]},
                     "charge counting, Li +1 in a neutral cell")
        o_mean = {m: sum(v[m]["o"]) / v["n_o"] for m in ("r3tb", "flow",
                                                         "greedy")}
        tab.add(rec, "mean OS(O)", f"{o_mean['r3tb']:+.4f}",
                f"{o_mean['flow']:+.4f}", f"{o_mean['greedy']:+.4f}",
                f"{v['ref_o_mean']:+.4f}",
                abs(o_mean["r3tb"] - v["ref_o_mean"]) < 1e-6,
                "reference: -n(Li)/n(O)")
        tab.add(rec, "OS(Li), all atoms", counts(v["r3tb"]["li"]),
                counts(v["flow"]["li"]), counts(v["greedy"]["li"]),
                f"+1 x{v['n_li']}", set(v["r3tb"]["li"]) == {1})
        tab.add(rec, "O atoms by OS", counts(v["r3tb"]["o"]),
                counts(v["flow"]["o"]), counts(v["greedy"]["o"]))
        holes = v["n_o"] - v["n_li"]
        n_units = len(v["pairs"])
        ref_units = sorted([-2] * (n_units - holes) + [-1] * holes) \
            if not v["free"] and 0 <= holes <= n_units else None
        tab.add(rec, f"charges of the {n_units} O2 units (O-O < "
                f"{PAIR_CUTOFF} A)", counts(v["r3tb"]["units"]),
                counts(v["flow"]["units"]), counts(v["greedy"]["units"]),
                counts(ref_units) if ref_units else "-",
                None if ref_units is None
                else sorted(v["r3tb"]["units"]) == ref_units,
                "reference: one hole per removed Li, each on one O2 unit"
                + (f"; unpaired O atoms {[i + 1 for i in v['free']]}: "
                   f"r3tb {listed(v['r3tb']['free'])}, flow "
                   f"{listed(v['flow']['free'])}" if v["free"] else ""))
        if rec["study"] in ("delithiation", "delithiation-other"):
            ok = (abs(o_mean["r3tb"] - v["ref_o_mean"]) < 1e-6
                  and set(v["r3tb"]["li"]) == {1}
                  and sorted(v["r3tb"]["units"]) == ref_units)
            judge(rec["study"], key, ok,
                  f"[structure {rec['structure_status']}] "
                  f"O mean {o_mean['r3tb']:+.4f} (ref "
                  f"{v['ref_o_mean']:+.4f}), Li {counts(v['r3tb']['li'])}, "
                  f"O atoms {counts(v['r3tb']['o'])}, O2 units "
                  f"{counts(v['r3tb']['units'])}; flow O atoms "
                  f"{counts(v['flow']['o'])}; greedy O atoms "
                  f"{counts(v['greedy']['o'])}; stray {rec['stray']}")
        if rec["study"] == "polaron":
            order = sorted(range(n_units), key=lambda k: v["pairs"][k][2])
            k0 = order[0]
            i, j, d0 = v["pairs"][k0]
            d1 = v["pairs"][order[1]][2]
            short = d0 < 1.45           # a unit of superoxide length exists
            where = lambda m: "; ".join(                     # noqa: E731
                [f"O{v['pairs'][k][0] + 1}-O{v['pairs'][k][1] + 1} "
                 f"{v['pairs'][k][2]:.3f} A {signed(v[m]['units'][k])}"
                 for k in order if v[m]["units"][k] != -2]
                + [f"O{a + 1} unpaired {signed(q)}"
                   for a, q in zip(v["free"], v[m]["free"])]) or "-"
            tab.add(rec, f"charge of the shortest O2 unit (O{i + 1}-O{j + 1},"
                    f" {d0:.3f} A; next {d1:.3f} A)",
                    signed(v["r3tb"]["units"][k0]),
                    signed(v["flow"]["units"][k0]),
                    f"{signed(v['greedy']['units'][k0])} (electron count "
                    f"{v['greedy']['units_ne'][k0]:+.2f})",
                    "-1" if short else "-",
                    v["r3tb"]["units"][k0] == -1 if short else None,
                    "reference: the hole sits on the short-bond (superoxide)"
                    " unit" if short else "no unit below 1.45 A")
            tab.add(rec, "O2 units / unpaired O with charge other than -2",
                    where("r3tb"), where("flow"), where("greedy"))
            judge("polaron", key, v["r3tb"]["units"][k0] == -1 if short
                  else None,
                  f"shortest unit {d0:.3f} A: r3tb "
                  f"{signed(v['r3tb']['units'][k0])}, flow "
                  f"{signed(v['flow']['units'][k0])}; r3tb units not -2: "
                  f"{where('r3tb')}")

    # ---- Li15O16 MD (trajectory EOSCAR)
    for key, rec in runs.items():
        if rec.get("tested") != "md":
            continue
        steps = rec["md"]["steps"]
        fake = dict(rec, tested=False)
        for tag in ("odd", "even"):
            def span(name, f=lambda q: sum(q) / len(q)):
                vals = [f(s[tag][name]) for s in steps]
                return f"{min(vals):+.4f}..{max(vals):+.4f}"

            def unit_span(name):
                u = [q for s in steps for q in s[tag][name]]
                return f"{signed(min(u))}..{signed(max(u))}"
            stray = [s[tag]["stray"] for s in steps]
            note = (f"{len(steps)} steps, {steps[0][tag]['n_wf']} WFs each, "
                    f"stray {min(stray)}..{max(stray)}, surplus e "
                    f"{steps[0][tag]['ghost']}; step-to-geometry map not "
                    f"recorded: XDATCAR configurations {tag}")
            tab.add(fake, f"mean OS(O) over the steps, min..max [{tag}]",
                    span("r3tb_o"), span("flow_o"), span("greedy_o"),
                    "-0.9375", all(abs(sum(s[tag]["r3tb_o"]) / 16 + 0.9375)
                                   < 1e-6 for s in steps), note)
            tab.add(fake, f"OS(Li) over all steps and atoms, min..max [{tag}]",
                    unit_span("r3tb_li"), unit_span("flow_li"),
                    unit_span("greedy_li"), "+1..+1",
                    all(set(s[tag]["r3tb_li"]) == {1} for s in steps), note)
            tab.add(fake, f"charge of an O2 unit, min..max [{tag}]",
                    unit_span("r3tb_units"), unit_span("flow_units"), "-",
                    "-2..-1", None,
                    f"r3-T-B undetermined in "
                    f"{sum(s[tag]['undetermined'] for s in steps)} of "
                    f"{len(steps)} steps")
        ok = all(abs(sum(s[t]["r3tb_o"]) / 16 + 0.9375) < 1e-6
                 and set(s[t]["r3tb_li"]) == {1}
                 for s in steps for t in ("odd", "even"))
        stray = [s["odd"]["stray"] for s in steps]
        units = [q for s in steps for t in ("odd", "even")
                 for q in s[t]["r3tb_units"]]
        judge("li15o16-md", key, ok and set(units) <= {-2, -1},
              f"O mean -0.9375 and Li +1 in all {len(steps)} steps: "
              f"{'yes' if ok else 'no'}; O2-unit charges "
              f"{signed(min(units))}..{signed(max(units))} (ref -2..-1); "
              f"{min(stray)}..{max(stray)} stray WFs of "
              f"{steps[0]['odd']['n_wf']} per step")

    # ---- methane on Co4/CeO2
    views = {}
    d_fs = methane_view(runs["ce18o36co4ch4_fs/v3"])["d_ch"]
    k_transferred = max(range(4), key=lambda q: d_fs[q])  # longest C-H, final
    for key, rec in runs.items():
        if rec["study"] != "methane" or rec.get("tested") is not True:
            continue
        v = methane_view(rec)
        views[key] = v
        image = rec["case"].rsplit("_", 1)[1]
        k = k_transferred
        ht = v["h_atoms"][k]
        ref_h = {"is": 1, "ts": -1, "fs": -1}[image]
        general_rows(tab, rec)
        element_rows(tab, rec, {"C": -4.0, "O": -2.0,
                                "H": (3 + ref_h) / 4},
                     "fig. 5(a): C -4, O -2, three H +1 and H_t")
        tab.add(rec, f"OS of the transferred H (atom {ht}, C-H "
                f"{v['d_ch'][k]:.3f} A)", signed(v["r3tb"]["h"][k]),
                signed(v["flow"]["h"][k]), signed(v["greedy"]["h"][k]),
                signed(ref_h), v["r3tb"]["h"][k] == ref_h,
                "reference: fig. 5(a), +1 before and -1 after C-H cleavage")
        tab.add(rec, f"OS of C (atom {v['c_atom']})", signed(v["r3tb"]["c"]),
                signed(v["flow"]["c"]), signed(v["greedy"]["c"]), "-4",
                v["r3tb"]["c"] == -4, "reference: fig. 5(a)")
        tab.add(rec, "Co4 total OS (per-atom values)",
                f"{signed(v['r3tb']['co4'])} ({listed(v['r3tb']['co'])})",
                f"{signed(v['flow']['co4'])} (fragment; per-atom run "
                f"{listed(v['flow']['co'])})",
                f"{signed(v['greedy']['co4'])} ({listed(v['greedy']['co'])})",
                "-", None, "fig. 5(a) supports only the step: +2 from the "
                "initial to the transition and final state")
        tab.add(rec, "Ce atoms by OS", counts(v["r3tb"]["ce"]),
                counts(v["flow"]["ce"]), counts(v["greedy"]["ce"]), "-", None,
                "WF occupations indicate three Ce(+3) (fig5_caption.md)")
        tab.add(rec, "O atoms by OS", counts(v["r3tb"]["o"]),
                counts(v["flow"]["o"]), counts(v["greedy"]["o"]),
                f"-2 x{len(v['r3tb']['o'])}", set(v["r3tb"]["o"]) == {-2})
        ok = (v["r3tb"]["h"][k] == ref_h and v["r3tb"]["c"] == -4
              and set(v["r3tb"]["o"]) == {-2}
              and set(v["r3tb"]["ce"]) <= {3, 4})
        judge("methane", key, ok,
              f"H_t {signed(v['r3tb']['h'][k])} (ref {signed(ref_h)}), C "
              f"{signed(v['r3tb']['c'])} (ref -4), Co4 "
              f"{signed(v['r3tb']['co4'])}, Ce {counts(v['r3tb']['ce'])}, O "
              f"{counts(v['r3tb']['o'])}; flow: H_t "
              f"{signed(v['flow']['h'][k])}, C {signed(v['flow']['c'])}, Co4 "
              f"{signed(v['flow']['co4'])}; stray {rec['stray']}")
    for name, images in METHANE_SETS.items():
        keys = [f"ce18o36co4ch4_{im}/{images[im]}" for im in ("is", "ts",
                                                              "fs")]
        rec = dict(runs[keys[0]], case="ce18o36co4ch4_is/ts/fs",
                   run="/".join(images[im] for im in ("is", "ts", "fs")),
                   tested=False)
        co4 = {m: [views[q][m]["co4"] for q in keys]
               for m in ("r3tb", "flow", "greedy")}
        tab.add(rec, f"Co4 total OS, initial/transition/final [{name}]",
                "/".join(signed(x) for x in co4["r3tb"]),
                "/".join(signed(x) for x in co4["flow"]),
                "/".join(signed(x) for x in co4["greedy"]),
                "step +2, then constant",
                co4["r3tb"][1] - co4["r3tb"][0] == 2
                and co4["r3tb"][2] == co4["r3tb"][1],
                "reference: fig. 5(a), the step the two data sets agree on")

    # ---- CeO2 with and without an O vacancy
    for key, rec in runs.items():
        if rec["study"] != "ceo2":
            continue
        el = rec["elements"]
        ce, ox = idx(el, "Ce"), idx(el, "O")
        n_vac = 2 * len(ce) - len(ox)
        general_rows(tab, rec)
        element_rows(tab, rec, {"O": -2.0, "Ce": 2 * len(ox) / len(ce)},
                     "charge counting with O -2")
        g = lambda m: [rec[m][i] for i in ce]               # noqa: E731
        tab.add(rec, "OS of " + ", ".join(f"Ce{i + 1}" for i in ce),
                listed(g("r3tb")), listed(g("flow")), listed(g("greedy")),
                "-", None, "the four Ce atoms are equivalent neighbours of "
                "the vacant site (2.343 A)" if n_vac else "")
        tab.add(rec, "number of Ce at +3", g("r3tb").count(3),
                g("flow").count(3), g("greedy").count(3), 2 * n_vac,
                g("r3tb").count(3) == 2 * n_vac
                and set(g("r3tb")) <= {3, 4},
                "reference: two electrons per O vacancy, O at -2")
        tab.add(rec, "O atoms by OS", counts([rec["r3tb"][i] for i in ox]),
                counts([rec["flow"][i] for i in ox]),
                counts([rec["greedy"][i] for i in ox]), f"-2 x{len(ox)}",
                {rec["r3tb"][i] for i in ox} == {-2})
        ok = (g("r3tb").count(3) == 2 * n_vac and set(g("r3tb")) <= {3, 4}
              and {rec["r3tb"][i] for i in ox} == {-2})
        judge("ceo2", key, ok,
              f"Ce {listed(g('r3tb'))} (flow {listed(g('flow'))}, greedy "
              f"{listed(g('greedy'))}), Ce(+3) {g('r3tb').count(3)} (ref "
              f"{2 * n_vac}), O {counts([rec['r3tb'][i] for i in ox])}; "
              f"stray {rec['stray']}")

    # ---- MH4 hydrides against Comas-Vila et al.
    for name, metal in MH4:
        key = f"{name}/v3"
        rec = runs[key]
        v = mh4_view(rec, metal)
        lit = v["lit"]
        lit_mean = (2 * lit[1] + 2 * lit[2]) / 4
        src = f"Comas-Vila et al. 2025, notes N{lit[3]}, p. {lit[4]}"
        general_rows(tab, rec)
        element_rows(tab, rec, {metal: float(lit[0]), "H": lit_mean}, src)
        tab.add(rec, f"OS({metal})", counts(v["r3tb"]["m"]),
                counts(v["flow"]["m"]), counts(v["greedy"]["m"]),
                signed(lit[0]), set(v["r3tb"]["m"]) == {lit[0]}, src)
        tab.add(rec, "mean OS(H)", f"{v['r3tb']['h_mean']:+.4f}",
                f"{v['flow']['h_mean']:+.4f}", f"{v['greedy']['h_mean']:+.4f}",
                f"{lit_mean:+.4f}", abs(v["r3tb"]["h_mean"] - lit_mean) < 1e-6,
                src)
        if v["two_sites"]:
            def sites(m):
                return ("H_b " + listed(v[m]["hb"]) + "; H_a-H_a "
                        + " ".join("(" + listed(p) + ")"
                                   for p in v[m]["ha_pairs"]))
            exact = (set(v["r3tb"]["hb"]) == {lit[1]}
                     and all(set(p) == {lit[2]}
                             for p in v["r3tb"]["ha_pairs"]))
            tab.add(rec, "OS of the H atoms by site", sites("r3tb"),
                    sites("flow"), sites("greedy"),
                    f"H_b {signed(lit[1])}; H_a {num(lit[2])} each", exact,
                    f"I4/mmm cell with sites H_b and H_a; H_a-H_a "
                    f"{v['d_ha']:.3f} A, H_b-H {v['d_hb']:.3f} A; " + src)
        else:
            def sites(m):
                return "pairs " + " ".join("(" + listed(p) + ")"
                                           for p in v[m]["pairs"])
            flat = sorted(x for p in v["r3tb"]["pairs"] for x in p)
            exact = flat == sorted([lit[1]] * 2 + [lit[2]] * 2)
            tab.add(rec, "OS of the H atoms (two equivalent pairs)",
                    sites("r3tb"), sites("flow"), sites("greedy"),
                    f"H_b {signed(lit[1])}; H_a {num(lit[2])} each", exact,
                    f"5-atom cell with one H site (H-H {v['d_h']:.3f} A), "
                    "not the two-site structure of the literature; " + src)
        tab.add(rec, "shortest H-H distance (A)", f"{v['d_min']:.3f}", "-",
                "-", "-", None, "H2 unit" if v["h2"] else "no H2 unit")
        judge("mh4", key,
              set(v["r3tb"]["m"]) == {lit[0]}
              and abs(v["r3tb"]["h_mean"] - lit_mean) < 1e-6,
              f"{metal} {counts(v['r3tb']['m'])} (lit {signed(lit[0])}, "
              f"greedy {counts(v['greedy']['m'])}), H mean "
              f"{v['r3tb']['h_mean']:+.2f} (lit {lit_mean:+.2f}), H atoms "
              f"{sites('r3tb')} [every atom as lit: "
              f"{'yes' if exact else 'no'}]; flow {sites('flow')}; stray "
              f"{rec['stray']}/{v['n_wf']}")

    # ---- the other cases of the Comas-Vila set
    for key, rec in runs.items():
        if rec["study"] != "comas-other":
            continue
        ref, line, remark = LIT_COMAS[rec["case"]]
        general_rows(tab, rec)
        element_rows(tab, rec, {e: float(x) for e, x in ref.items()},
                     f"Comas-Vila et al. 2025, notes N{line}"
                     + (f" ({remark})" if remark else ""))
        el = rec["elements"]
        same = all(abs(sum(rec["r3tb"][i] for i in idx(el, e))
                       - sum(rec["flow"][i] for i in idx(el, e))) < 1e-9
                   for e in set(el))
        lit_ok = all(abs(sum(rec["r3tb"][i] for i in idx(el, e))
                         / len(idx(el, e)) - x) < 1e-3
                     for e, x in ref.items()) if ref else None
        judge("comas-other", key, lit_ok,
              "; ".join(f"{e} {mmm([rec['r3tb'][i] for i in idx(el, e)])}"
                        for e in dict.fromkeys(el))
              + f" | element means {'equal to' if same else 'differ from'} "
              f"the flow | stray {rec['stray']}/{sum(rec['n_wf'])}")

    # ---- Fe3O4 (Verwey cell)
    for key, rec in runs.items():
        if rec["study"] != "fe3o4":
            continue
        el = rec["elements"]
        fe = idx(el, "Fe")
        general_rows(tab, rec)
        element_rows(tab, rec, {"O": -2.0, "Fe": 8 / 3},
                     "charge counting with O -2")
        ref = [2] * (len(fe) // 3) + [3] * (2 * len(fe) // 3)
        g = lambda m: [rec[m][i] for i in fe]               # noqa: E731
        tab.add(rec, "Fe atoms by OS", counts(g("r3tb")), counts(g("flow")),
                counts(g("greedy")), counts(ref), sorted(g("r3tb")) == ref,
                "reference: Fe(+2) : Fe(+3) = 1 : 2")
        judge("fe3o4", key, sorted(g("r3tb")) == ref,
              f"Fe {counts(g('r3tb'))} (ref {counts(ref)}; flow "
              f"{counts(g('flow'))}; greedy {counts(g('greedy'))}); stray "
              f"{rec['stray']}")

    # ---- held-out homonuclear set
    homo = data["homonuclear"]
    tot = {m: [0, 0, 0] for m in ("r3tb", "flow", "greedy")}
    for name, h in homo.items():
        key = h["dir"].split("cases/", 1)[1]
        rec = runs[key]
        general_rows(tab, rec)
        element_rows(tab, rec, {e: float(x) for e, x in h["ref"].items()},
                     "pyleos/sets/test_set.json")
        n = h["r3tb"]["atoms_labelled"]
        for m in tot:
            tot[m][0] += h[m]["atoms_ok"]
            tot[m][1] += h[m]["atoms_labelled"]
            tot[m][2] += h[m]["atoms_ok"] == h[m]["atoms_labelled"]
        tab.add(rec, "atoms correct / labelled",
                f"{h['r3tb']['atoms_ok']}/{n}", f"{h['flow']['atoms_ok']}/{n}",
                f"{h['greedy']['atoms_ok']}/{n}",
                f"{n}/{n}", h["r3tb"]["atoms_ok"] == n,
                "labels: pyleos/sets/test_set.json (bench.site_labels)")
        tab.add(rec, "wrong sites (element, label -> assigned x count)",
                h["r3tb"]["wrong"] or "-", h["flow"]["wrong"] or "-",
                h["greedy"]["wrong"] or "-")
        judge("homonuclear", key, h["r3tb"]["atoms_ok"] == n,
              f"{h['r3tb']['atoms_ok']}/{n} {h['r3tb']['wrong']} (flow "
              f"{h['flow']['atoms_ok']}/{n}, greedy "
              f"{h['greedy']['atoms_ok']}/{n}); stray {rec['stray']}")
    total = {"study": "homonuclear", "case": "TOTAL (13 cases)", "run": "v3",
             "role": "test-homonuclear", "elements": [], "tested": False}
    tab.add(total, "atoms correct / labelled",
            *(f"{tot[m][0]}/{tot[m][1]}" for m in ("r3tb", "flow", "greedy")),
            "131/134 (frozen one-shot run of r3-T-B)",
            (tot["r3tb"][0], tot["r3tb"][1]) == (131, 134),
            "frozen result: results/testset_results_summary.tsv "
            "(r3-T-B 131/134, baseline 71/134)")
    tab.add(total, "cases with every atom correct",
            *(f"{tot[m][2]}/{len(homo)}" for m in ("r3tb", "flow", "greedy")),
            "12/13 (frozen one-shot run of r3-T-B)", tot["r3tb"][2] == 12,
            "frozen result: results/testset_results_summary.tsv "
            "(r3-T-B 12/13, baseline 0/13)")

    # ---- not tested
    for key, rec in runs.items():
        if rec.get("tested") is False:
            tab.add(rec, "not tested", "-", "-", "-", "-", None, rec["why"])
    return tab, verdict, tot


def check_frozen(data, tot):
    """The homonuclear rows against the frozen one-shot result file."""
    frozen = [r for r in read_tsv(REPO / "results" /
                                  "testset_results_cases.tsv")
              if r["method"] == "r3-T-B-cap-repair"
              and r["set"] == "homonuclear"]
    diff = []
    for r in frozen:
        h = data["homonuclear"][r["case"]]["r3tb"]
        mine = (h["status"], str(h["atoms_ok"]), str(h["atoms_labelled"]),
                h["wrong"] or "-")
        theirs = (r["status"], r["atoms_ok"], r["atoms_labelled"], r["wrong"])
        if mine != theirs:
            diff.append((r["case"], mine, theirs))
    return len(frozen), diff


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--recompute", action="store_true",
                    help="ignore paper/_cache/r3tb/applications.json")
    args = ap.parse_args()
    tee = Tee(sys.stdout)
    sys.stdout = tee
    try:
        data = load(args.recompute)
        tab, verdict, tot = build_table(data)
        n_tested = sum(r.get("tested") is True for r in data["runs"].values())
        header = [
            "r3-T-B (pyleos/methods/frozen/r3-T-B-cap-repair.py) on the "
            "application cases",
            f"produced by: {COMMAND}   (run from the repo root; "
            f"{time.strftime('%Y-%m-%d')})",
            "r3tb   = plugin.assign(bench.scrub(load_case("
            "cases/<case>/<run>))) with plugin = bench.load_plugin(...), "
            "as `pyleos bench clean` calls it",
            "flow   = pyleos.transport.assign_case(load_case(...)), beta = 10,"
            " Allen (methane Co4 total: fragments=[Co atoms])",
            "greedy = OS(int) column of the EOSCAR (per-WF greedy rule)",
            "n_stray = len(assign_case(load_case(...)).stray); n_wf = WFs of "
            "all spin channels",
            "values 'a/b/c' in the OS(El) rows are min/max/mean over the atoms"
            " of the element; 'v xN' = N atoms (or units) at v",
            "r3tb_vs_ref compares the r3tb column with the reference column "
            "(OS(El) rows: the mean)",
            f"{n_tested} runs tested; tables read by r3-T-B: "
            f"{data['tables_used_r3tb'] or 'none'}; r3-T-B time "
            f"{data['seconds_r3tb']} s (assign + diagnose)",
        ]
        OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
        tab.write(OUT_TSV, header)
        print(f"\n{len(tab.rows)} rows -> {rel(OUT_TSV)}  ({n_tested} runs "
              f"tested; tables read by r3-T-B: "
              f"{data['tables_used_r3tb'] or 'none'}; element literals in "
              f"the plugin source: {data['element_literals'] or 'none'})")
        for study in sorted(verdict, key=STUDY_ORDER.index):
            items = verdict[study]
            good = sum(ok is True for _, ok, _ in items)
            scored = sum(ok is not None for _, ok, _ in items)
            print(f"\n== {study}: r3-T-B equals the reference in {good} of "
                  f"{scored} runs with a reference ({len(items)} runs)")
            for label, ok, text in items:
                mark = {True: "ok  ", False: "DIFF", None: "    "}[ok]
                print(f"  {mark} {label:<28} {text}")
        n_frozen, diff = check_frozen(data, tot)
        print(f"\n== homonuclear totals: r3-T-B {tot['r3tb'][0]}/"
              f"{tot['r3tb'][1]} atoms, {tot['r3tb'][2]}/13 cases | flow "
              f"{tot['flow'][0]}/{tot['flow'][1]}, {tot['flow'][2]}/13 | "
              f"greedy {tot['greedy'][0]}/{tot['greedy'][1]}, "
              f"{tot['greedy'][2]}/13")
        print(f"   frozen one-shot rows (results/testset_results_cases.tsv, "
              f"{n_frozen} rows): "
              + ("all reproduced" if not diff else f"DIFFER: {diff}"))
    finally:
        sys.stdout = tee.stream
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        LOG.write_text(f"$ {COMMAND}\n" + tee.buf.getvalue())


if __name__ == "__main__":
    main()
