"""Case registry: cases/INDEX.tsv lists every directory under cases/,
cases/STRUCTURES.tsv records where each distinct structure comes from.

INDEX.tsv, one row per run directory cases/<case>/<run>/:
  case       chemical system (directory name)
  run        run label: program version + what differs (v1, v3, v31-wp, ...)
  role       what the run is used for (bench74, comas, test-regression, ...)
  bench_set  74-case benchmark set (v1_mol / v1_oxide / wave1 / wave2) or empty
  ref_key    key of the case in pyleos/sets/reference_os.json, or empty
  origin     path the run had before the cases/ layout (provenance)
  structure  id of the run's structure in STRUCTURES.tsv

STRUCTURES.tsv, one row per distinct structure (runs whose cells match share
a row; the columns are described in cases/README.md). check_registry() keeps
the two tables consistent with the files on disk.
"""
from __future__ import annotations

import csv
import itertools
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
STATUS = {"ok", "suspect", "wrong", "unchecked"}
KIND = {"registered", "database", "literature", "experimental-cif", "hand-built",
        "derived", "generated", "legacy", "copy", "unknown"}
CONFIDENCE = {"verified", "phase-matched", "inferred", "unrecorded"}
SCORED = {"bench74", "comas", "test-regression", "test-homonuclear", "clean-override"}
MIN_DIST = 0.5        # Angstrom; anything closer is a broken cell, not chemistry


def load_index(root: Path = ROOT) -> list:
    """All INDEX.tsv rows, each with an absolute 'path' added."""
    with (root / "cases" / "INDEX.tsv").open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        r["path"] = root / "cases" / r["case"] / r["run"]
    return rows


def bench74_dirs(root: Path = ROOT) -> dict:
    """{(bench_set, ref_key): run dir} for the 74-case benchmark."""
    return {(r["bench_set"], r["ref_key"]): r["path"]
            for r in load_index(root) if r["role"] == "bench74"}


def load_structures(root: Path = ROOT) -> dict:
    """{structure id: STRUCTURES.tsv row}."""
    with (root / "cases" / "STRUCTURES.tsv").open(newline="") as fh:
        return {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def composition(elements) -> str:
    """'Li:16 O:16' - element counts, alphabetical, as STRUCTURES.tsv stores them."""
    return " ".join(f"{e}:{n}" for e, n in sorted(Counter(elements).items()))


def run_structure(run_dir: Path):
    """The run's cell as VASP reads it (POSCAR), or the Materials Studio file of
    a structure-only directory. None when the directory holds neither."""
    from .parsers import parse_eoscar, parse_poscar, Structure
    from .xsd2poscar import parse_xsd
    if (run_dir / "POSCAR").is_file():
        st = parse_poscar(run_dir / "POSCAR")
        if st.elements and st.elements[0].startswith("X") and (run_dir / "EOSCAR").is_file():
            st.elements = list(parse_eoscar(run_dir / "EOSCAR").elements)   # VASP4 POSCAR
        return st
    xsd = sorted(run_dir.glob("*.xsd"))
    if xsd:
        latt, atoms = parse_xsd(xsd[0])
        return Structure(latt, [el for el, _ in atoms], [list(f) for _, f in atoms])
    return None


def min_distance(st) -> float:
    """Shortest interatomic distance under periodic images (Angstrom)."""
    best = math.inf
    shifts = list(itertools.product((-1, 0, 1), repeat=3))
    lat = st.lattice
    for i, j in itertools.combinations(range(len(st.elements)), 2):
        df = [st.frac[i][k] - st.frac[j][k] for k in range(3)]
        df = [x - round(x) for x in df]
        for n in shifts:
            f = [df[k] + n[k] for k in range(3)]
            d2 = sum(sum(f[k] * lat[k][c] for k in range(3)) ** 2 for c in range(3))
            best = min(best, d2)
    return math.sqrt(best)


def check_registry(root: Path = ROOT) -> list:
    """Problems that make INDEX.tsv, STRUCTURES.tsv and cases/ disagree; [] when consistent."""
    problems = []
    index, structs = load_index(root), load_structures(root)
    used = set()
    on_disk = {(c.name, r.name) for c in (root / "cases").iterdir() if c.is_dir()
               for r in c.iterdir() if r.is_dir()}
    listed = {(r["case"], r["run"]) for r in index}
    problems += [f"cases/{c}/{r}: directory not in INDEX.tsv" for c, r in sorted(on_disk - listed)]
    problems += [f"cases/{c}/{r}: INDEX.tsv row without directory" for c, r in sorted(listed - on_disk)]
    for sid, s in structs.items():
        for col, allowed in (("status", STATUS), ("kind", KIND), ("confidence", CONFIDENCE)):
            if s[col] not in allowed:
                problems.append(f"STRUCTURES {sid}: {col} {s[col]!r} not one of {sorted(allowed)}")
        if not s["source"].strip():
            problems.append(f"STRUCTURES {sid}: empty source")
    for r in index:
        where = f"cases/{r['case']}/{r['run']}"
        sid = r.get("structure", "")
        if not sid:
            problems.append(f"{where}: no structure id")
            continue
        if sid not in structs:
            problems.append(f"{where}: structure {sid!r} not in STRUCTURES.tsv")
            continue
        used.add(sid)
        s = structs[sid]
        if s["case"] != r["case"]:
            problems.append(f"{where}: structure {sid} belongs to case {s['case']}")
        if r["role"] in SCORED and s["status"] == "unchecked":
            problems.append(f"{where}: scored run on an unchecked structure {sid}")
        st = run_structure(r["path"])
        if st is None:
            continue
        comp = composition(st.elements)
        if comp != s["composition"]:
            problems.append(f"{where}: composition {comp} differs from {sid} ({s['composition']})")
        if len(st.elements) > 1 and min_distance(st) < MIN_DIST and s["status"] != "wrong":
            problems.append(f"{where}: atoms closer than {MIN_DIST} A but {sid} is not marked wrong")
    problems += [f"STRUCTURES {sid}: not used by any run" for sid in sorted(set(structs) - used)]
    return problems
