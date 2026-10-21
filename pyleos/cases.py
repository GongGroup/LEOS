"""Case registry: cases/INDEX.tsv lists every directory under cases/.

One row per run directory cases/<case>/<run>/. Columns:
  case       chemical system (directory name)
  run        run label: program version + what differs (v1, v3, v31-wp, ...)
  role       what the run is used for (bench74, comas, test-regression, ...)
  bench_set  74-case benchmark set (v1_mol / v1_oxide / wave1 / wave2) or empty
  ref_key    key of the case in pyleos/sets/reference_os.json, or empty
  origin     path the run had before the cases/ layout (provenance)
"""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


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
