#!/usr/bin/env python3
"""
Helper script to organize downloaded Bader outputs per material and call
bader_stats to generate the Markdown summary automatically.
"""
from __future__ import annotations

import argparse
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import bader_stats


def guess_material_name(xsd_path: Path) -> str:
    """Derive a material folder name from the XSD filename."""
    base = xsd_path.stem
    lowered = base.lower()
    suffixes = ["-charge", "_charge", "-chg", "_chg", "-bader", "_bader"]
    for suffix in suffixes:
        if lowered.endswith(suffix):
            base = base[: -len(suffix)]
            break
    sanitized = base.strip().replace(" ", "_")
    return sanitized or "material"


def ensure_dir(path: Path) -> Path:
    """Create the directory if needed and return it."""
    path.mkdir(parents=True, exist_ok=True)
    return path


def relocate_file(source: Path, target_dir: Path) -> Path:
    """Move source into target_dir unless it is already there."""
    target_dir = ensure_dir(target_dir)
    destination = target_dir / source.name
    if source.resolve() == destination.resolve():
        return destination
    if destination.exists():
        # If a file with the same name exists, overwrite it to keep latest download.
        try:
            destination.unlink()
        except PermissionError:
            # If locked, overwrite via copy later.
            pass
    try:
        shutil.move(str(source), str(destination))
    except PermissionError:
        # Fallback to copy when the file is locked by another process.
        shutil.copy2(str(source), str(destination))
    return destination


def autodetect_avf(xsd_file: Path) -> Optional[Path]:
    """Return AVF path located alongside the XSD if available."""
    candidates = [
        xsd_file.with_name("AVF.dat"),
        xsd_file.with_suffix(".avf"),
        xsd_file.with_suffix(".AVF"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


@dataclass
class RunConfig:
    xsd: Path
    avf: Optional[Path]
    material: str
    material_dir: Path
    summary_dir: Optional[Path]


def build_config(args: argparse.Namespace) -> RunConfig:
    workspace_root = Path(__file__).resolve().parent.parent
    xsd_path = Path(args.xsd).resolve()
    if not xsd_path.exists():
        raise FileNotFoundError(f"XSD file not found: {xsd_path}")
    avf_path = Path(args.avf).resolve() if args.avf else autodetect_avf(xsd_path)
    if args.avf and not avf_path.exists():
        raise FileNotFoundError(f"AVF file not found: {avf_path}")
    material_name = args.material or guess_material_name(xsd_path)
    root_dir = Path(args.root).resolve() if args.root else workspace_root / "EOS"
    material_dir = root_dir / material_name
    summary_dir = ensure_dir(Path(args.summary_dir).resolve()) if args.summary_dir else None
    return RunConfig(
        xsd=xsd_path,
        avf=avf_path,
        material=material_name,
        material_dir=material_dir,
        summary_dir=summary_dir,
    )


def run(args: argparse.Namespace) -> Path:
    cfg = build_config(args)
    relocated_xsd = relocate_file(cfg.xsd, cfg.material_dir)
    relocated_avf = relocate_file(cfg.avf, cfg.material_dir) if cfg.avf else None

    summary_path = (
        cfg.summary_dir / f"{cfg.material}_bader_summary.md"
        if cfg.summary_dir
        else cfg.material_dir / f"{cfg.material}_bader_summary.md"
    )
    bader_stats.run(
        argparse.Namespace(
            xsd=str(relocated_xsd),
            avf=str(relocated_avf) if relocated_avf else None,
            out=str(summary_path),
            quiet=True,
        )
    )
    return summary_path


def main() -> None:
    default_root = Path(__file__).resolve().parent.parent / "EOS"
    parser = argparse.ArgumentParser(
        description="Organize downloaded Bader files into EOS and generate Markdown summaries."
    )
    parser.add_argument("--xsd", required=True, help="The downloaded XSD file.")
    parser.add_argument("--avf", help="The matching AVF.dat; if omitted, it is looked for in the same directory.")
    parser.add_argument("--material", help="Name of the material directory (default: guessed from the XSD file name).")
    parser.add_argument(
        "--root",
        default=str(default_root),
        help="Root of the material directories (default: the EOS directory of this repository).",
    )
    parser.add_argument(
        "--summary-dir",
        help="Directory for the Markdown summary (default: the material directory).",
    )
    args = parser.parse_args()
    summary_path = run(args)
    print(f"Done; summary written to {summary_path}")


if __name__ == "__main__":
    main()
