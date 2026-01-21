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
    summary_dir: Path


def build_config(args: argparse.Namespace) -> RunConfig:
    xsd_path = Path(args.xsd).resolve()
    if not xsd_path.exists():
        raise FileNotFoundError(f"找不到 XSD 文件: {xsd_path}")
    avf_path = Path(args.avf).resolve() if args.avf else autodetect_avf(xsd_path)
    if args.avf and not avf_path.exists():
        raise FileNotFoundError(f"找不到 AVF 文件: {avf_path}")
    material_name = args.material or guess_material_name(xsd_path)
    material_dir = Path(args.root).resolve() / material_name
    summary_dir = ensure_dir(Path(args.summary_dir).resolve())
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

    summary_path = cfg.summary_dir / f"{cfg.material}_bader_summary.md"
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
    parser = argparse.ArgumentParser(
        description="Organize downloaded Bader files and generate Markdown summaries."
    )
    parser.add_argument("--xsd", required=True, help="刚下载的 XSD 文件路径。")
    parser.add_argument("--avf", help="对应的 AVF.dat 路径；若省略将尝试在同目录中自动查找。")
    parser.add_argument("--material", help="手动指定物质文件夹名称，默认依据 XSD 文件名猜测。")
    parser.add_argument("--root", default=".", help="物质子目录的根路径，默认为当前目录。")
    parser.add_argument("--summary-dir", default="summary", help="汇总 Markdown 输出目录。")
    args = parser.parse_args()
    summary_path = run(args)
    print(f"整理完成，汇总写入: {summary_path}")


if __name__ == "__main__":
    main()
