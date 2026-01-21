#!/usr/bin/env python3
"""
Utility to summarize Bader charge results from a Materials Studio XSD file
and the accompanying AVF.dat mapping produced by the Henkelman group tools.

The script only relies on the Python standard library so it can run inside
restricted environments where installing extra packages is difficult.
"""
from __future__ import annotations

import argparse
import math
from collections import defaultdict
from pathlib import Path
import re
from statistics import mean, pstdev
from typing import Dict, Iterable, List, Sequence, Tuple


def infer_element(name: str) -> str:
    """Infer the chemical symbol from a label such as 'As12' or 'O_3'."""
    cleaned = []
    for char in name:
        if char.isalpha():
            cleaned.append(char)
        elif cleaned:
            break
    return "".join(cleaned) or name


def parse_xsd(path: Path) -> List[Dict[str, object]]:
    """Return ordered atom data extracted from the Materials Studio XSD."""
    atom_pattern = re.compile(r"<Atom3d\b([^>]*)/?>", re.IGNORECASE)
    attr_pattern = re.compile(r'(\w+)\s*=\s*"([^"]*)"')
    text = path.read_text(encoding="latin1")
    atoms = []
    for match in atom_pattern.finditer(text):
        attr_block = match.group(1)
        attributes = {key: value for key, value in attr_pattern.findall(attr_block)}
        charge_str = attributes.get("Charge", "")
        try:
            charge = float(charge_str)
        except ValueError:
            charge = float("nan")
        atoms.append(
            {
                "index": len(atoms) + 1,
                "name": attributes.get("Name", f"Atom{len(atoms)+1}"),
                "element": attributes.get("Components") or infer_element(attributes.get("Name", "")),
                "charge": charge,
            }
        )
    if not atoms:
        raise ValueError(f"No <Atom3d> entries found in {path}")
    return atoms


def parse_avf(path: Path) -> Dict[int, List[int]]:
    """
    Parse AVF.dat, mapping atom indices (1-based) to lists of basin volume IDs.
    The '...' column present in the original Henkelman output is ignored.
    """
    mapping: Dict[int, List[int]] = defaultdict(list)
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("-") or stripped.startswith("Atom"):
                continue
            tokens = stripped.split()
            try:
                idx = int(tokens[0])
            except (IndexError, ValueError):
                continue
            for token in tokens[1:]:
                if token.isdigit():
                    mapping[idx].append(int(token))
    return mapping


def format_markdown_table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    """Render a Markdown table using the provided headers and rows."""
    table_lines = []
    header_line = "| " + " | ".join(headers) + " |"
    divider = "| " + " | ".join("---" for _ in headers) + " |"
    table_lines.append(header_line)
    table_lines.append(divider)
    for row in rows:
        table_lines.append("| " + " | ".join(row) + " |")
    return "\n".join(table_lines)


def safe_format(value: float, precision: int = 3) -> str:
    """Format floats while handling NaN gracefully."""
    if math.isnan(value):
        return "nan"
    return f"{value:+.{precision}f}"


def build_atom_table(atoms: List[Dict[str, object]], basin_map: Dict[int, List[int]]) -> Tuple[str, List[float]]:
    """Return Markdown table for per-atom information plus the list of valid charges."""
    rows = []
    charges = []
    for atom in atoms:
        idx = atom["index"]
        charge = float(atom["charge"])
        if not math.isnan(charge):
            charges.append(charge)
        volumes = basin_map.get(idx, [])
        rows.append(
            [
                str(idx),
                str(atom["name"]),
                str(atom["element"]),
                safe_format(charge),
                ", ".join(map(str, volumes)) if volumes else "-",
            ]
        )
    table = format_markdown_table(
        ["Idx", "Atom", "Element", "Charge (e)", "Volume IDs"],
        rows,
    )
    return table, charges


def build_element_table(atoms: List[Dict[str, object]]) -> str:
    """Return Markdown table aggregating charge statistics per element."""
    bucket: Dict[str, List[float]] = defaultdict(list)
    for atom in atoms:
        charge = float(atom["charge"])
        if math.isnan(charge):
            continue
        bucket[str(atom["element"])].append(charge)
    rows = []
    for element in sorted(bucket.keys()):
        charges = bucket[element]
        rows.append(
            [
                element,
                str(len(charges)),
                safe_format(sum(charges)),
                safe_format(mean(charges)),
                safe_format(min(charges)),
                safe_format(max(charges)),
                safe_format(pstdev(charges) if len(charges) > 1 else 0.0),
            ]
        )
    if not rows:
        return "*(No per-element statistics available, all charges are NaN.)*"
    return format_markdown_table(
        [
            "Element",
            "Count",
            "Total Charge",
            "Mean Charge",
            "Min",
            "Max",
            "σ (pop)",
        ],
        rows,
    )


def summarize_charges(charges: Iterable[float]) -> str:
    """Return a short textual summary (total, positive/negative balance, std dev)."""
    charges = list(charges)
    if not charges:
        return "未检测到可用的净电荷数据。"
    total = sum(charges)
    positive = sum(c for c in charges if c > 0)
    negative = sum(c for c in charges if c < 0)
    dispersion = pstdev(charges) if len(charges) > 1 else 0.0
    return (
        f"总净电荷 {total:+.4f} e；得电子的原子累积 {negative:+.4f} e，"
        f"失电子的原子累积 {positive:+.4f} e；总体分布标准差 {dispersion:.4f} e。"
    )


def run(args: argparse.Namespace) -> None:
    xsd_path = Path(args.xsd).resolve()
    avf_path = Path(args.avf).resolve() if args.avf else None
    atoms = parse_xsd(xsd_path)
    basin_map = parse_avf(avf_path) if avf_path else {}
    atom_table, charges = build_atom_table(atoms, basin_map)
    element_table = build_element_table(atoms)

    report_lines = [
        f"# Bader 统计报告",
        f"- 数据文件: `{xsd_path.name}`",
        f"- AVF 映射: `{avf_path.name if avf_path else '无'}`",
        f"- 原子总数: {len(atoms)}",
        f"- 拆分盆的原子数: {sum(1 for vols in basin_map.values() if len(vols) > 1)}",
        "",
        "## 全局统计",
        summarize_charges(charges),
        "",
        "## 原子层级统计",
        atom_table,
        "",
        "## 元素层级统计",
        element_table,
        "",
        "## 注释",
        "- 当 AVF.dat 提供多个体积编号时，可视为该原子拥有多个不连通的 Bader 盆。",
    ]
    report = "\n".join(report_lines).strip() + "\n"

    output_path = Path(args.out).resolve()
    output_path.write_text(report, encoding="utf-8")
    if not args.quiet:
        print(report)
        print(f"\n已写入报告: {output_path}")


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Quick statistics + oxidation-state estimates from Bader data.")
    parser.add_argument("--xsd", required=True, help="Materials Studio XSD 文件路径，需包含 <Atom3d> Charge 属性。")
    parser.add_argument("--avf", help="AVF.dat 文件路径（可选，用于展示体积编号）。")
    parser.add_argument(
        "--out",
        default="bader_summary.md",
        help="输出 Markdown 报告路径，默认写入 bader_summary.md。",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="仅写文件，不在终端输出表格。",
    )
    args = parser.parse_args(argv)
    run(args)


if __name__ == "__main__":
    main()
