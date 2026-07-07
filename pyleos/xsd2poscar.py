#!/usr/bin/env python3
"""Convert Materials Studio XSD structures to POSCAR. Stdlib only.

Usage: python3 pyleos/xsd2poscar.py <in.xsd> <out/POSCAR>
Reads the P1 SpaceGroup lattice (AVector/BVector/CVector) and all Atom3d
entries that are not periodic-image copies (no ImageOf attribute).
"""
from __future__ import annotations

import re
import sys
from collections import OrderedDict
from pathlib import Path


def parse_xsd(path: Path):
    text = path.read_text(encoding="latin1")
    sg = re.search(r"<SpaceGroup\b[^>]*", text)
    if not sg:
        raise ValueError(f"{path}: no SpaceGroup element (non-periodic XSD?)")
    attrs = dict(re.findall(r'(\w+)="([^"]*)"', sg.group(0)))
    latt = []
    for key in ("AVector", "BVector", "CVector"):
        latt.append([float(x) for x in attrs[key].split(",")])

    atoms = []
    for m in re.finditer(r"<Atom3d\b([^>]*?)/?>", text):
        a = dict(re.findall(r'(\w+)="([^"]*)"', m.group(1)))
        if "ImageOf" in a:          # periodic display copy, not a real atom
            continue
        if "XYZ" not in a:
            continue
        el = a.get("Components", "").split(",")[0].strip()
        if not el:
            continue
        frac = [float(x) for x in a["XYZ"].split(",")]
        atoms.append((el, frac))
    if not atoms:
        raise ValueError(f"{path}: no Atom3d entries with XYZ")
    return latt, atoms


def write_poscar(latt, atoms, out: Path, title: str) -> None:
    groups: "OrderedDict[str, list]" = OrderedDict()
    for el, frac in atoms:
        groups.setdefault(el, []).append(frac)
    lines = [title, "1.0"]
    for v in latt:
        lines.append(f"  {v[0]:18.10f} {v[1]:18.10f} {v[2]:18.10f}")
    lines.append(" ".join(groups))
    lines.append(" ".join(str(len(v)) for v in groups.values()))
    lines.append("Direct")
    for el, fracs in groups.items():
        for f in fracs:
            lines.append(f"  {f[0]:16.10f} {f[1]:16.10f} {f[2]:16.10f}")
    out.write_text("\n".join(lines) + "\n")


def main():
    src, dst = Path(sys.argv[1]), Path(sys.argv[2])
    latt, atoms = parse_xsd(src)
    dst.parent.mkdir(parents=True, exist_ok=True)
    write_poscar(latt, atoms, dst, src.stem)
    from collections import Counter
    print(f"{src.name}: {len(atoms)} atoms {dict(Counter(el for el, _ in atoms))}")


if __name__ == "__main__":
    main()
