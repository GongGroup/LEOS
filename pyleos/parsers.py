"""Parsers for VASP/Wannier90/LEOS output files. Stdlib only."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def inv3(a: List[List[float]]) -> List[List[float]]:
    """Inverse of a 3x3 matrix (row-vector lattice convention)."""
    d = (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
         - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
         + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))
    cof = [
        [(a[1][1] * a[2][2] - a[1][2] * a[2][1]),
         -(a[0][1] * a[2][2] - a[0][2] * a[2][1]),
         (a[0][1] * a[1][2] - a[0][2] * a[1][1])],
        [-(a[1][0] * a[2][2] - a[1][2] * a[2][0]),
         (a[0][0] * a[2][2] - a[0][2] * a[2][0]),
         -(a[0][0] * a[1][2] - a[0][2] * a[1][0])],
        [(a[1][0] * a[2][1] - a[1][1] * a[2][0]),
         -(a[0][0] * a[2][1] - a[0][1] * a[2][0]),
         (a[0][0] * a[1][1] - a[0][1] * a[1][0])],
    ]
    return [[cof[i][j] / d for j in range(3)] for i in range(3)]


@dataclass
class Structure:
    lattice: List[List[float]]      # rows are lattice vectors (Angstrom)
    elements: List[str]             # per atom, POSCAR order
    frac: List[List[float]]         # fractional coordinates


@dataclass
class WannierChannel:
    centres: List[List[float]]      # cartesian, Angstrom
    spreads: List[float]            # Angstrom^2


@dataclass
class WfDetail:
    spin: int
    index: int
    occ: float                                # electrons carried (RSPIN-scaled)
    centre: Tuple[float, float, float]        # cartesian, Angstrom
    spread: float
    at1: int                                  # greedy primary atom (1-based)
    e1: float
    at2: int
    e2: float


@dataclass
class EoscarData:
    zval: Dict[str, float]                    # element -> ZVAL
    greedy_ne: List[float]                    # per atom, n_e assigned by greedy
    greedy_os: List[int]                      # per atom, integer OS
    elements: List[str]                       # per atom
    nelect_wann: Optional[float] = None
    dne: Optional[float] = None
    n_lost_states: Optional[int] = None       # EOSCAR v2 loss diagnostics
    n_lost_e: Optional[float] = None
    wf_detail: List[WfDetail] = field(default_factory=list)   # v2 #WF block
    wp_detail: Dict[Tuple[int, int], List[Tuple[int, float]]] = field(
        default_factory=dict)   # v3 #WP block: (spin, wf) -> [(atom, weight)]


@dataclass
class CaseData:
    path: Path
    structure: Structure
    channels: List[WannierChannel]  # 1 (ISPIN=1) or 2 (ISPIN=2) entries
    eoscar: Optional[EoscarData]
    ispin: int = 1

    @property
    def name(self) -> str:
        return self.path.name


def parse_poscar(path: Path) -> Structure:
    lines = path.read_text().splitlines()
    scale = float(lines[1].split()[0])
    lattice = [[float(x) * scale for x in lines[i].split()[:3]] for i in (2, 3, 4)]
    maybe_counts = lines[5].split()
    if all(tok.isdigit() for tok in maybe_counts):
        # VASP4-style POSCAR without the species line; caller must map
        # elements from EOSCAR/POTCAR order (load_case handles this).
        species = [f"X{i}" for i in range(len(maybe_counts))]
        counts = [int(x) for x in maybe_counts]
        row = 6
    else:
        species = maybe_counts
        counts = [int(x) for x in lines[6].split()]
        row = 7
    if lines[row].strip().lower().startswith("s"):  # Selective dynamics
        row += 1
    mode = lines[row].strip().lower()
    cartesian = mode.startswith(("c", "k"))
    if not cartesian and not mode.startswith("d"):
        raise ValueError(f"{path}: unknown coordinate mode {mode!r}")
    row += 1
    elements, coords = [], []
    for sp, n in zip(species, counts):
        for _ in range(n):
            coords.append([float(x) for x in lines[row].split()[:3]])
            elements.append(sp)
            row += 1
    if cartesian:
        binv = inv3(lattice)
        coords = [[sum(c[i] * scale * binv[i][j] for i in range(3)) for j in range(3)]
                  for c in coords]
    return Structure(lattice, elements, coords)


def parse_wout(path: Path) -> WannierChannel:
    """Extract the Final State WF centres and spreads from wannier90 .wout."""
    centres, spreads = [], []
    in_final = False
    for line in path.read_text(errors="replace").splitlines():
        if "Final State" in line:
            in_final = True
            centres, spreads = [], []   # keep the last Final State block
            continue
        if in_final:
            if "WF centre and spread" in line:
                body = line.split("(", 1)[1]
                xyz, tail = body.split(")", 1)
                centres.append([float(x) for x in xyz.split(",")])
                spreads.append(float(tail.split()[0]))
            elif centres:
                in_final = False
    if not centres:
        raise ValueError(f"{path}: no Final State block found")
    return WannierChannel(centres, spreads)


def parse_eoscar(path: Path) -> EoscarData:
    zval: Dict[str, float] = {}
    elements: List[str] = []
    ne: List[float] = []
    osi: List[int] = []
    nelect_wann = dne = None
    n_lost_states = n_lost_e = None
    wf_detail: List[WfDetail] = []
    wp_detail: Dict[Tuple[int, int], List[Tuple[int, float]]] = {}
    for line in path.read_text().splitlines():
        s = line.strip()
        if s.startswith("#"):
            if "NELECT(Wannier)" in s:
                nelect_wann = float(s.split("=")[1])
            elif "|dNe|" in s:
                dne = float(s.split("=")[1])
            elif s.startswith("# N_LOST_STATES"):
                n_lost_states = int(s.split("=")[1])
            elif s.startswith("# N_LOST_E"):
                n_lost_e = float(s.split("=")[1])
            elif s.startswith("#WP") and len(s.split()) >= 11 and s.split()[1].isdigit():
                q = s.split()
                pairs = []
                for k in range(3, 11, 2):
                    at, w = int(q[k]), float(q[k + 1])
                    if at > 0 and w > 1e-4:
                        pairs.append((at, w))
                wp_detail[(int(q[1]), int(q[2]))] = pairs
            elif s.startswith("#WF") and len(s.split()) >= 11 and s.split()[1].isdigit():
                q = s.split()
                wf_detail.append(WfDetail(
                    int(q[1]), int(q[2]), float(q[3]),
                    (float(q[4]), float(q[5]), float(q[6])),
                    float(q[7]), int(q[8]), float(q[9]),
                    int(q[10]), float(q[11]) if len(q) > 11 else 0.0))
            continue
        parts = s.split()
        if len(parts) >= 6 and parts[0].isdigit():
            el = parts[1]
            elements.append(el)
            zval[el] = float(parts[2])
            ne.append(float(parts[3]))
            osi.append(int(parts[5]))
    if not elements:
        raise ValueError(f"{path}: no atom rows parsed")
    return EoscarData(zval, ne, osi, elements, nelect_wann, dne,
                      n_lost_states, n_lost_e, wf_detail, wp_detail)


def load_case(case_dir: Path) -> CaseData:
    """Load one calculation directory (POSCAR + wout channel(s) + EOSCAR)."""
    case_dir = Path(case_dir)
    structure = parse_poscar(case_dir / "POSCAR")
    spin_files = [case_dir / "wannier90.1.wout", case_dir / "wannier90.2.wout"]
    if all(p.is_file() for p in spin_files):
        channels = [parse_wout(p) for p in spin_files]
        ispin = 2
    else:
        channels = [parse_wout(case_dir / "wannier90.wout")]
        ispin = 1
    eoscar_path = case_dir / "EOSCAR"
    eoscar = parse_eoscar(eoscar_path) if eoscar_path.is_file() else None
    if structure.elements and structure.elements[0].startswith("X") and eoscar:
        if len(eoscar.elements) == len(structure.elements):
            structure.elements = list(eoscar.elements)   # VASP4 POSCAR fallback
    return CaseData(case_dir, structure, channels, eoscar, ispin)
