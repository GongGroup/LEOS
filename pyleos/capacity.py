"""Element capacity bounds (wide/physical) and electronegativity scales.

Wide capacity encodes only what is physically impossible, not what is
chemically unusual: n_e must satisfy  ZVAL - OS_MAX <= n_e <= ZVAL - OS_MIN.
Bounds are deliberately loose (rare-but-real states must stay feasible);
they exist to prevent catastrophic over-assignment (e.g. S receiving 14 e),
not to steer the answer. Tight/chemical-prior mode is a separate diagnostic.
"""
from __future__ import annotations

from typing import Dict, Tuple

# Physical OS bounds (OS_MIN, OS_MAX) per element. Main group: anion limited
# by valence-shell closure, cation by group valence. TM/f: generous margins
# (Pt(-2), Au(-3) exist; Fe(+6), Re(+7) exist).
OS_BOUNDS: Dict[str, Tuple[int, int]] = {
    "H": (-1, 1), "Li": (-1, 1), "Na": (-1, 1), "K": (-1, 1),
    "Rb": (-1, 1), "Cs": (-1, 1),
    "Be": (0, 2), "Mg": (0, 2), "Ca": (-1, 2), "Sr": (-1, 2), "Ba": (-1, 2),
    "B": (-5, 3), "Al": (-2, 3), "Ga": (-5, 3), "In": (-5, 3),
    "C": (-4, 4), "Si": (-4, 4), "Ge": (-4, 4), "Sn": (-4, 4), "Pb": (-4, 4),
    "N": (-3, 5), "P": (-3, 5), "As": (-3, 5), "Sb": (-3, 5), "Bi": (-3, 5),
    "O": (-2, 2), "S": (-2, 6), "Se": (-2, 6), "Te": (-2, 6),
    "F": (-1, 0), "Cl": (-1, 7), "Br": (-1, 7), "I": (-1, 7),
    "Sc": (-1, 3), "Ti": (-2, 4), "V": (-3, 5), "Cr": (-4, 6),
    "Mn": (-3, 7), "Fe": (-4, 6), "Co": (-3, 5), "Ni": (-2, 4),
    "Cu": (-2, 4), "Zn": (-2, 2),
    "Y": (-1, 3), "Zr": (-2, 4), "Nb": (-3, 5), "Mo": (-4, 6),
    "Ru": (-4, 8), "Rh": (-3, 6), "Pd": (0, 4), "Ag": (-2, 3), "Cd": (-2, 2),
    "La": (0, 3), "Ce": (0, 4), "Hf": (-2, 4), "Ta": (-3, 5),
    "W": (-4, 6), "Re": (-3, 7), "Os": (-4, 8), "Ir": (-3, 8),
    "Pt": (-3, 6), "Au": (-3, 5), "Hg": (-2, 4),
}
DEFAULT_BOUNDS = (-4, 8)

# Allen electronegativity (default; matches IUPAC 2014/2016 recommendation).
ENEG_ALLEN: Dict[str, float] = {
    "H": 2.300, "Li": 0.912, "Be": 1.576, "B": 2.051, "C": 2.544,
    "N": 3.066, "O": 3.610, "F": 4.193, "Na": 0.869, "Mg": 1.293,
    "Al": 1.613, "Si": 1.916, "P": 2.253, "S": 2.589, "Cl": 2.869,
    "K": 0.734, "Ca": 1.034, "Sc": 1.19, "Ti": 1.38, "V": 1.53,
    "Cr": 1.65, "Mn": 1.75, "Fe": 1.80, "Co": 1.84, "Ni": 1.88,
    "Cu": 1.85, "Zn": 1.588, "Ga": 1.756, "Ge": 1.994, "As": 2.211,
    "Se": 2.424, "Br": 2.685, "Rb": 0.706, "Sr": 0.963, "Y": 1.12,
    "Zr": 1.32, "Nb": 1.41, "Mo": 1.47, "Ru": 1.54, "Rh": 1.56,
    "Pd": 1.58, "Ag": 1.87, "Cd": 1.52, "In": 1.656, "Sn": 1.824,
    "Sb": 1.984, "Te": 2.158, "I": 2.359, "Cs": 0.659, "Ba": 0.881,
    "La": 1.09, "Ce": 1.09, "Hf": 1.16, "Ta": 1.34, "W": 1.47,
    "Re": 1.60, "Os": 1.65, "Ir": 1.68, "Pt": 1.72, "Au": 1.92,
    "Hg": 1.76, "Pb": 1.854, "Bi": 2.01,
}

# Pauling scale kept for comparison with the current Fortran greedy.
ENEG_PAULING: Dict[str, float] = {
    "H": 2.20, "Li": 0.98, "Be": 1.57, "B": 2.04, "C": 2.55, "N": 3.04,
    "O": 3.44, "F": 3.98, "Na": 0.93, "Mg": 1.31, "Al": 1.61, "Si": 1.90,
    "P": 2.19, "S": 2.58, "Cl": 3.16, "K": 0.82, "Ca": 1.00, "Sc": 1.36,
    "Ti": 1.54, "V": 1.63, "Cr": 1.66, "Mn": 1.55, "Fe": 1.83, "Co": 1.88,
    "Ni": 1.91, "Cu": 1.90, "Zn": 1.65, "Ga": 1.81, "Ge": 2.01, "As": 2.18,
    "Se": 2.55, "Br": 2.96, "Rb": 0.82, "Sr": 0.95, "Y": 1.22, "Zr": 1.33,
    "Nb": 1.6, "Mo": 2.16, "Ru": 2.2, "Rh": 2.28, "Pd": 2.20, "Ag": 1.93,
    "Cd": 1.69, "In": 1.78, "Sn": 1.96, "Sb": 2.05, "Te": 2.1, "I": 2.66,
    "Cs": 0.79, "Ba": 0.89, "La": 1.10, "Ce": 1.12, "Hf": 1.3, "Ta": 1.5,
    "W": 2.36, "Re": 1.9, "Os": 2.2, "Ir": 2.20, "Pt": 2.28, "Au": 2.54,
    "Hg": 2.00, "Pb": 2.33, "Bi": 2.02,
}


def capacity_bounds(element: str, zval: float) -> Tuple[int, int]:
    """Return (n_min, n_max) electron capacity from the wide physical bounds."""
    os_min, os_max = OS_BOUNDS.get(element, DEFAULT_BOUNDS)
    z = round(zval)
    n_min = max(0, z - os_max)
    n_max = z - os_min
    return n_min, n_max


def eneg(element: str, scale: str = "allen") -> float:
    table = ENEG_ALLEN if scale == "allen" else ENEG_PAULING
    return table.get(element, 0.0)
