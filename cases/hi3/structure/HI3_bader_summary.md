# Bader charge report
- Data file: `HI3-charge.xsd`
- AVF map: `AVF.dat`
- Atoms: 8
- Atoms with split basins: 0

## Totals
Net charge +0.0000 e; atoms that gained electrons sum to -0.1880 e, atoms that lost electrons to +0.1880 e; standard deviation over all atoms 0.0523 e.

## Per atom
| Idx | Atom | Element | Charge (e) | Volume IDs |
| --- | --- | --- | --- | --- |
| 1 | H1 | H | +0.037 | 8 |
| 2 | H2 | H | +0.037 | 7 |
| 3 | I3 | I | -0.078 | 2 |
| 4 | I4 | I | +0.057 | 3 |
| 5 | I5 | I | -0.016 | 4 |
| 6 | I6 | I | -0.016 | 5 |
| 7 | I7 | I | +0.057 | 6 |
| 8 | I8 | I | -0.078 | 1 |

## Per element
| Element | Count | Total Charge | Mean Charge | Min | Max | σ (pop) |
| --- | --- | --- | --- | --- | --- | --- |
| H | 2 | +0.074 | +0.037 | +0.037 | +0.037 | +0.000 |
| I | 6 | -0.074 | -0.012 | -0.078 | +0.057 | +0.055 |

## Notes
- Several volume IDs for one atom in AVF.dat mean that the atom has several disconnected Bader basins.
