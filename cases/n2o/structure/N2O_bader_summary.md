# Bader charge report
- Data file: `N2O-charge.xsd`
- AVF map: `AVF.dat`
- Atoms: 12
- Atoms with split basins: 0

## Totals
Net charge +0.0020 e; atoms that gained electrons sum to -1.4120 e, atoms that lost electrons to +1.4140 e; standard deviation over all atoms 0.2825 e.

## Per atom
| Idx | Atom | Element | Charge (e) | Volume IDs |
| --- | --- | --- | --- | --- |
| 1 | N1 | N | +0.339 | 10 |
| 2 | N2 | N | +0.015 | 8 |
| 3 | N3 | N | +0.015 | 2 |
| 4 | N4 | N | +0.339 | 4 |
| 5 | N5 | N | +0.338 | 3 |
| 6 | N6 | N | +0.015 | 5 |
| 7 | N7 | N | +0.015 | 11 |
| 8 | N8 | N | +0.338 | 7 |
| 9 | O9 | O | -0.353 | 1 |
| 10 | O10 | O | -0.353 | 9 |
| 11 | O11 | O | -0.353 | 12 |
| 12 | O12 | O | -0.353 | 6 |

## Per element
| Element | Count | Total Charge | Mean Charge | Min | Max | σ (pop) |
| --- | --- | --- | --- | --- | --- | --- |
| N | 8 | +1.414 | +0.177 | +0.015 | +0.339 | +0.162 |
| O | 4 | -1.412 | -0.353 | -0.353 | -0.353 | +0.000 |

## Notes
- Several volume IDs for one atom in AVF.dat mean that the atom has several disconnected Bader basins.
