# Bader charge report
- Data file: `CrO5-charge.xsd`
- AVF map: `AVF.dat`
- Atoms: 6
- Atoms with split basins: 3

## Totals
Net charge +0.0000 e; atoms that gained electrons sum to -1.5280 e, atoms that lost electrons to +1.5280 e; standard deviation over all atoms 0.6838 e.

## Per atom
| Idx | Atom | Element | Charge (e) | Volume IDs |
| --- | --- | --- | --- | --- |
| 1 | O1 | O | -0.337 | 6, 8, 9 |
| 2 | O2 | O | -0.289 | 2, 3 |
| 3 | O3 | O | -0.338 | 7, 10 |
| 4 | O4 | O | -0.283 | 1 |
| 5 | O5 | O | -0.281 | 4 |
| 6 | Cr6 | Cr | +1.528 | 5 |

## Per element
| Element | Count | Total Charge | Mean Charge | Min | Max | σ (pop) |
| --- | --- | --- | --- | --- | --- | --- |
| Cr | 1 | +1.528 | +1.528 | +1.528 | +1.528 | +0.000 |
| O | 5 | -1.528 | -0.306 | -0.338 | -0.281 | +0.026 |

## Notes
- Several volume IDs for one atom in AVF.dat mean that the atom has several disconnected Bader basins.
