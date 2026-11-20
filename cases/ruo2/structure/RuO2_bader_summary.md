# Bader charge report
- Data file: `RuO2-charge.xsd`
- AVF map: `AVF.dat`
- Atoms: 6
- Atoms with split basins: 0

## Totals
Net charge -0.0010 e; atoms that gained electrons sum to -3.2390 e, atoms that lost electrons to +3.2380 e; standard deviation over all atoms 1.1449 e.

## Per atom
| Idx | Atom | Element | Charge (e) | Volume IDs |
| --- | --- | --- | --- | --- |
| 1 | O1 | O | -0.810 | 2 |
| 2 | O2 | O | -0.809 | 5 |
| 3 | O3 | O | -0.810 | 3 |
| 4 | O4 | O | -0.810 | 4 |
| 5 | Ru5 | Ru | +1.619 | 6 |
| 6 | Ru6 | Ru | +1.619 | 1 |

## Per element
| Element | Count | Total Charge | Mean Charge | Min | Max | σ (pop) |
| --- | --- | --- | --- | --- | --- | --- |
| O | 4 | -3.239 | -0.810 | -0.810 | -0.809 | +0.000 |
| Ru | 2 | +3.238 | +1.619 | +1.619 | +1.619 | +0.000 |

## Notes
- Several volume IDs for one atom in AVF.dat mean that the atom has several disconnected Bader basins.
