# Bader charge report
- Data file: `(Cu5I7)2--charge.xsd`
- AVF map: `AVF.dat`
- Atoms: 12
- Atoms with split basins: 2

## Totals
Net charge -2.0000 e; atoms that gained electrons sum to -3.3390 e, atoms that lost electrons to +1.3390 e; standard deviation over all atoms 0.3760 e.

## Per atom
| Idx | Atom | Element | Charge (e) | Volume IDs |
| --- | --- | --- | --- | --- |
| 1 | Cu1 | Cu | +0.246 | 12 |
| 2 | Cu2 | Cu | +0.249 | 8 |
| 3 | Cu3 | Cu | +0.274 | 4 |
| 4 | Cu4 | Cu | +0.274 | 9 |
| 5 | Cu5 | Cu | +0.296 | 13 |
| 6 | I6 | I | -0.563 | 6 |
| 7 | I7 | I | -0.535 | 2 |
| 8 | I8 | I | -0.310 | 10 |
| 9 | I9 | I | -0.532 | 5 |
| 10 | I10 | I | -0.314 | 11 |
| 11 | I11 | I | -0.540 | 1, 3, 15, 16 |
| 12 | I12 | I | -0.545 | 7, 14 |

## Per element
| Element | Count | Total Charge | Mean Charge | Min | Max | σ (pop) |
| --- | --- | --- | --- | --- | --- | --- |
| Cu | 5 | +1.339 | +0.268 | +0.246 | +0.296 | +0.018 |
| I | 7 | -3.339 | -0.477 | -0.563 | -0.310 | +0.105 |

## Notes
- Several volume IDs for one atom in AVF.dat mean that the atom has several disconnected Bader basins.
