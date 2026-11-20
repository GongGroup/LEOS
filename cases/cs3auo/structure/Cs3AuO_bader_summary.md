# Bader charge report
- Data file: `Cs3AuO-charge.xsd`
- AVF map: `AVF.dat`
- Atoms: 10
- Atoms with split basins: 0

## Totals
Net charge -0.0000 e; atoms that gained electrons sum to -4.3300 e, atoms that lost electrons to +4.3300 e; standard deviation over all atoms 0.9011 e.

## Per atom
| Idx | Atom | Element | Charge (e) | Volume IDs |
| --- | --- | --- | --- | --- |
| 1 | O1 | O | -1.360 | 1 |
| 2 | O2 | O | -1.360 | 3 |
| 3 | Cs3 | Cs | +0.720 | 5 |
| 4 | Cs4 | Cs | +0.724 | 2 |
| 5 | Cs5 | Cs | +0.721 | 9 |
| 6 | Cs6 | Cs | +0.721 | 6 |
| 7 | Cs7 | Cs | +0.724 | 4 |
| 8 | Cs8 | Cs | +0.720 | 10 |
| 9 | Au9 | Au | -0.805 | 7 |
| 10 | Au10 | Au | -0.805 | 8 |

## Per element
| Element | Count | Total Charge | Mean Charge | Min | Max | σ (pop) |
| --- | --- | --- | --- | --- | --- | --- |
| Au | 2 | -1.610 | -0.805 | -0.805 | -0.805 | +0.000 |
| Cs | 6 | +4.330 | +0.722 | +0.720 | +0.724 | +0.002 |
| O | 2 | -2.720 | -1.360 | -1.360 | -1.360 | +0.000 |

## Notes
- Several volume IDs for one atom in AVF.dat mean that the atom has several disconnected Bader basins.
