# Bader charge report
- Data file: `TiO-charge.xsd`
- AVF map: `AVF.dat`
- Atoms: 6
- Atoms with split basins: 0

## Totals
Net charge +0.0010 e; atoms that gained electrons sum to -3.5600 e, atoms that lost electrons to +3.5610 e; standard deviation over all atoms 1.1910 e.

## Per atom
| Idx | Atom | Element | Charge (e) | Volume IDs |
| --- | --- | --- | --- | --- |
| 1 | O1 | O | -1.188 | 5 |
| 2 | O2 | O | -1.188 | 2 |
| 3 | O3 | O | -1.184 | 6 |
| 4 | Ti4 | Ti | +0.988 | 1 |
| 5 | Ti5 | Ti | +1.294 | 3 |
| 6 | Ti6 | Ti | +1.279 | 4 |

## Per element
| Element | Count | Total Charge | Mean Charge | Min | Max | σ (pop) |
| --- | --- | --- | --- | --- | --- | --- |
| O | 3 | -3.560 | -1.187 | -1.188 | -1.184 | +0.002 |
| Ti | 3 | +3.561 | +1.187 | +0.988 | +1.294 | +0.141 |

## Notes
- Several volume IDs for one atom in AVF.dat mean that the atom has several disconnected Bader basins.
