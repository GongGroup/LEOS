# Bader charge report
- Data file: `WCl3-charge.xsd`
- AVF map: `AVF.dat`
- Atoms: 4
- Atoms with split basins: 4

## Totals
Net charge -0.0010 e; atoms that gained electrons sum to -1.3620 e, atoms that lost electrons to +1.3610 e; standard deviation over all atoms 0.7937 e.

## Per atom
| Idx | Atom | Element | Charge (e) | Volume IDs |
| --- | --- | --- | --- | --- |
| 1 | Cl1 | Cl | -0.326 | 1, 2, 3, 4, 5, 6, 7, 8, 9, 22, 23, 24, 25, 26, 27, 29, 30, 31, 32 |
| 2 | Cl2 | Cl | -0.407 | 10, 11, 13, 14, 17, 18, 19, 21, 34 |
| 3 | Cl3 | Cl | -0.629 | 15, 16, 20, 28, 33, 35, 36, 37, 39, 40, 41, 42, 43, 44 |
| 4 | W4 | W | +1.361 | 12, 38 |

## Per element
| Element | Count | Total Charge | Mean Charge | Min | Max | σ (pop) |
| --- | --- | --- | --- | --- | --- | --- |
| Cl | 3 | -1.362 | -0.454 | -0.629 | -0.326 | +0.128 |
| W | 1 | +1.361 | +1.361 | +1.361 | +1.361 | +0.000 |

## Notes
- Several volume IDs for one atom in AVF.dat mean that the atom has several disconnected Bader basins.
