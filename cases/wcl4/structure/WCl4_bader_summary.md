# Bader charge report
- Data file: `WCl4-charge.xsd`
- AVF map: `AVF.dat`
- Atoms: 4
- Atoms with split basins: 1

## Totals
Net charge -0.0010 e; atoms that gained electrons sum to -0.8840 e, atoms that lost electrons to +0.8830 e; standard deviation over all atoms 0.5106 e.

## Per atom
| Idx | Atom | Element | Charge (e) | Volume IDs |
| --- | --- | --- | --- | --- |
| 1 | Cl1 | Cl | -0.274 | 3 |
| 2 | Cl2 | Cl | -0.273 | 4, 5, 6, 7 |
| 3 | Cl3 | Cl | -0.337 | 2 |
| 4 | W4 | W | +0.883 | 1 |

## Per element
| Element | Count | Total Charge | Mean Charge | Min | Max | σ (pop) |
| --- | --- | --- | --- | --- | --- | --- |
| Cl | 3 | -0.884 | -0.295 | -0.337 | -0.273 | +0.030 |
| W | 1 | +0.883 | +0.883 | +0.883 | +0.883 | +0.000 |

## Notes
- Several volume IDs for one atom in AVF.dat mean that the atom has several disconnected Bader basins.
