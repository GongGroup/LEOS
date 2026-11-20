# Bader charge report
- Data file: `WCl5-charge.xsd`
- AVF map: `AVF.dat`
- Atoms: 11
- Atoms with split basins: 10

## Totals
Net charge +0.0010 e; atoms that gained electrons sum to -2.6840 e, atoms that lost electrons to +2.6850 e; standard deviation over all atoms 0.6906 e.

## Per atom
| Idx | Atom | Element | Charge (e) | Volume IDs |
| --- | --- | --- | --- | --- |
| 1 | Cl1 | Cl | -0.197 | 1, 8, 10, 11, 12 |
| 2 | Cl2 | Cl | -0.233 | 14, 17, 21, 24, 40, 41, 45, 46, 47 |
| 3 | Cl3 | Cl | -0.235 | 3, 4, 6, 49, 51, 60, 68, 70, 74 |
| 4 | Cl4 | Cl | +0.255 | 37, 38, 52, 53, 54, 55, 56, 59, 63, 64, 65, 66, 67 |
| 5 | Cl5 | Cl | -0.569 | 28, 50 |
| 6 | Cl6 | Cl | -0.228 | 16, 18, 19, 22, 23, 25, 32, 35, 36, 42, 48 |
| 7 | Cl7 | Cl | +1.975 | 7, 9, 13, 61, 71, 73, 77 |
| 8 | Cl8 | Cl | -0.246 | 2, 5, 15, 39, 43, 44, 62, 69, 72, 75 |
| 9 | Cl9 | Cl | -0.569 | 26, 27, 29, 30, 31 |
| 10 | W10 | W | +0.455 | 20 |
| 11 | W11 | W | -0.407 | 33, 34, 57, 58, 76 |

## Per element
| Element | Count | Total Charge | Mean Charge | Min | Max | σ (pop) |
| --- | --- | --- | --- | --- | --- | --- |
| Cl | 9 | -0.047 | -0.005 | -0.569 | +1.975 | +0.736 |
| W | 2 | +0.048 | +0.024 | -0.407 | +0.455 | +0.431 |

## Notes
- Several volume IDs for one atom in AVF.dat mean that the atom has several disconnected Bader basins.
