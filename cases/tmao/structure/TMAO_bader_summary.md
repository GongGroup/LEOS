# Bader charge report
- Data file: `TMAO-charge.xsd`
- AVF map: `AVF.dat`
- Atoms: 14
- Atoms with split basins: 1

## Totals
Net charge +0.0000 e; atoms that gained electrons sum to -1.0060 e, atoms that lost electrons to +1.0060 e; standard deviation over all atoms 0.2148 e.

## Per atom
| Idx | Atom | Element | Charge (e) | Volume IDs |
| --- | --- | --- | --- | --- |
| 1 | H1 | H | +0.118 | 8 |
| 2 | H2 | H | +0.098 | 14 |
| 3 | H3 | H | -0.016 | 4 |
| 4 | H4 | H | +0.084 | 1 |
| 5 | H5 | H | +0.097 | 12 |
| 6 | H6 | H | +0.095 | 15 |
| 7 | H7 | H | +0.004 | 7 |
| 8 | H8 | H | +0.141 | 3 |
| 9 | H9 | H | +0.083 | 11 |
| 10 | C10 | C | +0.120 | 2 |
| 11 | C11 | C | +0.043 | 13 |
| 12 | C12 | C | +0.123 | 6 |
| 13 | N13 | N | -0.339 | 9 |
| 14 | O14 | O | -0.651 | 5, 10 |

## Per element
| Element | Count | Total Charge | Mean Charge | Min | Max | σ (pop) |
| --- | --- | --- | --- | --- | --- | --- |
| C | 3 | +0.286 | +0.095 | +0.043 | +0.123 | +0.037 |
| H | 9 | +0.704 | +0.078 | -0.016 | +0.141 | +0.048 |
| N | 1 | -0.339 | -0.339 | -0.339 | -0.339 | +0.000 |
| O | 1 | -0.651 | -0.651 | -0.651 | -0.651 | +0.000 |

## Notes
- Several volume IDs for one atom in AVF.dat mean that the atom has several disconnected Bader basins.
