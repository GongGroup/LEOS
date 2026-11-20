# Bader charge report
- Data file: `H2S5-charge.xsd`
- AVF map: `AVF.dat`
- Atoms: 7
- Atoms with split basins: 5

## Totals
Net charge +0.0000 e; atoms that gained electrons sum to -0.2150 e, atoms that lost electrons to +0.2150 e; standard deviation over all atoms 0.0736 e.

## Per atom
| Idx | Atom | Element | Charge (e) | Volume IDs |
| --- | --- | --- | --- | --- |
| 1 | H1 | H | +0.045 | 10 |
| 2 | H2 | H | +0.100 | 7 |
| 3 | S3 | S | -0.089 | 8, 15, 16 |
| 4 | S4 | S | -0.003 | 5, 6, 9 |
| 5 | S5 | S | +0.049 | 12, 13, 14, 17 |
| 6 | S6 | S | -0.123 | 1, 2, 11 |
| 7 | S7 | S | +0.021 | 3, 4 |

## Per element
| Element | Count | Total Charge | Mean Charge | Min | Max | σ (pop) |
| --- | --- | --- | --- | --- | --- | --- |
| H | 2 | +0.145 | +0.073 | +0.045 | +0.100 | +0.028 |
| S | 5 | -0.145 | -0.029 | -0.123 | +0.049 | +0.066 |

## Notes
- Several volume IDs for one atom in AVF.dat mean that the atom has several disconnected Bader basins.
