# Bader charge report
- Data file: `H2S3-charge.xsd`
- AVF map: `AVF.dat`
- Atoms: 5
- Atoms with split basins: 3

## Totals
Net charge +0.0000 e; atoms that gained electrons sum to -0.0490 e, atoms that lost electrons to +0.0490 e; standard deviation over all atoms 0.0251 e.

## Per atom
| Idx | Atom | Element | Charge (e) | Volume IDs |
| --- | --- | --- | --- | --- |
| 1 | H1 | H | +0.003 | 6 |
| 2 | H2 | H | +0.036 | 3 |
| 3 | S3 | S | -0.008 | 5, 8 |
| 4 | S4 | S | +0.010 | 4, 7 |
| 5 | S5 | S | -0.041 | 1, 2 |

## Per element
| Element | Count | Total Charge | Mean Charge | Min | Max | σ (pop) |
| --- | --- | --- | --- | --- | --- | --- |
| H | 2 | +0.039 | +0.019 | +0.003 | +0.036 | +0.016 |
| S | 3 | -0.039 | -0.013 | -0.041 | +0.010 | +0.021 |

## Notes
- Several volume IDs for one atom in AVF.dat mean that the atom has several disconnected Bader basins.
