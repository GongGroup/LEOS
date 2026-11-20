# Bader charge report
- Data file: `H2S2_hand-charge.xsd`
- AVF map: `AVF.dat`
- Atoms: 4
- Atoms with split basins: 2

## Totals
Net charge +0.0000 e; atoms that gained electrons sum to -1.8470 e, atoms that lost electrons to +1.8470 e; standard deviation over all atoms 0.9258 e.

## Per atom
| Idx | Atom | Element | Charge (e) | Volume IDs |
| --- | --- | --- | --- | --- |
| 1 | H1 | H | -0.869 | 2, 3 |
| 2 | H2 | H | -0.978 | 5, 6 |
| 3 | S3 | S | +0.849 | 1 |
| 4 | S4 | S | +0.998 | 4 |

## Per element
| Element | Count | Total Charge | Mean Charge | Min | Max | σ (pop) |
| --- | --- | --- | --- | --- | --- | --- |
| H | 2 | -1.847 | -0.923 | -0.978 | -0.869 | +0.054 |
| S | 2 | +1.847 | +0.923 | +0.849 | +0.998 | +0.075 |

## Notes
- Several volume IDs for one atom in AVF.dat mean that the atom has several disconnected Bader basins.
