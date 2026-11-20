# Bader charge report
- Data file: `H2S2-charge.xsd`
- AVF map: `AVF.dat`
- Atoms: 4
- Atoms with split basins: 0

## Totals
Net charge -0.0000 e; atoms that gained electrons sum to -0.0410 e, atoms that lost electrons to +0.0410 e; standard deviation over all atoms 0.0247 e.

## Per atom
| Idx | Atom | Element | Charge (e) | Volume IDs |
| --- | --- | --- | --- | --- |
| 1 | H1 | H | +0.007 | 3 |
| 2 | H2 | H | +0.025 | 4 |
| 3 | S3 | S | -0.041 | 1 |
| 4 | S4 | S | +0.009 | 2 |

## Per element
| Element | Count | Total Charge | Mean Charge | Min | Max | σ (pop) |
| --- | --- | --- | --- | --- | --- | --- |
| H | 2 | +0.032 | +0.016 | +0.007 | +0.025 | +0.009 |
| S | 2 | -0.032 | -0.016 | -0.041 | +0.009 | +0.025 |

## Notes
- Several volume IDs for one atom in AVF.dat mean that the atom has several disconnected Bader basins.
