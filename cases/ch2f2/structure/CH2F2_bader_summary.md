# Bader charge report
- Data file: `CH2F2-charge.xsd`
- AVF map: `AVF.dat`
- Atoms: 5
- Atoms with split basins: 1

## Totals
Net charge +0.0000 e; atoms that gained electrons sum to -1.2260 e, atoms that lost electrons to +1.2260 e; standard deviation over all atoms 0.6015 e.

## Per atom
| Idx | Atom | Element | Charge (e) | Volume IDs |
| --- | --- | --- | --- | --- |
| 1 | H1 | H | +0.099 | 5 |
| 2 | H2 | H | +0.110 | 2 |
| 3 | C3 | C | +1.017 | 3, 4 |
| 4 | F4 | F | -0.636 | 1 |
| 5 | F5 | F | -0.590 | 6 |

## Per element
| Element | Count | Total Charge | Mean Charge | Min | Max | σ (pop) |
| --- | --- | --- | --- | --- | --- | --- |
| C | 1 | +1.017 | +1.017 | +1.017 | +1.017 | +0.000 |
| F | 2 | -1.226 | -0.613 | -0.636 | -0.590 | +0.023 |
| H | 2 | +0.209 | +0.105 | +0.099 | +0.110 | +0.005 |

## Notes
- Several volume IDs for one atom in AVF.dat mean that the atom has several disconnected Bader basins.
