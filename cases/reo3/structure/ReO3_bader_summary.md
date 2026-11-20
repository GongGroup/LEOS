# Bader charge report
- Data file: `ReO3-charge.xsd`
- AVF map: `AVF.dat`
- Atoms: 4
- Atoms with split basins: 2

## Totals
Net charge +0.0000 e; atoms that gained electrons sum to -2.7100 e, atoms that lost electrons to +2.7100 e; standard deviation over all atoms 1.5646 e.

## Per atom
| Idx | Atom | Element | Charge (e) | Volume IDs |
| --- | --- | --- | --- | --- |
| 1 | O1 | O | -0.904 | 5, 6, 7 |
| 2 | O2 | O | -0.903 | 3, 4 |
| 3 | O3 | O | -0.903 | 2 |
| 4 | Re4 | Re | +2.710 | 1 |

## Per element
| Element | Count | Total Charge | Mean Charge | Min | Max | σ (pop) |
| --- | --- | --- | --- | --- | --- | --- |
| O | 3 | -2.710 | -0.903 | -0.904 | -0.903 | +0.000 |
| Re | 1 | +2.710 | +2.710 | +2.710 | +2.710 | +0.000 |

## Notes
- Several volume IDs for one atom in AVF.dat mean that the atom has several disconnected Bader basins.
