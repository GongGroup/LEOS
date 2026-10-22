# Fig. 4 – caption draft and provenance

Script: `fig4_delithiation.py` (run `python3 fig4_delithiation.py` in this
directory; `--recompute` ignores `_cache/fig4_*.json`). Output:
`fig4_delithiation.pdf` / `.png`, page 17.1 cm x 12.2 cm. Layout: (a) above
(b) in the left column, (c) over the full height of the right column.

## Caption

**Fig. 4** Delithiation: which species is oxidised when Li is removed.
Circles: min-cost flow (this work); squares: per-WF greedy rule; dashed line
and open diamonds: charge counting. (a) Li2O2 (32-atom Li16O16 cell, PBE):
element-mean oxidation state of O (top) and of Li (bottom) after removal of
0, 1 and 2 Li atoms. Charge counting takes Li as +1 in a neutral cell, i.e.
-n(Li)/16 for O. (b) Location of the hole in Li15O16. Each pair of symbols
joined by a line is one O2 unit, placed at its O-O distance: the circle is
the sum of the flow oxidation states of its two atoms, the square is the
charge obtained from the electrons that the greedy rule counts on them. The
distance axis is broken between 1.60 and 2.11 Å. The six calculations start
from a cell in which one peroxide unit was shortened by hand: PBE+U with U =
0, 4, 6 and 8 eV on O 2p, each on its own relaxed geometry, and PBE and HSE06
single points on the U = 8 eV geometry (PBE // U = 8, HSE06 // U = 8). Units
beyond 2.1 Å are peroxide units pulled apart by U; two O atoms without a
partner in the U = 6 eV geometry are not shown. (c) Oxidation states along
50 ps molecular-dynamics trajectories (MACE-MPA-0, about 300 K) of Li-rich
cathodes, analysed every 2.5 ps (20 frames per trajectory, PBE+U single
points with U on the transition metals). Columns: fraction x of Li removed.
Rows: mean oxidation state of the named element in Li2(1-x)MnO3, in
Li2(1-x)RuO3 and in disordered Li1.2Ni0.13Mn0.54Co0.13O2 (NMC, special
quasi-random structure). The dashed line is the charge-counting value
4 + 2x; the three NMC rows share one vertical scale, which differs from that
of the two upper rows. Once Li is +1 and O is -2, the mean oxidation state of
O in (a) and of the transition metal in the two upper rows of (c) follows
from electron counting; which O2 unit (b) and which metal (c, NMC rows)
carries the hole does not.

## What is fixed by stoichiometry and what is not

- Fixed: the flow assigns +1 to every Li in (a), (b) and (c) and -2 to every
  O in (c). With that, the mean O state in Li2O2 (-1, -0.9375, -0.875) and
  the mean transition-metal state in Li2MnO3 / Li2RuO3 (4 + 2x) are forced.
  The test passed here is that Li and O are counted correctly; the greedy
  rule fails it (Li at -1 in Li2O2, O at -1.67 to -1.99 in the cathodes).
- Not fixed, and not delivered reliably: the site. In (b) the flow gives -2
  to the short-bond unit in all five calculations that have one and places
  the -1 on a unit of peroxide length; R_site is 0.05-0.71 and the smallest
  site margin 0.0005-0.0071 cost units. In the NMC rows of (c) the partition
  among Ni, Mn and Co changes from frame to frame.

## Notes for the text

- The greedy rule is the per-atom table that the VASP patch writes to each
  EOSCAR; it is not a pure nearest-atom rule. Element means use its OS(int)
  column. The unit charges in (b) use its n_e(assigned) column (12 minus the
  electrons on the two atoms), as the table in RESULTS.md does. Sums of the
  per-atom OS(int) differ by up to 1 per unit because a shared pair counted
  half and half is rounded on each atom (U = 8 eV: -1, -2, -3, -2, -2, -2,
  -2, -1); they are printed in the provenance table.
- (b), greedy rule: the short-bond unit carries -1.21 (U = 4 eV, 1.405 Å),
  -1.13 (U = 6 eV, 1.377 Å), -1.11 (U = 8 eV, 1.401 Å), -1.22 (PBE // U = 8)
  and -1.11 (HSE06 // U = 8). With U = 0 no bond is shorter than 1.48 Å and
  no unit is above -1.73.
- (b), flow: the unit at -1 is O16-O20 (1.530 Å, the longest bond) for
  U = 0, O25-O29 (1.550 Å) for U = 4 eV, O26-O30 (1.527 Å) for U = 6 eV,
  O18-O22 (1.566 Å) for U = 8 eV and HSE06 // U = 8, O17-O21 (1.541 Å) for
  PBE // U = 8. Unit charges outside -2 / -1 also occur: 0 and -3 (U = 6 eV,
  3 surplus electrons); 0, 0, 0 and -4, -4, -4 (PBE // U = 8, 7 surplus
  electrons).
- (c): Flow: Li = +1 and O = -2 in all 180 frames (Mn and Ru: 4.0 / 4.5 /
  5.0 with zero spread; NMC transition-metal mean 3.526 / 3.895 / 4.263).
  In NMC the standard deviation of the element mean over 20 frames is 0.34-
  0.53 (Ni), 0.36-0.52 (Co) and 0.17-0.19 (Mn), and the time averages
  (x = 0: Ni +1.87, Mn +4.53, Co +0.85) are not the conventional Ni2+, Mn4+,
  Co3+. The greedy rule is steadier for Ni and Co (0.27-0.39 and 0.18-0.27)
  but puts part of the charge on O (mean -1.88 to -1.78) and a little on Li.

## Provenance

All plotted values are printed by the script in a provenance table.

| Panel | Quantity | Source |
|---|---|---|
| a | flow OS per atom | `pyleos.transport.assign_case(load_case(dir))`, defaults (beta = 10, Allen), `cases/{li16o16,li15o16,li14o16}/v3` (EOSCAR + POSCAR) |
| a | greedy OS per atom | OS(int) column of the same EOSCAR files |
| a | charge counting | Li/O counts of the POSCAR files |
| b | O2 units and O-O distances | minimum-image distances from `cases/li15o16/polaron-{u0,u4,u6,u8,sp-pbe,sp-hse}/POSCAR` (labelled u0, u4, u6, u8, sp_pbe, sp_hse below); atoms are paired by repeatedly taking the shortest remaining O-O distance up to 2.4 Å |
| b | flow charge of a unit | sum of the two per-atom values of `assign_case(load_case(dir))` |
| b | greedy charge of a unit | 2 x ZVAL(O) minus the two n_e(assigned) values of the EOSCAR |
| b | facet labels | U from `LDAUU`, hybrid from `LHFCALC` of each INCAR (asserted); "// U = 8" when the POSCAR equals that of `u8` |
| b | R_site, R_elem, site margins (not plotted) | `pyleos.transport.reliability_margins` |
| b | MD check (not plotted) | `cases/li15o16/v3-md2` (EOSCAR with 20 steps, XDATCAR) |
| c | element means per frame, flow and greedy | `results/cathode_md/timeseries_{li2mno3,li2ruo3,nmc}_x{0.00,0.25,0.50}.tsv`, columns `<El>_flow`, `<El>_greedy`, `t_ps` (written by `tools/timeseries.py`, which runs the flow with `--tier2 arbiter`) |
| c | 4 + 2x line | (2 n(O) - n(Li)) / n(TM) from `cases/<system>/md-x*/POSCAR` |

pyleos results are cached in `_cache/fig4_li2o2.json`,
`_cache/fig4_polaron.json` and `_cache/fig4_li2o2_md.json` with the md5 of
every input file and of `pyleos/{parsers,transport,capacity}.py`. No number
in the script is typed by hand.

## Verified against docs/RESULTS.md

| Statement in RESULTS.md | Recomputed | Status |
|---|---|---|
| Li2O2: O mean -1.00 / -0.94 / -0.88 for both flow and greedy (cases/<x>/v31, Foppl cells) | -1.0000 / -0.9375 / -0.8750 for both; no Li away from +1 | agrees |
| U = 0: target O16-O20 back to 1.530 Å, -1.98; hole spread over 6 units (-1.7 to -2.0) | 1.530 Å, -1.98; the eight units are -1.73, -1.73, -1.87, -1.89, -1.89, -1.91, -1.98, -1.98 | agrees. Not stated there: O18, O22, O27, O31 have changed partners in this geometry (O18-O31 1.511 Å, O22-O27 1.516 Å), as for U = 4 eV |
| U = 4: hole not on the target unit; two O2 pairs dissociated | target O16-O20 1.541 Å, -2.00; the short bond is O18-O31, 1.405 Å, -1.21, formed from atoms of two former units (the other two atoms form O22-O27, 1.537 Å) | agrees; RESULTS.md gives no numbers |
| U = 6: hole on O25-O29 (1.377 Å, -1.13); 5 pairs elongated | 1.377 Å, -1.13; five former pairs at 2.14-2.33 Å | agrees. Not stated there: O20 and O21 form a new bond (1.505 Å, -1.72) and O16, O17 are left without a partner |
| U = 8: 1.401 Å, -1.11; other peroxides -1.93 to -2.00; 3 pairs at about 2.2 Å | 1.401 Å, -1.11; -1.93 to -2.00; 2.203, 2.222, 2.286 Å | agrees |
| PBE on the U = 8 geometry: -1.22, about 0.3 e diffuse, ghost = 7 | -1.22; 7 surplus electrons; other bonded units -1.57 to -1.84, elongated pairs -2.24 to -2.32 | -1.22 and 7 agree; "0.3 e" is not reproducible as one number |
| HSE06 on the U = 8 geometry: -1.11, same as PBE+U | -1.11 on the target unit; two other bonded units at -1.50 each, and one electron counted on Li11 | target agrees; the rest of the cell is not the same as PBE+U |
| "In all calculations the O mean of flow and greedy is -0.9375" | flow: -0.9375 in all six; greedy OS(int): -0.9375 in five, -0.875 for HSE06 (Li11 at 0) | differs for HSE06 |
| The unit charges are counted from the greedy occupations | reproduced with the n_e(assigned) column; sums of OS(int) differ (see Notes) | agrees |
| Moving the hole to a neighbouring unit costs about 0.005; R_site 0.3-0.5; R_elem about 91 | smallest site margin 0.0005 / 0.0071 / 0.0058 / 0.0031 / 0.0028 / 0.0053; R_site 0.05 / 0.71 / 0.58 / 0.31 / 0.28 / 0.53 (u0, u4, u6, u8, sp_pbe, sp_hse); R_elem 90.6-90.9 | agrees for the three calculations on the U = 8 geometry; the range over all six is wider |
| The flow cannot say which unit carries the hole | the flow gives -2 to the short-bond unit in all five calculations that have one | agrees; RESULTS.md does not mention the unit charges of 0, -3 and -4 |
| MD of Li15O16: 40 steps, analysed every 2 (20 frames), flow O mean -0.94 in all; data in `li15o16_md` | -0.9375 in all 20 steps, for either possible step-to-geometry map | value agrees; the 20 steps are in `cases/li15o16/v3-md2` (`v3-md`, the former `li15o16_md`, holds one step), which has no INCAR; `v3-md/INCAR` has TEBEG = 400, not room temperature |
| Li2MnO3, Li2RuO3: +4.0 / +4.5 / +5.0, drift 0.000 | 4.0000 / 4.5000 / 5.0000 in all frames | agrees |
| NMC: Li +1, O -2 throughout | agrees | agrees |
| greedy Mn mean 3.984 +- 0.047 | 3.984 for Li2MnO3 at x = 0; s.d. 0.027, half-range 0.031 | mean agrees; +-0.047 not reproduced from the table |
| NMC: Ni fluctuation +-0.9 (x = 0) to +-1.7 (x = 0.50) | element-mean Ni: half-range 0.67 / 1.17 / 1.17, s.d. 0.34 / 0.52 / 0.53 | not reproducible from the tsv files (they hold element means only) |

## Caveats

1. **Panel (b), greedy quantity.** The squares are electron counts, not
   integer oxidation states; the flow circles are integers. This follows
   the table in RESULTS.md and is the only form in which the greedy rule
   gives a clean unit charge (see Notes).
2. **Panel (b), geometries.** The O 2p U elongates peroxide units that
   do not carry the hole (three at 2.20-2.29 Å for U = 8 eV, three plus
   further rearrangement for U = 6 eV) and, for U = 0 and 4 eV, two units
   have exchanged partners. The U = 8 eV geometry is cross-checked by the
   PBE and HSE06 single points; none of the geometries is a relaxed HSE06
   structure.
3. **Panel (b), molecular dynamics not drawn.** The 20 analysed steps of
   `li15o16/v3-md2` carry 76 Wannier functions per spin (41 surplus
   electrons). The flow O mean is -0.9375 in every step, which is fixed by
   Li = +1; the flow charges of the O2 units scatter between -4 and +2. The
   greedy OS(int) mean is also -0.9375 in every step. All O-O distances
   stay at 1.51-1.54 Å, and the step-to-geometry map is not recorded, so
   there is no bond-length signal to plot.
4. **Panel (a), site resolution.** The flow charges of the eight O2 units
   are [-2 x 8] (Li16O16), [-2 x 7, -1] (Li15O16) and [-2 x 6, -1, -1]
   (Li14O16), the charge-counting pattern. Surplus electrons: 0 / 1 / 2;
   |dNe| below 2e-10 in all three. The cells are the Foppl 2x2x1 supercells
   of cases/<x>/v31; the earlier v3 cells had the two O columns swapped
   (8 Li at O-O bond midpoints), see docs/testset_findings.md.
5. **Panel (c), functional.** With U on the transition metals only, all
   holes created by delithiation are assigned to the metals and O stays at
   -2; which element is oxidised depends on the functional (RESULTS.md
   section 5). The frames carry 18-67 surplus electrons and up to 27 stray
   Wannier functions (NMC). The x values of NMC are nominal (7 and 14 of 29
   Li removed, 0.24 and 0.48). The flow values in the tsv files were
   produced with the projection-assisted option (`--tier2 arbiter`), unlike
   panels (a) and (b), which use the default cost.
