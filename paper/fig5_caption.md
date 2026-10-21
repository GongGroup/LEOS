# Fig. 5 – caption draft and provenance

Script: `fig5_case_studies.py` (run `python3 fig5_case_studies.py` in this
directory; `--recompute` ignores `_cache/fig5_*.json`, about 1 min). Output:
`fig5_case_studies.pdf` / `.png`, page 17.1 cm x 11.8 cm. Layout: (a) fills
the left column, (b) above (c) in the right column. The tables printed by the
script are kept in `_cache/fig5_case_studies.run.log`.

## Caption

**Fig. 5** Three cases in which the oxidation states are not fixed by the
composition. Blue circles: min-cost flow (this work); orange squares: per-WF
greedy rule; grey plus signs: electron count weighted by the occupations of
the Wannier functions (WFs), without the cap on the electron number (see
below); diamonds: Comas-Vilà et al.; black dots in the bottom rows of (a) and
(b): distances. (a) C–H activation of methane on a Co4 cluster on CeO2(111)
(Ce18O36Co4CH4, PBE+U). Top: the cluster, the four surface O atoms bound to
it and the molecule in the initial, transition and final state; the other
atoms are omitted. H_t (atom 60) is the H atom that is transferred to the
cluster; the dotted line is the C···H_t contact of the transition state.
Below, from top to bottom: total oxidation state of the Co4 cluster treated
as one fragment; number of Ce atoms at +3 (of 18); oxidation states of H_t
and of C; distances from H_t to C and to the nearest Co atom. Filled symbols
and solid lines: Wannierisation with the number of WFs set for each spin
channel, in the magnetic ground state; open symbols and dotted lines: the
same number of WFs in both spin channels. For the plus signs the same flow
problem is solved without the cap on the total electron number and every WF
is counted with its occupation; in the second row they give the number of Ce
atoms that hold more than 8.5 electrons. The greedy rule is not drawn in (a);
its values are listed in Table Sx (ESI). (b) The tetrahydrides MH4. Top:
oxidation state of the metal; for the greedy rule only the metal is drawn,
where it differs from the flow, with one square for each of the two Th atoms
of the ThH4 cell (+12 and +8). Middle: oxidation state of each of the four H
atoms of one formula unit; two circles joined by a line are the two atoms of
an H2 unit (H–H below 1.0 Å). The diamonds are the literature values for the
two H sites of the I4/mmm structure, H_b (left) and H_a (right); the two
sites are labelled once, above CeH4. Bottom: shortest H–H distance in the
cell. Grey bands: cells in which the four H atoms are equivalent, which is
not the two-site structure of the literature. There the four atoms form two
equivalent pairs (H–H 1.37 to 1.81 Å, no H2 unit), drawn as the left and the
right two circles with the pair of lower total oxidation state on the left;
this order is a drawing convention and does not assign the pairs to H_b and
H_a. † unlocalised WFs present: 3 of 8 WFs (MgH4), 3 of 9 (ScH4), 8 of 20
(CeH4) and 8 of 18 (ThH4); ‡ more than 0.01 electron missing in the Wannier
representation. (c) O vacancy in CeO2: oxidation states of the four Ce atoms
of the Ce4O8 cell and of the same cell with one O atom removed (Ce4O7).
Right: the four Ce atoms around the vacancy V_O, viewed along [001]; each is
2.34 Å from the vacant site (the dashed lines are the projections of these
distances).

## What the data support and what they do not

- (a) Supported by both Wannierisations, by the two replica runs and by the
  WF occupations: on C–H cleavage the Co4 cluster loses two electrons, the
  transferred atom H_t (atom 60, written H60 in the tables and sections
  below) goes from +1 to -1, C stays at -4, the other three H stay at +1 and
  all O at -2. H60 is also the only H atom that moves
  (C–H 1.10 -> 2.35 -> 2.92 Å, Co–H 2.82 -> 1.60 -> 1.66 Å).
- (a) Not supported: the absolute oxidation state of Co4 and the number of
  Ce(+3). The flow gives 0 / +2 / +2 with no Ce(+3) for one data set and
  +2 / +4 / +4 with two Ce(+3) for the other; the WF occupations give +2.6 to
  +3.0 / +4.5 to +4.8 / +4.5 to +4.9 with three Ce(+3) for both. See
  "Methane: what differs between the two data sets".
- (b) Supported: the metal oxidation state (+2, +2, +2, +3, +3, +4, +4, +4)
  and the mean H oxidation state (-0.50, -0.50, -0.50, -0.75, -0.75, -1, -1,
  -1) equal the literature values for all eight compounds; for MgH4 and SrH4
  the H2 unit is neutral as a whole and the other two H are -1; for ZrH4,
  CeH4 and ThH4 all H are -1.
- (b) Not supported: values of single atoms inside a unit of equivalent
  atoms. The two atoms of an H2 unit come out as -1 and +1 (R_site 10.7 for
  MgH4, 0.0 for SrH4), and in CaH4, ScH4 and YH4 one of two equivalent H–H
  pairs carries the odd charge (R_site 0.0, 0.5, 0.0). For CaH4, ScH4, YH4
  and ZrH4 the comparison is also not made on the structure of the
  literature.
- (b) Qualified: the agreement for MgH4, ScH4, CeH4 and ThH4. In these four
  cells 3 of 8, 3 of 9, 8 of 20 and 8 of 18 WFs (33 to 44 %) are flagged as
  unlocalised; `CONTEXT.md` defines this flag as a sign that the
  Wannierisation should be redone. For CeH4 and ThH4 the flow equals the
  literature on every atom, but on Wannierisations that carry this flag.
- (c) Supported: two Ce(+3) and two Ce(+4) after removal of one O, with O at
  -2 and the electron number conserved. Once O is -2 this count follows from
  the composition.
- (c) Not supported: which Ce atoms are reduced. The four Ce atoms are
  equivalent in this cell; the flow reduces Ce1 and Ce3, the greedy rule Ce2
  and Ce3, and the WF occupations put the same charge on all four.

## Numbers behind the panels

### (a) Methane on Co4/CeO2, 500 valence electrons

Columns: initial / transition / final state. Atoms: Ce 1–18, O 19–54, Co
55–58, H 59–62, C 63. "No cap" = the same flow with every WF shipped.

| Quantity | N_WF set per spin channel (filled) | Equal N_WF in both channels (open) |
|---|---|---|
| Directory under `cases/` | `ce18o36co4ch4_is/v31-mag6`, `ce18o36co4ch4_ts/v31-mag6`, `ce18o36co4ch4_fs/v31-mag10` | `ce18o36co4ch4_is/v3`, `ce18o36co4ch4_ts/v3`, `ce18o36co4ch4_fs/v3` |
| N_WF, spin 1 / spin 2 | 253/248, 252/249, 253/248 | 253/253, 254/254, 253/253 |
| Sum of WF occupations, spin 1 / spin 2 | 253.00/247.00, 252.00/248.00, 253.00/247.00 | 253.00/246.98, 254.00/246.00, 253.00/246.81 |
| Moment (μB) | 6.00, 4.00, 6.00 | 6.02, 8.00, 6.19 |
| NELECT(Wannier) | 499.999, 499.998, 499.999 | 499.979, 499.990, 499.807 |
| \|dNe\| | 6.3e-4, 1.5e-3, 1.1e-3 | 2.1e-2, 9.6e-3, 1.9e-1 |
| Surplus electrons (N_WF - 500, "ghost") | 1, 1, 1 | 6, 8, 6 |
| Stray WFs | 0, 0, 0 | 2, 0, 1 |
| MAGMOM in INCAR | Ce all 0; Co 4 x 1.5 (IS, TS), 4 x 2.5 (FS) | IS, TS: Ce6 -1, Ce7 +1, Ce8 +1, Co 4 x 2.0; FS: not set |
| Flow: electrons on Co4 -> oxidation state | 58 -> +2, 56 -> +4, 56 -> +4 | 60 -> 0, 58 -> +2, 58 -> +2 |
| Flow: Ce at +3 | Ce6, Ce8 in all three | none |
| Flow: H60 / C | +1/-4, -1/-4, -1/-4 | +1/-4, -1/-4, -1/-4 |
| WFs left out by the cap: nearest atom (occupation) | Ce7 (1.00) in all three | IS: Ce6, Ce7, Ce8 (1.00 each), Co57 (0.67, 0.70), Co58 (1.00). TS: Ce3 (0.99, 1.00), Ce9 (1.00, 1.00), Ce6 (0.96), Ce7 (1.00), Ce8 (1.00), Co57 (0.83). FS: Ce6 (0.99), Ce7 (1.00), Ce8 (1.00), Co57 (0.87, 0.93), Co58 (0.83) |
| Least occupied WF: atom, occupation (all spin 2) | Co55 0.751, Co56 0.858, Co55 0.763 | Co57 0.625, Ce9 0.308, Co56 0.539 |
| No cap: Ce atoms with nine WFs, their occupation (e) | Ce6, Ce7, Ce8: 9.00/9.00/9.00, 8.85/9.00/9.00, 8.96/9.00/9.00 | Ce6, Ce7, Ce8: 8.98/9.00/9.00, 8.83/8.99/9.00, 8.93/9.00/9.00 (TS also Ce3 and Ce9 with ten WFs holding 8.12 and 8.06 e) |
| No cap: WFs on Co4 -> oxidation state | 58 -> +2, 56 -> +4, 56 -> +4 | 63 -> -3, 59 -> +1, 61 -> -1 |
| No cap: WF occupation on Co4 -> oxidation state (plus signs) | 57.03 -> +2.97, 55.25 -> +4.75, 55.11 -> +4.89 | 57.41 -> +2.59, 55.49 -> +4.51, 55.48 -> +4.52 |
| No cap: occupation missing on Co4 / Ce / O / CH4 (e) | 0.97/0.00/0.03/0.00, 0.75/0.16/0.08/0.01, 0.89/0.04/0.04/0.03 | 5.59/0.03/0.40/0.01, 3.51/3.99/0.47/0.04, 5.52/0.08/0.48/0.11 |
| No cap: H60 / C by occupation | +1.00/-4.00, -0.99/-4.00, -1.00/-3.97 | +1.00/-3.99, -0.97/-3.99, -0.98/-3.91 |
| WFs nearest to Ce6, Ce7, Ce8, spin 1 + spin 2 | 4+5 for each, in all three | IS 4+5 each; TS 4+5, 5+4, 5+4; FS 4+6, 4+5, 4+5 |
| Greedy: Ce not at +4 | IS: Ce5 +5, Ce6 +3, Ce7 +2, Ce8 +3. TS: Ce6, Ce7, Ce8 +3. FS: Ce6, Ce7, Ce12 +3 | IS, TS: Ce6, Ce7, Ce8 +3. FS: Ce6, Ce7, Ce8 +3, Ce14 +5 |
| Greedy: electrons on Ce beyond 18 x 8 | 2.97, 2.80, 2.92 | 2.68, 2.72, 2.61 |
| Greedy: sum over the four Co | +3, +7, +7 | +3, +8, +7 |
| Greedy: H60 / C | +1/-4, -3/-4, -3/-4 | +1/-4, -4/-4, -2/-6 |
| R_elem / R_site with the Co4 fragment (`cli.py`) | 21/2, 22/0, 23/21 | 21/21, 22/22, 23/23 |
| C–H60 / shortest Co–H60 (Å) | 1.101/2.824, 2.350/1.599, 2.915/1.658 | same geometries |
| Nearest Co of Ce6, Ce7, Ce8 (Å) | 3.00/3.94/4.19, 2.87/3.90/4.25, 2.99/4.13/3.87 | same geometries |

Replica runs with a second initial moment (`ce_is_mag10`, `ce_ts_mag10`)
reproduce the filled series: Co4 +2 / +4, Ce6 and Ce8 at +3, moments 6.00 /
4.00, occupation-weighted Co4 +2.97 / +4.73. `ce_fs_mag6_4gpu` holds no
EOSCAR and a truncated `wannier90.2.wout`; it was not analysed.

### (b) MH4 hydrides: this work and Comas-Vilà et al.

H values are per atom of one formula unit. Literature: M and H_b per atom,
H_a per atom of the H_a–H_a unit. Stray WFs are given as a count out of all
WFs of the cell (spin-unpolarised runs, one WF per electron pair).

| Compound | Cell | Shortest H–H (Å) | M: flow / lit. / greedy | H, flow | H, lit. | Mean H: flow / lit. | Greedy, all H of the cell | Stray WFs / N_WF | Surplus e | \|dNe\| | R_elem / R_site |
|---|---|---|---|---|---|---|---|---|---|---|---|
| MgH4 † | I4/mmm, Z = 2, sites H_b and H_a | 0.757 (H_a–H_a) | +2 / +2 / +2 | H_b -1, -1; H2 unit -1, +1 (total 0) | H_b -1; H_a–H_a neutral H2 | -0.50 / -0.50 | -1, -1, -1, -1, -1, +1, +1, -1 | 3 / 8 | 4 | 1.8e-3 | 83.2 / 10.7 |
| CaH4 | 5 atoms, one H site | 1.366 | +2 / +2 / +4 | pairs (-1, -1) and (-1, +1) | as MgH4 | -0.50 / -0.50 | -1, -1, -2, 0 | 0 / 8 | 2 | 3.3e-4 | 91.5 / 0.0 |
| SrH4 | I4/mmm, Z = 2 | 0.765 (H_a–H_a) | +2 / +2 / +10 | H_b -1, -1; H2 unit -1, +1 (total 0) | as MgH4 | -0.50 / -0.50 | -1, -9, -5, -5, 0, 0, 0, 0 | 0 / 14 | 0 | 2.5e-11 | 92.6 / 0.0 |
| ScH4 † | 5 atoms, one H site | 1.686 | +3 / +3 / +5 | pairs (-1, -1) and (-1, 0) | H_b -1; H_a -1/2 each | -0.75 / -0.75 | -2, 0, -2, -1 | 3 / 9 | 3 | 4.4e-5 | 85.4 / 0.5 |
| YH4 ‡ | 5 atoms, one H site | 1.806 | +3 / +3 / +11 | pairs (-1, -1) and (-1, 0) | as ScH4 | -0.75 / -0.75 | -3, -3, -4, -1 | 0 / 8 | 1 | 8.2e-2 | 91.0 / 0.0 |
| ZrH4 ‡ | 5 atoms, one H site | 1.776 | +4 / +4 / +10 | -1, -1, -1, -1 | all -1 | -1.00 / -1.00 | -4, -1, -4, -1 | 0 / 9 | 2 | 1.2e-1 | 81.7 / 81.7 |
| CeH4 † | I4/mmm, Z = 2 | 1.566 (H_a–H_a) | +4 / +4 / +9 | -1, -1, -1, -1 | all -1 | -1.00 / -1.00 | -4, -4, -1, -1, -1, -1, -3, -3 | 8 / 20 | 8 | 1.1e-5 | 84.4 / 84.4 |
| ThH4 † | I4/mmm, Z = 2 | 1.664 (H_a–H_a) | +4 / +4 / +12 (Th9) and +8 (Th10) | -1, -1, -1, -1 | all -1 | -1.00 / -1.00 | -1, -3, -3, -3, -3, -1, -3, -3 | 8 / 18 | 4 | 1.4e-3 | 82.3 / 82.3 |

H_b–H_b distances in the four I4/mmm cells: 2.023, 2.372, 1.926, 1.946 Å.
Site reliability of the H atoms in MgH4: 67 (H_b) and 11 (H_a).

Order of the two pairs in the one-site cells (CaH4, ScH4, YH4, ZrH4): the
pair with the lower total oxidation state is listed first and drawn on the
left (`load_hydrides`, sort key (sum, values)); within a pair the lower value
is drawn first. The left pair is not H_b and the right pair is not H_a; the
two pairs are equivalent. Taking the two pairs by lowest atom index instead
gives the same left/right order for these four cells.

Literature values (the only hand-typed numbers; table `LIT_MH4` of the
script): G. Comas-Vilà et al., *J. Chem. Theory Comput.* 2025, 21,
7075–7086, as transcribed in `docs/comas_vila_2025_notes.md`, lines 352
(Mg, Ca, Sr; p. 7081), 353 (Sc, Y; pp. 7081–7082) and 354 (Zr, Ce, Th;
p. 7082).

### (c) O vacancy in CeO2

| Quantity | Ce1 | Ce2 | Ce3 | Ce4 |
|---|---|---|---|---|
| Distance to the vacant O site (Å) | 2.343 | 2.343 | 2.343 | 2.343 |
| O neighbours in Ce4O7 (all at 2.343 Å) | 7 | 7 | 7 | 7 |
| Ce4O8: flow / greedy / occupation-weighted | +4 / +4 / +4.00 | +4 / +4 / +4.00 | +4 / +4 / +4.00 | +4 / +4 / +4.00 |
| Ce4O7: flow | +3 | +4 | +3 | +4 |
| Ce4O7: greedy | +4 | +3 | +3 | +4 |
| Ce4O7: electrons counted by the greedy rule | 8.47 | 8.58 | 8.53 | 8.43 |
| Ce4O7: occupation-weighted | +3.42 | +3.42 | +3.42 | +3.43 |

| Electron count | Ce4O8 | Ce4O7 |
|---|---|---|
| N_WF, spin 1 / spin 2 | 48 / 48 | 48 / 48 |
| Sum of WF occupations, spin 1 / spin 2 | 48.00 / 48.00 | 45.00 / 45.00 |
| NELECT(Wannier); \|dNe\| | 96.000000; 2.2e-12 | 90.000000; 8.7e-9 |
| Surplus electrons; stray WFs | 0; 0 | 6; 0 |
| O atoms, flow and greedy | all -2 | all -2 |
| Sum of the flow oxidation states | 0 | 0 |
| R_elem / R_site | 91.5 / 91.5 | 87.5 / 0.3 |

The two POSCAR files differ only by the missing O atom at (1/4, 1/4, 1/4);
no atom is displaced.

## Methane: what differs between the two data sets

Both sets use the same three geometries. They differ in the Wannierisation
and, for one image, in the magnetic state.

1. Number of WFs. The legacy runs use the same number of WFs in both spin
   channels (253 or 254), 6 / 8 / 6 more than the 500 electrons. The reruns
   set the number per spin channel (253/248, 252/249, 253/248) and carry one
   WF too many, in the minority channel: 248 (249) WFs hold 247.00 (248.00)
   electrons. No WF is empty there; the missing electron is spread over
   Co-centred minority WFs, the least occupied of which holds 0.75 to 0.86
   electrons.
2. Quality flags. Reruns: no stray WF, |dNe| at most 1.5e-3. Legacy: 2 / 0 /
   1 stray WFs, |dNe| = 2.1e-2 / 9.6e-3 / 1.9e-1 (the final state lacks 0.19
   electrons; the warning threshold of the patch is 1e-2).
3. Magnetic state. Reruns: 6.00 / 4.00 / 6.00 μB, as in
   `groundstate_verdict.tsv` (two initial moments converge to the same
   solution for each image). Legacy: 6.02 / 8.00 / 6.19 μB. In the legacy
   transition state the ninth WF of Ce7 and of Ce8 lies in spin channel 1
   instead of 2, which accounts for the difference of 4 μB; this run is not
   in the ground-state moment.
4. What is common. Without the cap the flow assigns nine WFs to each of
   Ce6, Ce7 and Ce8 in all eight analysed runs, holding 8.83 to 9.00
   electrons, and eight WFs to every other Ce atom (except Ce3 and Ce9 in
   the legacy transition state: ten WFs holding 8.12 and 8.06 electrons).
   The legacy INCAR files of the initial and transition state place initial
   moments on exactly these three atoms (Ce6 -1, Ce7 and Ce8 +1). The WF
   occupation on Co4 is 57.0 / 55.2 / 55.1 electrons in the reruns and
   57.4 / 55.5 / 55.5 in the legacy runs.
5. What the cap does. The flow has to leave out as many electrons as there
   are surplus WFs and leaves out the costliest ones; the cost contains
   beta (1 - chi/chi_max), which is 6.98 for Ce and 4.90 for Co. In the
   reruns the one electron is taken from a fully occupied WF at Ce7
   (occupation 1.00). This turns Ce7 from +3 into +4 and leaves all WFs of
   Co4 counted as full electrons (initial state: 58 WFs holding 57.03
   electrons). In the legacy runs the electrons left out are those of Ce6,
   Ce7 and Ce8 (occupations 0.96 to 1.00), four more on Ce3 and Ce9 in the
   transition state, and three / one / three on Co4, where 5.6 / 3.5 / 5.5
   electrons of occupation are missing.
6. Greedy column. It counts occupations and shows the three reduced Ce in
   both sets: Ce6, Ce7, Ce8 at +3 in six of the eight runs; in the rerun of
   the initial state one WF of Ce5 is counted on Ce7 (Ce5 +5, Ce7 +2), in
   the rerun of the final state one WF of Ce8 is counted on Ce12 (Ce12 +3,
   Ce8 +4). The electrons on Ce beyond 18 x 8 are 2.6 to 3.0 in every run.
   It does not agree with the two Ce(+3) of the reruns' flow, nor with the
   zero of the legacy flow.

Conclusions that the data allow:

- The step is robust: +2 on Co4 between the initial and the transition state
  in both sets (flow), +1.78 and +1.92 by occupation; H60 at -1 afterwards
  (-0.97 to -1.00 by occupation); C at -4 (-3.91 to -4.00).
- The absolute values are not. "Co4 0, all 18 Ce +4" (legacy) and "Co4 +2,
  two Ce +3" (reruns) both result from the cap removing occupied Ce
  electrons. The occupations, the greedy column and the initial moments of
  the legacy input indicate three Ce(+3) in every image and a Co4 charge
  near +3 before and near +5 after C–H cleavage.
- Not tested: that a Wannierisation with exactly 253/247 (252/248) WFs
  would make the flow return +3 / +5 / +5 with three Ce(+3). No such run
  exists in `cases/`.
- R_elem is 21 to 23 in all six runs: the flow itself reports that moving
  one electron between elements costs little here.

## Notes for the text

- The greedy rule is the per-atom table that the VASP patch writes to each
  EOSCAR (column OS(int)); it is not a pure nearest-atom rule. In (a) it is
  not drawn; its values are in the table above (rows "Greedy: ..."). For the
  six runs of that table: sum over the four Co +3 / +7 / +7 (reruns) and
  +3 / +8 / +7 (legacy); H60 +1 before the transfer and -3, -3 (reruns), -4,
  -2 (legacy) after it; C -4 except -6 in the legacy final state. Ce6, Ce7
  and Ce8 are at +3 in six of the eight analysed runs. The replica of the
  transition state gives a Co sum of +12 and H60 at -8 (run log only). The
  caption refers to these values as "Table Sx (ESI)"; that table has to be
  made from these rows and its number inserted. In (b) its H values (-9 to
  +1) are not drawn either.
- The plus signs are a diagnostic computed in the figure script, not an
  output of pyleos: `assign_case` is called with the cap removed
  (`nelect_wann = None`), and the occupations of the #WF block of the EOSCAR
  are summed over the WFs assigned to each atom or fragment.
- Style: the plus signs use the marker `MARKER["aux2"]` of `style.py` with
  the grey `C["muted"]`, not with `C["aux2"]`. This departs from the
  one-colour-one-marker mapping of `style.py` on purpose: `C["aux2"]` is an
  orange close to that of the greedy rule, next to which the plus signs are
  drawn in (c), and grey sets the diagnostic apart from the two methods.
- Symbols: circles are the flow in every oxidation-state and count row,
  filled or open. "Filled" and "open" distinguish only the two
  Wannierisations of (a) (circles and plus signs alike). The black circles in
  the key of (a) stand for that distinction; the black dots in the bottom
  rows of (a) and (b) are distances.
- (a) The reaction energy is not part of the figure. In
  `groundstate_verdict.tsv` the image called transition state (-489.554 eV)
  lies below the initial state (-489.458 eV); the final state is at
  -490.234 eV. The text should not quote a barrier from these single points.
- (a) The view shows the atoms at azimuth 110° about the surface normal and
  15° above the surface plane; a stick is drawn for C–H below 1.30 Å, Co–H
  and Co–O below 2.05 Å, C–Co below 2.30 Å and Co–Co below 2.90 Å.
- (b) In CaH4, ScH4 and YH4 the H atoms form two equivalent pairs of mutual
  nearest neighbours (1.37, 1.69, 1.81 Å). One pair comes out as (-1, -1),
  the other as (-1, +1) or (-1, 0). Read as units this is the composition of
  the literature (2 H- plus a neutral H2; 2 H- plus an H2 unit of charge
  -1), but the two pairs are equivalent and the cells contain no H–H
  distance below 1.37 Å. In the figure the pair with the lower total is on
  the left, under the H_b diamond; this placement is a drawing convention
  (see the note under the table of (b)), and in a talk CaH4 must not be read
  as "H_b at -1 and an H2 unit at H_a, as in the literature".
- (b) The site names H_b and H_a are written once, above CeH4, where all
  symbols sit at -1 and leave room; they apply to the diamonds of every
  compound and to the circles of the four cells without a grey band.
- (b) ThH4 has two greedy squares, +12 (Th9) and +8 (Th10), one for each Th
  atom of the cell. They are not joined: they are two atoms, not a range.
  In the other three two-formula-unit cells both metal atoms have the same
  greedy value.
- (c) The sketch is a projection along [001]. The dashed lines are 1.91 Å
  long on the page; the label "4 × 2.34 Å" is the Ce–vacancy distance in
  three dimensions.
- (c) With 0.43 to 0.58 extra electrons on each Ce and no net spin, the
  calculation describes a vacancy whose two electrons are shared by the four
  Ce neighbours, not two localised Ce(+3).

## Provenance

All plotted values are printed by the script in a provenance table.

| Panel | Quantity | Source |
|---|---|---|
| a | Co4 fragment, Ce, H, C (flow) | `pyleos.transport.assign_case(load_case(dir), fragments=[[54, 55, 56, 57]])`, defaults (beta = 10, Allen), i.e. `cli.py case <dir> --fragment 55,56,57,58`; filled: `cases/ce18o36co4ch4_{is,ts}/v31-mag6`, `cases/ce18o36co4ch4_fs/v31-mag10`; open: `cases/ce18o36co4ch4_{is,ts,fs}/v3` |
| a | number of Ce at +3 (flow) | count of Ce atoms with flow oxidation state +3 |
| a | plus signs | the same call with `eoscar.nelect_wann = None`; occupations from the #WF block of the EOSCAR; Co4: 60 minus the summed occupation; Ce: atoms with more than 8.5 electrons |
| a | distances | minimum-image distances from the POSCAR files (identical in both sets, asserted); the transferred H is the H atom with the longest C–H distance in the final state (atom 60) |
| a | views of the site | POSCAR of the three reruns; atoms 55–63 and the O atoms within 2.05 Å of a Co atom in all three images (O20, O22, O23, O25) |
| a | moments, ground state | sums of the #WF occupations per spin; `results/methane_groundstate_verdict.tsv` |
| a | R_elem, R_site (table only) | `python3 pyleos/cli.py case <dir> --fragment 55,56,57,58` for the reruns; `results/legacy_rerun_pyleos.txt` for the legacy runs |
| b | metal and H (flow), R values | `reliability_margins(case, res=assign_case(case))`, `cases/{mgh4,cah4,srh4,sch4,yh4,zrh4,ceh4,thh4}/v3` |
| b | greedy | OS(int) column of the same EOSCAR files |
| b | H sites, pairs, H–H distances | POSCAR; H atoms are grouped by their H–H and M–H neighbour distances, pairs are mutual nearest neighbours |
| b | † and ‡ | † `len(res.stray) > 0` (spread above 3 Å² or centre more than 1.6 Å from every atom); the fractions in the caption are `len(res.stray)` of `len(res.wfs)`; ‡ \|dNe\| above `EOS_NELECT_WARN` of the EOSCAR header (1e-2) |
| b | literature | `LIT_MH4` in the script, from `docs/comas_vila_2025_notes.md` lines 352–354 |
| c | Ce (flow, greedy), R values | `reliability_margins(case, res=assign_case(case))` and OS(int), `cases/{ce4o8,ce4o7}/v3` |
| c | vacancy position and distances | the O atom of the Ce4O8 POSCAR without counterpart in the Ce4O7 POSCAR; minimum-image vectors to the four Ce (2.34 Å); the sketch draws their x and y components (1.91 Å) |

pyleos results are cached in `_cache/fig5_methane.json`,
`_cache/fig5_hydrides.json` and `_cache/fig5_ceo2.json` with the md5 of every
input file and of `pyleos/{parsers,transport,capacity}.py`; `--recompute`
reproduces the three files unchanged. `_cache/fig5_case_studies_SHOW_OCC_False.png`
is the same figure without the plus signs (`SHOW_OCC = False` in the script).

## Verified against docs/RESULTS.md

| Statement in the docs | Recomputed | Status |
|---|---|---|
| Section 6: Co4 electrons 60 -> 58 -> 58, oxidation state 0 -> +2 -> +2 | reproduced from the legacy runs (`cases/ce18o36co4ch4_*/v3`); the reruns of 2026-07-11 give 58 -> 56 -> 56, +2 -> +4 -> +4 | agrees for the legacy data only; the rerun values are not in RESULTS.md |
| Section 6: initial state with 18 Ce at +4; Ce and O unchanged in the other two images | legacy flow: yes. Reruns: Ce6 and Ce8 at +3. WF occupations of all runs: Ce6, Ce7, Ce8 hold 8.83–9.00 electrons | not supported beyond the legacy flow result |
| Section 6: initial state stray 2, ghost 6 | 2 and 6 (legacy); transition and final state: stray 0 and 1, surplus 8 and 6 | agrees; the final state also has \|dNe\| = 0.19, not mentioned |
| Section 6: one H is transferred as H-; CH3(-1) + H-(-1) + Co4(+2) | H60 -1, C -4, three H +1 in both sets | agrees |
| Section 6: C–H 1.10 -> 2.35 -> 2.92 Å; Co–H 1.60 / 1.66 Å | 1.101 -> 2.350 -> 2.915; 1.599 / 1.658 | agrees |
| Section 6: greedy gives C -6 (final), H between -4 and +1 (transition), Ce mean 3.83–3.89 | legacy: C -6, H60 -4, Ce mean 3.83 / 3.83 / 3.89. Reruns: C -4 in all three, H60 -3 in the transition and final state, Ce mean 3.83 | agrees for the legacy data |
| Section 6: +2 on Co57 in the transition state, +1 on Co57 and Co58 in the final state | legacy per-atom flow: Co 0, 0, +2, 0 and 0, 0, +1, +1; reruns: 0, 0, +2, +2 and 0, +1, +1, +2 | agrees for the legacy data |
| Section 6: R of Ce and O 80–92, of the Co atoms 4–22 | per-atom values not recomputed; case values R_elem 21–23 | not checked |
| Section 6: ground state moments 6 / 4 / 6, final state -490.23423 eV for both initial moments | `groundstate_verdict.tsv`; occupations of the reruns give 6.00 / 4.00 / 6.00 | agrees. The oxidation states quoted in the same section come from runs with 6.02 / 8.00 / 6.19 |
| Section 4: Mg, Ca, Sr: M +2 + 2 H- + H2; Sc, Y: M +3 + 3 H- + 1/2 H2; Zr, Ce, Th: M +4 + 4 H- | metal and mean H agree for all eight; an H2 unit (0.76 Å) exists only in MgH4 and SrH4 | agrees as a count; for CaH4, ScH4 and YH4 it is not a site-resolved result |
| Section 4: Sc, Y differ at the site level (literature: 2 H- and 2 H at -1/2) | flow -1, -1, -1, 0 | agrees |
| Section 4: greedy gives +10 / +9 / +10 for the metal in SrH4, CeH4, ThH4 | +10, +9, and +12 / +8 for the two Th atoms; also +4 (Ca), +5 (Sc), +11 (Y), +10 (Zr); only MgH4 is +2 | agrees; incomplete |
| Section 4: exact structures of the MH4 hydrides in `docs/structures_fix/` | Mg, Sr, Ce, Th only (`MANIFEST.tsv`); the Ca, Sc, Y, Zr cells have 5 atoms and one H site | incomplete |
| Section 4: stray WFs ceh4 8, thh4 8, mgh4 3, sch4 3 | 8 of 20, 8 of 18, 3 of 8, 3 of 9 WFs | agrees; the totals are not given there; \|dNe\| of 0.082 (YH4) and 0.120 (ZrH4) is not mentioned |
| Section 5: Ce4O8 [4, 4, 4, 4]; Ce4O7 [3, 3, 4, 4] | [4, 4, 4, 4]; flow [3, 4, 3, 4] in atom order, greedy [4, 3, 3, 4] | agrees as a count |
| Section 5: electron number conserved to three decimals, no warning, for both cells | 96.000000 and 90.000000; Ce4O7 has 6 surplus electrons (reported by pyleos) and R_site 0.3 | the electron number agrees; "no warning" does not hold for Ce4O7 |
| Section 5: DFT+U was used to localise the Ce 4f electrons | no INCAR is stored in `cases/{ce4o8,ce4o7}/v3`; the occupations show 0.43–0.58 extra electrons on each of the four Ce and no net spin | not verifiable; the stored result is not localised |

## Caveats

1. **Panel (a), absolute values.** Neither series is a result that the WF
   occupations confirm (see the section on methane). The figure shows both
   and the occupation-weighted count; the text should claim the step and the
   hydride, not "Co4(0)" or "Co4(+2)".
2. **Panel (a), plus signs.** They are non-integer, depend on the WF-to-atom
   map of the uncapped flow, and include the 0.02 to 0.19 electrons that the
   legacy Wannierisations lack. They are shown because they are the only
   quantity on which the two data sets agree in absolute terms.
3. **Panel (a), what is not drawn.** The two replica runs, the greedy rule
   and the reliability values are in the tables only. Panels (b) and (c)
   draw the greedy rule and (a) does not; the caption says so and points to
   an ESI table that still has to be made.
4. **Panel (b), structures.** MgH4, SrH4, CeH4 and ThH4 are the I4/mmm cells
   of Bi and Zurek (0 GPa, `docs/structures_fix/`). CaH4, ScH4, YH4 and ZrH4
   are 5-atom cells with a single H site (listed as 0 GPa approximations in
   `docs/testset_registry.tsv`); the literature values next to them refer to
   the two-site structure. In these cells the left/right placement of the
   two H pairs follows their total oxidation state, not a site; the
   coincidence of the left pair with the H_b diamond is produced by that
   convention.
5. **Panel (b), quality.** Stray WFs: 3 of 8 in MgH4, 3 of 9 in ScH4, 8 of
   20 in CeH4 and 8 of 18 in ThH4 (38, 33, 40 and 44 %); \|dNe\| = 0.082
   (YH4) and 0.120 (ZrH4); surplus electrons in seven of the eight cells.
   Only CaH4 and SrH4 carry none of the two marks, and CaH4 is a one-site
   cell. The agreement of CeH4 and ThH4 with the literature on every atom
   was obtained with two fifths of the WFs unlocalised and should not be
   quoted without that fraction. The H2 units of MgH4 and SrH4 are
   recognised by distance (below 1.0 Å, threshold `H2_CUT` in the script);
   the flow does not treat them as fragments.
6. **Panel (c), size of the test.** The cell has 11 atoms (CeO1.75), is not
   relaxed, has no net spin and makes all four Ce atoms equivalent
   neighbours of the vacancy. It tests the electron count, not the location
   of the reduced Ce atoms.
