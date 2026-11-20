# r3-T-B Fig. 3 – caption draft and provenance

Script: `r3tb_fig3_case_studies.py` (run `python3 r3tb_fig3_case_studies.py`
in this directory; `--recompute` ignores the cache, about 2 min). Output:
`r3tb_fig3_case_studies.pdf` / `.png`, page 17.1 cm x 8.42 cm. Layout: (a)
over the full height on the left, (b) above (c) on the right. Each panel is
a function: `draw_delithiation(ax)`, `draw_vacancy(ax)`, `draw_heldout(ax)`.
The printout of the script (every plotted number with its source) is kept in
`_cache/r3tb/fig3_case_studies.run.log`.

Data: `r3tb_applications.py` (same directory; `python3
paper/r3tb_applications.py` from the repo root) runs the frozen plugin
`pyleos/methods/frozen/r3-T-B-cap-repair.py` on every application case with
an EOSCAR in the repo, as the clean benchmark calls it (`bench.load_plugin`,
`bench.scrub`, `plugin.assign`), together with the min-cost flow and the
greedy column of the EOSCAR. It writes `results/r3tb_applications.tsv` (433
rows, 69 runs), caches the per-atom results in
`_cache/r3tb/applications.json` and keeps its printout in
`_cache/r3tb/applications.run.log`. The plugin reads no electronegativity
table and no oxidation-state table in any of the 69 runs.

Standard of the figure: every calculation shown is free of unlocalised
(stray) Wannier functions. The 18 calculations of the three panels have none
(`len(assign_case(...).stray)` = 0; the script asserts it).

## Caption

**Fig. 3** Case studies of the in-situ electronegativity assignment.
Triangles: in-situ electronegativity (this work); squares: per-WF greedy
rule; circles: min-cost flow; open diamonds and dashed line: charge
counting. (a) Delithiation of Li2O2 (32-atom Li16O16 cell, PBE) with 0, 1
and 2 Li atoms removed. Top and middle: mean oxidation state of O and of Li.
Charge counting takes Li as +1 in a neutral cell, i.e. −n(Li)/16 for O.
Bottom: oxidation state of every O atom; one symbol for each value, with the
number of atoms at that value printed next to it (the symbol grows with that
number). This work assigns −1 to all 16 O atoms of Li16O16; each removed Li
leaves one O2 unit with the charge −1, whose atoms are −1 and 0. (b) O
vacancy in CeO2: oxidation states of the four Ce atoms of the Ce4O8 cell and
of the same cell with one O atom removed (Ce4O7); all O atoms are −2 in both
cells. Removing the O atom reduces two Ce atoms to +3. Right: the four Ce
atoms around the vacancy V_O, viewed along [001]; each is 2.34 Å from the
vacant site. (c) The 13 held-out compounds with homonuclear bonds (134 atoms
with a reference oxidation state, registered before the calculations were
run): percentage of the atoms of each compound that receive the reference
oxidation state. Over all compounds: 131 of 134 atoms with this work, 71 of
134 with the min-cost flow, 93 of 134 with the per-WF greedy rule.

## Numbers behind the panels

### (a) Li2O2 delithiation (`cases/<cell>/v31`, 112 WFs each)

| Quantity | Li16O16 | Li15O16 | Li14O16 |
|---|---|---|---|
| Li removed | 0 | 1 | 2 |
| Mean OS(O), charge counting | −1.0000 | −0.9375 | −0.8750 |
| Mean OS(O), this work | −1.0000 | −0.9375 | −0.8750 |
| Mean OS(O), greedy; min-cost flow | −1.0000; −1.0000 | −0.9375; −0.9375 | −0.8750; −0.8750 |
| OS(Li), this work (greedy and flow the same) | +1 on all 16 | +1 on all 15 | +1 on all 14 |
| O atoms, this work | −1 x16 | −1 x15, 0 x1 | −1 x14, 0 x2 |
| O atoms, greedy | −3 x1, −2 x3, −1 x8, 0 x3, +1 x1 | −3 x2, −2 x3, −1 x3, 0 x8 | −2 x1, −1 x12, 0 x3 |
| O atoms, min-cost flow (not drawn) | −2 x4, −1 x8, 0 x4 | −2 x7, −1 x1, 0 x8 | −2 x4, −1 x6, 0 x6 |
| Charges of the eight O2 units, this work | −2 x8 | −2 x7, −1 x1 | −2 x6, −1 x2 |
| Same, min-cost flow; greedy (not drawn) | −2 x8; −2 x8 | −2 x7, −1 x1; −3 x2, −2 x3, −1 x3 | −2 x6, −1 x2; −2 x6, −1 x2 |

Symbol sizes in the bottom row: marker size 0.62 + 0.34 √n times the base
size, n = number of atoms; the printed number is n.

### (b) O vacancy in CeO2 (`cases/{ce4o8,ce4o7}/v3`, 96 WFs each)

| Quantity | Ce1 | Ce2 | Ce3 | Ce4 | Ce(+3) | O |
|---|---|---|---|---|---|---|
| Ce4O8, this work | +4 | +4 | +4 | +4 | 0 | −2 x8 |
| Ce4O7, this work | +4 | +3 | +4 | +3 | 2 | −2 x7 |
| Ce4O7, min-cost flow (not drawn) | +3 | +4 | +3 | +4 | 2 | −2 x7 |
| Ce4O7, greedy (not drawn) | +4 | +3 | +3 | +4 | 2 | −2 x7 |
| Distance to the vacant O site (Å) | 2.343 | 2.343 | 2.343 | 2.343 | | |

Reference: two electrons per O vacancy with O at −2, i.e. two Ce(+3). The
sketch is a projection along [001]; the dashed lines are 1.91 Å long on the
page, the label gives the distance in three dimensions.

### (c) Held-out compounds with homonuclear bonds (`cases/<name>/v3`)

Atoms that receive the reference oxidation state, of the labelled atoms
(labels: `pyleos/sets/test_set.json`).

| Compound | Labelled atoms | This work | Min-cost flow | Greedy |
|---|---|---|---|---|
| FeS2 | 12 | 12 | 4 | 12 |
| Na2S2 | 8 | 8 | 4 | 8 |
| BaS3 | 8 | 8 | 2 | 0 |
| Li2O2 | 8 | 8 | 4 | 8 |
| BaO2 | 3 | 3 | 1 | 0 |
| H2O2 | 16 | 16 | 8 | 16 |
| N2H4 | 12 | 12 | 8 | 12 |
| CaC2 | 3 | 3 | 1 | 3 |
| Hg2Cl2 | 4 | 4 | 2 | 4 |
| GaS | 8 | 8 | 4 | 8 |
| InSe | 4 | 1 | 1 | 2 |
| NaSi | 16 | 16 | 8 | 16 |
| S8 | 32 | 32 | 24 | 4 |
| Total | 134 | 131 | 71 | 93 |
| Compounds with every atom correct | 13 | 12 | 0 | 9 |

InSe with this work: both In at +1 (reference +2), one Se at 0 (reference
−2). The totals of this work and of the min-cost flow equal the frozen
one-shot result in `results/testset_results_summary.tsv` (131/134 and 12/13;
71/134 and 0/13); the 13 rows of this work equal those of
`results/testset_results_cases.tsv`. The greedy score is computed here from
the OS(int) column of the EOSCAR files against the same labels.

## Provenance

| Panel | Quantity | Source |
|---|---|---|
| all | this work | `plugin.assign(bench.scrub(load_case(dir)))` with `plugin = bench.load_plugin("pyleos/methods/frozen/r3-T-B-cap-repair.py")`, after `bench.static_check` and `bench.install_trackers` |
| all | greedy | OS(int) column of the EOSCAR of the same run |
| all | stray WFs | `len(pyleos.transport.assign_case(load_case(dir)).stray)` |
| a | charge counting | Li and O counts of the POSCAR files: Li +1, O −n(Li)/n(O) |
| a | O2 units | O atoms paired by repeatedly taking the shortest remaining O–O distance up to 2.4 Å (as in `fig4_delithiation.py`) |
| b | vacancy position and distances | the O atom of the Ce4O8 POSCAR without counterpart in the Ce4O7 POSCAR; minimum-image vectors to the four Ce |
| c | this work, min-cost flow | `bench.score_case` on the 13 cases of `bench.test_cases("homonuclear")`, as `python3 -m pyleos bench test --set homonuclear --method <plugin>` does (flow: `pyleos/methods/baseline.py`) |
| c | greedy | OS(int) of the EOSCAR against `bench.site_labels` |
| c | reference labels | `pyleos/sets/test_set.json`, `test_homonuclear` |

## Tested with r3-T-B, not in the figure

All rows are in `results/r3tb_applications.tsv`.

### MH4 hydrides (`cases/<name>/v3`)

Four of the eight calculations contain stray WFs. `draw_hydrides(ax)` in the
script draws this comparison; the figure does not call it. H values are per
atom of one formula unit. Literature: Comas-Vilà et al., *J. Chem. Theory
Comput.* 2025, 21, 7075–7086, table `LIT_MH4` of `r3tb_applications.py`
(`docs/comas_vila_2025_notes.md`, lines 352–354).

| Compound | Cell | Shortest H–H (Å) | M: this work / lit. / greedy | H atoms, this work | H, lit. | Mean H: this work / lit. | H atoms, min-cost flow | Stray WFs / N_WF | \|dNe\| |
|---|---|---|---|---|---|---|---|---|---|
| MgH4 | I4/mmm, sites H_b and H_a | 0.757 (H_a–H_a) | +2 / +2 / +2 | H_b −1, −1; H2 unit 0, 0 | H_b −1; H_a 0 | −0.50 / −0.50 | H_b −1, −1; H2 unit −1, +1 | 3 / 8 | 1.8e-3 |
| CaH4 | 5 atoms, one H site | 1.366 | +2 / +2 / +4 | pairs (−1, −1) and (−1, +1) | H_b −1; H_a 0 | −0.50 / −0.50 | pairs (−1, −1) and (−1, +1) | 0 / 8 | 3.3e-4 |
| SrH4 | I4/mmm | 0.765 (H_a–H_a) | +2 / +2 / +10 | H_b −1, −1; H2 unit 0, 0 | H_b −1; H_a 0 | −0.50 / −0.50 | H_b −1, −1; H2 unit −1, +1 | 0 / 14 | 2.5e-11 |
| ScH4 | 5 atoms, one H site | 1.686 | +3 / +3 / +5 | pairs (−1, −1) and (−1, 0) | H_b −1; H_a −1/2 | −0.75 / −0.75 | pairs (−1, −1) and (−1, 0) | 3 / 9 | 4.4e-5 |
| YH4 | 5 atoms, one H site | 1.806 | +3 / +3 / +11 | pairs (−1, −1) and (−1, 0) | H_b −1; H_a −1/2 | −0.75 / −0.75 | pairs (−1, −1) and (−1, 0) | 0 / 8 | 8.2e-2 |
| ZrH4 | 5 atoms, one H site | 1.776 | +4 / +4 / +10 | −1, −1, −1, −1 | all −1 | −1.00 / −1.00 | −1, −1, −1, −1 | 0 / 9 | 1.2e-1 |
| CeH4 | I4/mmm | 1.566 (H_a–H_a) | +4 / +4 / +9 | −1, −1, −1, −1 | all −1 | −1.00 / −1.00 | −1, −1, −1, −1 | 8 / 20 | 1.1e-5 |
| ThH4 | I4/mmm | 1.664 (H_a–H_a) | +4 / +4 / +12 and +8 | −1, −1, −1, −1 | all −1 | −1.00 / −1.00 | −1, −1, −1, −1 | 8 / 18 | 1.4e-3 |

Metal: this work equals the literature in 8 of 8 compounds. Mean OS(H):
equal in 8 of 8. Every H atom equal to the literature value of its site: 5
of 8 (MgH4, SrH4, ZrH4, CeH4, ThH4).

### Other cases

Values for methane: initial / transition / final state.

| Case study | Reference | r3-T-B | Min-cost flow |
|---|---|---|---|
| Methane on Co4/CeO2, transferred H (atom 60) | +1 / −1 / −1 (fig. 5(a)) | +1 / 0 / −1 in both data sets | +1 / −1 / −1 |
| Methane, C | −4 / −4 / −4 | −4 / −3 / −4 (`v3`); −4 / −2 / −4 (`v31-mag6`, `v31-mag6`, `v31-mag10`) | −4 / −4 / −4 |
| Methane, Co4 total | step +2, then constant | +3 / +1 / +5 (`v3`); +2 / +2 / +5 (reruns) | 0 / +2 / +2; +2 / +4 / +4 |
| Hole in Li15O16, charge of the short-bond O2 unit (five `polaron-*` runs that have one) | −1 | −1 for U = 4 eV (1.405 Å), U = 8 eV (1.401 Å), HSE06 // U = 8 (1.401 Å); −2 for U = 6 eV (1.377 Å) and PBE // U = 8 (1.401 Å) | −2 in all five |
| Li15O16 MD, `v3-md2`, 20 steps (18 to 21 stray WFs of 152 per step) | O2 units −2 or −1 | mean OS(O) −0.9375 and Li +1 in every step; O2 units −3 to 0 | O2 units −4 to +2 |
| Fe3O4, `v3-verwey` | Fe +2 x8, +3 x16 | Fe +1 x1, +2 x6, +3 x17 | Fe 0 x7, +1 x4, +2 x4, +4 x1, +6 x8 |
| Cathode MD, `cases/{li2mno3,li2ruo3,nmc}/md-*` | – | not tested: no EOSCAR in the repo | – |

Methane, Ce and O with r3-T-B: three Ce at +3 in seven of the nine runs,
four in two; single Ce atoms at +2, +5 or +6, or O atoms at −1, occur in six
of the nine runs. `li15o16/v3-u` gives the same result as `li15o16/v31`. The
`v3` runs of the three Li–O cells use the structures registered as wrong in
`cases/STRUCTURES.tsv`; `li15o16/v3-md` has 88 stray WFs of 152. Of the 22
other runs of the Comas-Vilà set, the element means of r3-T-B equal those of
the min-cost flow in 21 (not `xeo3_b`, 9 stray WFs) and the element values
of Comas-Vilà et al. in 17 of the 21 runs for which the notes give one (not
`hgf4_b`, `k_hp`, `xeo3_a`, `xeo3_b`).
