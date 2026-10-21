# Fig. 2 - caption draft and provenance

Script: `fig2_benchmark.py` (run `python3 fig2_benchmark.py` in this
directory; `--recompute` rebuilds the caches under `_cache/` first;
`--check-beta` repeats the *β* scan on the current data for the printed
comparison of panel (e), a few minutes). Output: `fig2_benchmark.pdf` (vector)
/ `.png` (300 dpi), page 17.1 cm x 17.9 cm (PDF page 484.7 x 507.1 pt). The
script prints every plotted number with its source file and the cross-checks
against `docs/`; the last run is kept in `_cache/fig2_benchmark.run.log`.

Layout: upper block (a) above (c) on the left, (b) over the full height on the
right; lower block (d) above (e) on the left, (f) on the right. Reference X
below is G. Comas-Vilà, L. Pujal, A. Otero-de-la-Roza, D. Tiana,
J. Contreras-García and P. Salvador, *J. Chem. Theory Comput.*, 2025, **21**,
7075-7086.

## Caption

**Fig. 2** Benchmark of the minimum-cost-flow assignment. Circles and blue:
min-cost flow (this work); squares and orange: per-WF greedy rule of the
original implementation, which assigns every Wannier function (WF) on its own;
triangles and green: min-cost flow without the electronegativity term
(*β* = 0; panel (c) only). All methods use the same WFs.

**(a)** Verdicts on the 72 scored cases of the 74-case benchmark (two cases
have no reference and are not scored), for the legacy set, the two extension
batches and all cases together. For each case, dev is the largest absolute
difference, over the elements present, between the element-mean oxidation
state and its reference value. Pass: dev ≤ 0.25; intermediate:
0.25 < dev ≤ 1; fail: dev > 1. Bars are normalised to the number of scored
cases *n* of each batch; the columns on the right give the numbers of cases.
The min-cost flow passes 64 of 72 cases (89%), the per-WF greedy rule 51 of 72
(71%).

**(b)** dev of both methods for the 23 scored cases that at least one method
does not pass; the other 49 cases pass with both. Top: 15 cases that only the
min-cost flow passes; middle: 2 cases that only the per-WF greedy rule passes;
bottom: 6 cases that neither passes. Within each group the cases are ordered
by the dev of the method that does not pass (in the bottom group, of the
min-cost flow), largest first. Dashed line: pass threshold, dev = 0.25. Tags
give the batch: L1 and L2, the two parts of the legacy set; E1 and E2,
extensions 1 and 2. ReO<sub>3</sub>, As<sub>4</sub>S<sub>4</sub> and
CSe<sub>2</sub> appear twice because the benchmark contains two independent
calculations of each. L = PPh<sub>3</sub>.

**(c)** Per-atom accuracy on the development subset of the clean benchmark (38
cases, 765 atoms with a per-atom reference). Left: fraction of these atoms
that receive the reference oxidation state, and fraction of cases in which
every one of them does, for the three methods; the numbers are counts. Right:
the ten development cases in which the min-cost flow is not fully correct,
with the number of incorrect atoms (incorrect/labelled). In every case the
incorrect atoms belong to a single element; the entry gives that element pair
and the shortest distance from an incorrect atom to another atom of the same
element. A dash marks the nine cases in which every incorrect atom shares a
Wannier centre with an atom of the same element (a homonuclear bond): the
reference divides the shared pair equally, whereas the assignment places the
whole pair on one atom. Dots mark Re<sub>2</sub>, in which the two Re atoms
are 4.80 Å apart and share no Wannier centre.

**(d)** Reliability index *R*<sub>elem</sub> of the min-cost-flow result for
the 72 scored cases, grouped by the verdict of panel (a): pass (*n* = 64), and
intermediate or fail (*n* = 8). *R*<sub>elem</sub> = 100 *m*/(*m* + 1), where
*m* (in Å<sup>2</sup>) is the smallest increase of the total cost when one
atom is forced to take one electron more or fewer and the assignment is solved
again, counting only solutions that change the electron count of an element;
the value of a case is the lowest over its atoms. Each symbol is one case;
crosses mark cases with at least one stray WF (spread > 3 Å<sup>2</sup>, or
centre farther than 1.6 Å from every atom). Black bars: medians (83.0 and
25.55); grey boxes: interquartile ranges. Every case that does not pass is
named, as are the two passing cases below the median of the other group. The
three bracketed cases (*R*<sub>elem</sub> = 25.9, 25.2 and 25.0, top to
bottom) coincide at this scale and share one leader line. Asterisks: cases
whose reference follows another convention than the electronegativity rule
(Ca<sub>3</sub>AuN: Au<sup>3−</sup> versus the electride description;
WAu<sub>12</sub>: 18-electron rule). Tags as in (b), given where a formula
occurs twice among the scored cases. Every case that does not pass has
*R*<sub>elem</sub> ≤ 25.9 or at least one stray WF; the converse does not
hold: 16 of the 64 passing cases contain stray WFs and two have
*R*<sub>elem</sub> < 25.55.

**(e)** Number of the 72 scored cases that pass as a function of the
electronegativity weight *β* in the cost *d*<sup>2</sup> +
*β*(1 − *χ*/*χ*<sub>max</sub>) of assigning a WF to an atom at distance *d*.
The ordinate ends at 72, the number of scored cases. Dashed line: default,
*β* = 10 Å<sup>2</sup>. Filled symbols: the scan stored with the benchmark
data, which gives 62 at the default. Open symbol: the 64 cases passed in the
present benchmark (panel (a)), which the stored scan predates.

**(f)** Element-mean oxidation states from ref. X (abscissa) and from this
work (ordinate) for the 40 solids of ref. X, one point per element per
compound (88 points from 39 solids); the marker area is proportional to the
number of coincident points. Filled symbols: both methods were applied to the
same crystal structure. Open symbols: the cell computed here differs from, or
could not be matched to, the structure analysed in ref. X (K, CsH<sub>3</sub>
and the three HgF<sub>4</sub> cells). Where points of both kinds coincide, the
marker is divided into a filled and an open sector in proportion to their
numbers. Horizontal bars: XeO<sub>3</sub> at 150 GPa, for which
ref. X gives site-resolved states (Xe +6 and +4) without site multiplicities,
so that only an interval for the element mean follows. Every off-diagonal
point is labelled with compound and element; HgF<sub>4</sub>* is the
HgF<sub>4</sub> cell that contains F<sub>2</sub> units. LiN<sub>5</sub> is not
shown because ref. X reports no assignment for it. In SrBe(OH)<sub>4</sub> the
hydroxide group enters as one unit, as in ref. X.

## Numbers behind the caption (all printed by the script)

### Panels (a) and (b)

| Batch | *n* | Min-cost flow, pass / intermediate / fail | Per-WF greedy rule, pass / intermediate / fail |
|---|---|---|---|
| Legacy set (L1 + L2) | 35 | 33 / 1 / 1 | 30 / 2 / 3 |
| Extension 1 (E1) | 12 | 10 / 0 / 2 | 9 / 1 / 2 |
| Extension 2 (E2) | 25 | 21 / 0 / 4 | 12 / 2 / 11 |
| All | 72 | 64 / 1 / 7 (88.9% pass) | 51 / 5 / 16 (70.8% pass) |

Passed by both 49; by the min-cost flow only 15 (10 of them in E2); by the
per-WF greedy rule only 2 (1 in E2); by neither 6.

| Group in (b) | Case (batch) | dev, min-cost flow | dev, per-WF greedy rule |
|---|---|---|---|
| flow only | Rb2CuCl4 (E2) | 0.00 | 8.00 |
| flow only | SrBe(OH)4 (E2) | 0.00 | 6.25 |
| flow only | RbPb (E2) | 0.06 | 5.88 |
| flow only | RhH (E2) | 0.00 | 5.00 |
| flow only | KMgF3 (E2) | 0.00 | 4.00 |
| flow only | WCl3 (E2) | 0.00 | 2.67 |
| flow only | As4S4 (E2) | 0.00 | 2.00 |
| flow only | S4N4 (E2) | 0.00 | 2.00 |
| flow only | ReO3 (L1) | 0.00 | 2.00 |
| flow only | ReO3 (L2) | 0.00 | 2.00 |
| flow only | CrO5 (E1) | 0.00 | 1.00 |
| flow only | H3PO3 (E2) | 0.00 | 1.00 |
| flow only | Li2O2 (L2) | 0.00 | 1.00 |
| flow only | ReF6 (L2) | 0.00 | 1.00 |
| flow only | MoCl2 (E2) | 0.00 | 0.71 |
| greedy only | CSe2 (L1) | 4.00 | 0.00 |
| greedy only | HI3 (E2) | 2.00 | 0.00 |
| neither | WAu12 (E2) | 6.00 | 6.00 |
| neither | WCl5 (E2) | 4.67 | 7.33 |
| neither | CSe2 (E1) | 4.00 | 6.50 |
| neither | RhCl(CO)(PPh3)2 (E2) | 4.00 | 6.00 |
| neither | Ca3AuN (E1) | 2.00 | 7.00 |
| neither | As4S4 (L1) | 1.00 | 2.12 |

### Panel (c), development split (38 cases, 765 labelled atoms)

| Quantity | Min-cost flow | Per-WF greedy rule | Min-cost flow, *β* = 0 |
|---|---|---|---|
| labelled atoms correct | 687 (89.8%) | 585 (76.5%) | 499 (65.2%) |
| cases with all labelled atoms correct | 28 (73.7%) | 27 (71.1%) | 19 (50.0%) |
| cases passing by element means (not plotted) | 38 | 31 | 25 |

Min-cost flow, the ten cases that are not fully correct (shared centre = Wannier
centre whose two nearest atoms are of the same element, nearer distance at
least half of the farther one):

| Case | Batch | Incorrect / labelled | Pair | Shortest like-atom distance (Å) | Shared centres | Incorrect atoms on one |
|---|---|---|---|---|---|---|
| MoCl2 | E2 | 16/72 | Mo-Mo | 2.600-2.606 | 84 | 16/16 |
| K2S2O3 | E2 | 16/56 | S-S | 2.033-2.041 | 8 | 16/16 |
| As4S4 | E2 | 16/32 | As-As | 2.614-2.630 | 8 | 16/16 |
| GaSe | E2 | 8/16 | Ga-Ga | 2.457-2.458 | 4 | 8/8 |
| NaN3 | E1 | 6/8 | N-N | 1.185 | 12 | 6/6 |
| H2S5 | E2 | 5/7 | S-S | 2.059-2.065 | 4 | 5/5 |
| CrO5 | E1 | 4/6 | O-O | 1.393-1.394 | 2 | 4/4 |
| H2S3 | E2 | 3/5 | S-S | 2.059 | 2 | 3/3 |
| H2S2 | E2 | 2/4 | S-S | 2.054 | 1 | 2/2 |
| Re2 | L2 | 2/2 | Re...Re | 4.800 | 0 | 0/2 |

Total 78 incorrect atoms = 765 − 687; 76 of them (nine cases) on homonuclear
bonds, 2 (Re2) not. The plotted distance is the lower end of each range.

### Panel (d)

| Group | *n* | Median *R*<sub>elem</sub> | Interquartile range | Range | Cases with ≥ 1 stray WF |
|---|---|---|---|---|---|
| pass | 64 | 83.0 | 74.58-87.48 | 14.5-100.0 | 16 |
| intermediate or fail | 8 | 25.55 | 24.18-58.55 | 3.6-75.6 | 4 |

| Named case (batch) | Verdict, dev | *R*<sub>elem</sub> | Stray WFs |
|---|---|---|---|
| Ca3AuN* (E1) | fail, 2.00 | 75.6 | 6 |
| WAu12* (E2) | fail, 6.00 | 74.0 | 9 |
| As4S4 (L1) | intermediate, 1.00 | 53.4 | 15 |
| CSe2 (E1) | fail, 4.00 | 25.9 | 0 |
| HI3 (E2) | fail, 2.00 | 25.2 | 0 |
| CSe2 (L1) | fail, 4.00 | 25.0 | 0 |
| RhCl(CO)(PPh3)2 (E2) | fail, 4.00 | 21.7 | 0 |
| WCl5 (E2) | fail, 4.67 | 3.6 | 316 |
| (CH3)3NO (E1) | pass | 21.1 | 0 |
| CH2F2 (L1) | pass | 14.5 | 0 |

### Panel (e)

| *β* (Å<sup>2</sup>) | 3 | 5 | 8 | 9 | 10 | 11 | 12 | 15 | 20 |
|---|---|---|---|---|---|---|---|---|---|
| cases passed, stored scan (plotted) | 55 | 58 | 63 | 62 | 62 | 62 | 62 | 61 | 61 |

Open symbol: 64 = passes of `results/bench74_golden.tsv` at the default.
Ordinate from 48 to 72 (labelled ticks every 8, the top of the frame is all 72
scored cases); the stored scan spans 55-63, that is 76.4-87.5% of the cases.
Status changes relative to *β* = 10 listed in the scan file: at 3, CH2F2, FeO,
(CH3)3NO, H2S2, H2S3, H2S5, RhH pass → fail, ReO3 pass → intermediate,
RhCl(CO)(PPh3)2 fail → pass; at 5, CH2F2, (CH3)3NO, H2S2, H2S3, H2S5 pass →
fail, RhCl(CO)(PPh3)2 fail → pass; at 8, RhCl(CO)(PPh3)2 fail → pass; at 15,
Cr(CO)6 pass → fail; at 20, Cr(CO)6 and (CH3)3NO pass → fail, CSe2 (L1) fail
→ pass.

### Panel (f), 88 species from 39 solids

| | on the diagonal | off the diagonal |
|---|---|---|
| same crystal structure (79) | 73 | 6: Ca3AuN (Ca, Au); XeO3 130 GPa (Xe, O); XeO3 150 GPa (Xe, O; intervals) |
| different or unmatched structure (9) | 6: CsH3 (Cs, H); HgF4 cells a and c (Hg, F) | 3: K; HgF4* = cell b (Hg, F) |

Off-diagonal values (ref. X → this work): Ca3AuN Au −3 → −1 and Ca +2 → +1.33
(sites +2, +2, 0); K +1 → 0 (sites −1 and +1); XeO3 130 GPa Xe +4 → +5 (sites
+4 and +6) and O −1.33 → −1.67; XeO3 150 GPa Xe between +4 and +6 → +6 and O
between −2 and −1.33 → −2; HgF4* Hg +4 → +3 (sites +2 and +4) and F −1 →
−0.75. Per solid (printed, not drawn): 31 agree, 3 differ (Ca3AuN and the two
XeO3 phases), 6 are not comparable (K, the three HgF4 cells, CsH3, LiN5).
Fifteen of the 40 calculations carry a stray-WF warning (ba3si4 36, xeo3_b 9,
ceh4 8, thh4 8, ca3aun 6, k_hp 6, mgxe 5, ybafe2o5_b 4, csh3 3, mgh4 3, sch4 3,
rb3auo 2, wcl4 2, lin5 1, xeo3_a 1), among them all three solids that differ.

## Provenance

| Panel | Source | What is read or computed |
|---|---|---|
| (a), (b) | `results/bench74_golden.tsv` (min-cost flow), `results/bench74_greedy.tsv` (per-WF greedy rule = `OS(int)` column of EOSCAR) | columns `batch`, `case`, `status`, `dev`. The two rows with status `run` (AuClB complex, CuCl4: no reference) are left out. Legacy set = `v1_mol` (L1) + `v1_oxide` (L2), extension 1 = `wave1` (E1), extension 2 = `wave2` (E2). Both tables come from `python3 -m pyleos bench full` (without and with `--greedy`), same dev definition. The script checks that every `status` agrees with its `dev` |
| (c), left | frozen evaluator, development split only: `python3 -m pyleos bench clean --split dev --json paper/_cache/flow.json`; `... --split dev --greedy --json .../nearest.json`; `... --split dev --method pyleos/methods/baseline_geo.py --json .../geo.json` | `summary` block of the three JSON files (the script refuses a file whose `split` is not `dev`, and runs the commands itself when a file is missing). No other split was scored |
| (c), right | the same `flow.json`; `_cache/flow_contacts.json` | counts from `flow.json`. The incorrect sites are located by re-running the frozen plugin `pyleos/methods/baseline.py` through the evaluator's own `scrub` and `site_labels` (`pyleos.bench`); distances from POSCAR with periodic images; Wannier centres from the `wannier90*.wout` files of the case |
| (d) | `results/bench74_golden.tsv` | columns `status`, `R_elem`, `n_stray`; pass = `OK`, intermediate or fail = `~`, `FAIL`. Asterisks are set in the script (`CONVENTION`) from `docs/RESULTS.md`, sections 3 and 8, not from a data column |
| (e) | `results/scan_beta.tsv` (from `tools/beta_scan.py`) | columns `beta`, `n_OK`; 72 = `n_OK + n_mid + n_FAIL` for every *β*, which also sets the upper end of the ordinate and the number in its label; default *β* from the column name `changed_vs_beta10`, checked against the default of `pyleos.transport.assign_case`. Open symbol: number of `OK` rows of `bench74_golden.tsv`; it is drawn only while that number differs from the scan value at the default |
| (e), printed check only | `_cache/beta_scan_check.tsv`, written by `python3 tools/beta_scan.py --out paper/_cache/beta_scan_check.tsv` (`fig2_benchmark.py --check-beta` runs it) | the same scan on the data tree as it is now. Not plotted and not under `results/`; the script prints it next to the stored scan |
| (f), ordinate | POSCAR, EOSCAR and `wannier90*.wout` of `cases/<case>/v3/`: 30 solids with role `comas` in `cases/INDEX.tsv`, and cro3, ca3aun (benchmark set `wave1`), kmgf3, nabf4, srbe_oh_4, cs3auo, rb3auo, sega, wcl4, rb2cucl4 (set `wave2`; 10 solids shared with the benchmark) | `pyleos.parsers.load_case` + `pyleos.transport.assign_case` at the defaults (*β* = 10, Allen electronegativity). Cache `_cache/fig3_leos.json`, keyed by the SHA-1 of the input files and shared with `figS1_comparison.py` |
| (f), abscissa | table `LIT` in the script: the only hand-typed numbers of the figure | copied unchanged from `figS1_comparison.py`; transcribed from `docs/comas_vila_2025_notes.md`, lines 341-365 (journal pages 7078-7084), each row with its line number, page and the basis of the element mean (stated, by charge neutrality, fragment charge per atom, or site average); see `figS1_caption.md` for the full list |

The script is standalone: the code of panels (d)-(f) is copied in, nothing is
imported from the `figS*` scripts. `LIT` and the functions that build the
points of (f) are identical to those of `figS1_comparison.py`.

Independent checks made when the figure was finished (outputs kept outside the
repository): the three evaluator runs repeated on the development split
reproduce the cached JSON files case by case; `assign_case` repeated on the 40
solids reproduces the cached oxidation states, stray and ghost counts; the
counts of (a), (b), (d) and (e) recomputed directly from the TSV files agree
with the printed provenance; the *β* scan repeated on the present data tree
(`--check-beta`) is the check scan described under Caveats.

## Verified against docs/RESULTS.md

- Section 2: pass counts 33/30, 10/9, 21/12, 64/51 of 35, 12, 25, 72 (89% and
  71%); 15 and 2 one-sided cases, 10 and 1 of them in extension 2; the eight
  cases the flow does not pass and their dev values; the six dev pairs of the
  table of corrected cases (SrBe(OH)4, KMgF3, Rb2CuCl4, RhH, RbPb, H3PO3).
- Section 3: As4S4, dev of the greedy rule 2.00 (E2) and 2.12 (L1); 8 As at +1
  and 8 at +3 in the E2 cell.
- Section 4: the 21 element means and the 12 stray counts quoted there are
  reproduced; n_stray and n_ghost of the 10 shared solids equal
  `bench74_golden.tsv`.
- Section 7: 55, 58, 63, 62 (*β* = 9-12) and 62 instead of 64 at the default.
- Section 8: 64 / 8 cases, medians 83.0 / 25.55, *R*<sub>elem</sub> of HI3
  (25.2), RhCl(CO)(PPh3)2 (21.7), CSe2 (25.0, 25.9), WCl5 (3.6), Ca3AuN (75.6),
  WAu12 (74.0).

## Where the data and the documentation disagree

1. The per-atom numbers of panel (c) (687, 585, 499 of 765; 28, 27, 19 of 38)
   are not in `docs/RESULTS.md` or `docs/PLAN.md`; the clean benchmark is not
   described there at all. Their only source is the evaluator output.
2. `docs/PLAN.md`, section "β = 10", states that the result does not change
   between *β* = 8 and 12. The scan file has 63 at *β* = 8 and 62 at 9-12
   (RhCl(CO)(PPh3)2 passes at 8). `docs/RESULTS.md`, section 7, has the
   numbers right.
3. `docs/RESULTS.md`, section 7, describes the changes at *β* = 3 as seven X-H
   and sulfur-chain cases. The scan file lists nine changes: seven pass → fail
   (among them FeO and RhH), ReO3 pass → intermediate, RhCl(CO)(PPh3)2 fail →
   pass. At *β* = 20 the file also lists CSe2 (L1) fail → pass, which the text
   does not mention. The counts agree.
4. The stored *β* scan and the benchmark table are not from the same state of
   the benchmark (62 against 64 at the default); documented in
   `docs/RESULTS.md`, section 7, and shown in panel (e). The reason given
   there (the scan was run before the tungsten chloride structures were
   replaced) is not recorded in any file under `results/`; it was checked in the
   git history: the scan file dates from commit `cad142d` (2026-07-12) and has
   not changed since, and `git diff cad142d 077f1d8 --
   data/v3/reliability_74.tsv` (the path of `results/bench74_golden.tsv` in
   those commits) changes the four tungsten chloride rows only
   (wcl4 FAIL 2.00 → OK, wcl6 FAIL 5.00 → OK, wcl5 FAIL 3.50 → FAIL 4.67,
   wcl3 OK before and after), which takes the count from 62 to 64.
5. `docs/RESULTS.md`, section 4, flags only `hgf4_b` as an unverified match.
   The POSCAR files show that none of the three HgF4 cells is the *I*4/*m*
   structure described in ref. X (cell a: square-planar molecules; cell b:
   F2 units, F-F 1.43 Å; cell c: cubic, 5 atoms), so all three are drawn as
   open symbols. The same section does not say that the K cell
   (`k_hp/v3`, *oP*8, 8 atoms) differs from the cell of ref. X (4 K
   and 2 non-nuclear attractors).
6. `docs/RESULTS.md`, section 4, takes the 10 shared solids from
   `bench74_golden.tsv`; panel (f) recomputes them from the case directories.
   The results are the same.

## Caveats

- dev tests element means only. A case can pass with wrong individual atoms
  (As4S4 in E2: dev = 0.00, 16 of 32 atoms off by one); panel (c) is the
  per-atom test.
- The 72 scored cases are not 72 distinct compounds: ReO3 (L1, L2), As4S4
  (L1, E2), CSe2 (L1, E1) and TiO2 (L2, E1) were each computed twice.
- The case labels are the benchmark names. Per `docs/RESULTS.md` and
  `docs/testset_findings.md`, the HI3 cell contains HI + I2 chains, the
  RhCl(CO)(PPh3)2 cell carries an extra SO2 ligand, and three passing cases
  are not the compound of their name (CO, SnO, Ti2O3). The reasons for the
  eight cases the flow does not pass (disputed references, near-degenerate
  electronegativities, a failed Wannierisation) are in `docs/RESULTS.md`,
  section 2, and belong in the text, not in the figure.
- Panel (c) uses the development split of the clean benchmark
  (`pyleos/sets/clean_set.json`). The clean set keeps 50 of the 74 cases
  after removing those with a missing, disputed or convention-dependent
  reference, suspect structures, a broken Wannierisation and duplicates; the
  38 of them without unlocalised WFs in the archived data form the development
  split. The other 12 are held out and were not scored, so these numbers are
  not a test-set result. On this subset the min-cost flow passes all 38 cases
  by element means, which is why panel (c) cannot be compared directly with
  panel (a). All three methods are scored on the same 38 cases.
- The per-atom reference divides the electrons of a homonuclear bond equally.
  An assignment that moves whole WFs cannot reproduce that, so the nine bonded
  cases are a property of the representation, not of the cost function.
- Re2 is not a bonded dimer in the computed cell. The POSCAR is in Cartesian
  coordinates with a scale factor of 12, which also multiplies the
  coordinates; the atoms end up 4.80 Å apart in the 12 Å box. All 12 WFs are
  centred on the atoms (within 0.02 Å). The +1/−1 split comes from
  distributing seven electron pairs (14 valence electrons) over two equivalent
  atoms; the per-WF greedy rule gives 0 and 0 for this case.
- "Shared Wannier centre" depends on one constant of the script
  (`SHARED_RATIO = 0.5`). In all cases but MoCl2 the distance ratios fall into
  two well-separated groups (at most 0.24, at least 0.69), so the
  classification does not depend on the threshold; in MoCl2 the ratios are
  spread continuously across 0.5 (0.493 and 0.502 are the nearest values).
- The per-WF greedy rule has almost as many fully correct cases as the
  min-cost flow (27 against 28 of 38); the difference between the two methods
  on this subset is in the number of atoms (585 against 687 of 765).
- Panel (d): *R*<sub>elem</sub> measures how stable the solution is for the
  given input; a high value does not mean agreement with the reference
  (Ca3AuN, WAu12), and a low value does not mean disagreement (CH2F2,
  (CH3)3NO). Median and quartiles of the group of eight are descriptive only.
  The bracket is drawn for named cases whose *R*<sub>elem</sub> differ by less
  than 1 (`R_TIE`).
- Panel (e) plots the stored scan, which is two cases behind the benchmark
  table. A check scan of the data tree as it is now
  (`_cache/beta_scan_check.tsv`, see Provenance; not plotted, not under
  `results/`) gives 57, 60, 65, 64, 64, 64, 64, 63, 63 for *β* = 3, 5, 8, 9, 10,
  11, 12, 15, 20: the stored curve shifted by +2 at every *β*, with an
  identical list of cases changing status, and 64 at the default as in panel
  (a). Regenerating `results/scan_beta.tsv` (and the numbers of
  `docs/RESULTS.md`, section 7) would make panels (a) and (e) refer to the
  same benchmark; the open symbol and its label then disappear by themselves,
  and the two caption sentences on filled and open symbols can be dropped.
  The check-scan values are deliberately not in the printed caption.
- The printed caption (e) says only that the stored scan predates the present
  benchmark. If the scan is not regenerated, the reason can be added: "The two
  counts differ because the scan was run before the structures of
  WCl<sub>4</sub> and WCl<sub>6</sub> were corrected." This is stated in
  `docs/RESULTS.md`, section 7, and was checked in the git history (item 4 of
  the previous section), but no file under `results/` records it.
- Panel (e): the ordinate runs from 48 to 72, not from zero. Its upper end is
  the number of scored cases, which the axis label also states, so that the
  change along the scan (55 to 63 of 72) is seen against the full count. On an
  axis from zero the single-case steps (62, 63, 61) would not be resolved at
  this panel height.
- Panel (f) compares element means only. Site-resolved patterns differ in
  several solids whose means agree (XeO at 200 GPa, both XeO2 phases, the
  mixed-valence YBaFe2O5 phase, ScH4 and YH4, GaSe, Rb2CuCl4; listed in
  `figS1_caption.md`). The three HgF4 cells are paired with solids 34-36 of
  ref. X by directory order only. Cs2Pt (c/a 16.6% below experiment) and MgXe
  (lattice parameter from an equation of state) are the same phase with an
  approximate lattice and are drawn filled.
- Panel (f) is drawn without gridlines (`GRID_F = False` in the script), in
  line with the style of the figure set; the integer ticks on all four sides
  locate the points, every off-diagonal point is labelled, and the values are
  listed above. `figS1_comparison.py` keeps a hairline grid at the integer
  oxidation states; `GRID_F = True` restores it here.
- Panel (f): the markers that are entirely open (K and the two
  HgF<sub>4</sub>* points) carry a 0.5 pt white rim, which separates the
  HgF<sub>4</sub>* (F) point at (−1, −0.75) from the large node at (−1, −1)
  that it touches. The legend marker "both kinds coincide" is an example with
  two thirds filled, not a data point.
- Panels (d), (e) and (f) repeat panels (a), (b) of `figS2_reliability.py` and
  panel (a) of `figS1_comparison.py`. If this figure is used, those panels
  should not appear a second time in the supplementary figures.
- Typography: page 17.1 cm wide, text 7-9 pt (Arial), formula labels through
  `style.formula`. Panel labels sit on one grid (two columns, each 0.10 cm
  above its axes).
