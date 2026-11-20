# r3-T-B Fig. 2 - caption draft and provenance

Script: `r3tb_fig2_benchmark.py` (run `python3 r3tb_fig2_benchmark.py` in this
directory; `--recompute` rebuilds the two caches under `_cache/r3tb/` first,
about 75 s on a laptop). Output: `r3tb_fig2_benchmark.pdf` (vector) / `.png`
(300 dpi), page 17.1 cm x 20.77 cm (PDF page 484.7 x 588.8 pt), and the two
halves used on slides: `_cache/r3tb/fig2_top.png` (panels a-c, 17.1 x 9.73 cm)
and `_cache/r3tb/fig2_bottom.png` (panels d-f, 17.1 x 11.04 cm), 300 dpi. The
script prints every plotted number with its source; the last run is kept in
`_cache/r3tb/fig2_benchmark.run.log`.

Single panels at poster scale:
`python3 r3tb_fig2_benchmark.py --panel margins --scale 3 --out FILE.png`
(panels: `accuracy`, `cases`, `elements`, `chi`, `constants`, `margins`), or
import the module and call `draw_accuracy(ax)`, `draw_cases(ax)`,
`draw_elements(ax)`, `draw_chi(ax)`, `draw_constants(ax)`, `draw_margins(ax)`
inside `with scaled(k):` - each function fetches its own data, and all
offsets, marker sizes and line widths follow `rcParams["font.size"]`.

Layout: upper block (a) across the full width, (b) below it on the left, (c) on
the right; lower block (d) above (e) on the left, (f) over the full height on
the right. Method names, colours and markers are those of `paper/style.py`
(`aux1`, `flow`, `nearest`, `geo`, `aux2`).

## Caption

**Fig. 2** Benchmark of the in-situ electronegativity assignment on clean38,
the 38 cases of the clean benchmark without unlocalised Wannier functions
(WFs): 793 atoms, 769 of them with a per-atom reference oxidation state, 2768
WFs, 34 elements. Down-pointing triangles and purple: in-situ
electronegativity (this work); circles and blue: min-cost flow with the Allen
table (*β* = 10); squares and orange: per-WF greedy rule; up-pointing
triangles and green: min-cost flow without the electronegativity term
(*β* = 0); plus signs and yellow: Lewis-first assignment with the Allen table.
All methods use the same WFs. This work reads no electronegativity table and
no oxidation-state table.

**(a)** Accuracy of the five methods: atoms that receive the reference
oxidation state (of 769), cases in which every labelled atom does (of 38), and
cases whose element-mean oxidation states are within 0.25 of the reference (of
38). Bars are percentages; the numbers are counts. This work: 758/769 atoms,
36/38 cases fully correct, 38/38 cases by element means.

**(b)** Incorrect atoms, in per cent of the labelled atoms, for this work, the
min-cost flow and the per-WF greedy rule in the 14 cases that at least one of
the three does not get fully right; the other 24 cases are fully correct with
all three. The columns give the numbers of incorrect atoms of the three
methods and the number of labelled atoms. Above the dashed line: the 12 cases
in which this work assigns every atom correctly, ordered by the fraction of
incorrect atoms of the min-cost flow, then of the greedy rule. Below: the two
cases in which this work has incorrect atoms (MoCl<sub>2</sub> 9 of 72,
Rb<sub>2</sub>CuCl<sub>4</sub> 2 of 28). Where markers coincide, the smaller
one is drawn on the larger.

**(c)** Element coverage of this work. A tile is filled from the bottom in
proportion to the labelled atoms of the element that are assigned correctly;
the number is the count of labelled atoms. Every labelled atom is correct for
31 of the 34 elements; the others are Cl (79/82), Mo (18/24) and Cu (11/13).
Empty tiles: elements that do not occur in clean38.

**(d)** In-situ electronegativity *χ* = *Z*<sub>eff</sub>/*ρ*<sup>2</sup>
computed by the method against the Allen electronegativity, for the eight
elements whose atoms carry lone pairs at a contested WF. *Z*<sub>eff</sub> is
the valence charge minus the electrons of the atom's own WFs inside the WF
distance minus 5/16; *ρ*<sup>2</sup> is the second moment of the atom's
lone-pair WFs about the nucleus. Symbols: median over the atoms of the
element; bars: interquartile range. Filled: p-block elements, Spearman rank
correlation *r*<sub>s</sub> = 0.96 (*n* = 7); open: Mo, the one d-block
element; all eight elements: *r*<sub>s</sub> = 0.62. The Allen values are
used for this comparison only.

**(e)** Atoms correct on clean38 (of 769, the upper end of the ordinate) when
one of the four constants of the method is varied and the other three keep
their defaults: the screening of a pair partner in *Z*<sub>eff</sub> (default
5/16 electron), the centre offset below which an own WF counts as a
dipole-free shell (1/8 rms radius), the extent of a WF (1 rms radius) and the
numerical tolerance (10<sup>−6</sup>). Triangle and dashed line: default.
Grey band: values for which every one of the 38 cases gives the same result
as at the default: 0.05-0.625, 0.075-0.175, 0.55-1.25 and
10<sup>−12</sup>-10<sup>−3</sup>.

**(f)** Winner/loser ratio of the deciding quantity for the 566 contested WFs
that have a partner atom of another element (a contested WF lies nearer a bond
midpoint than either nucleus). Filled circles: decided by the force
*Z*<sub>eff</sub>/*d*<sup>2</sup> (527 WFs); open triangles: decided by *χ*,
which applies when both ends carry lone pairs (39 WFs). Rows are bond classes,
winner first, ordered by the median ratio; the columns give the numbers of WFs
decided by force and by *χ*. Upper block: one atom takes the pair (466 WFs).
Lower block: atoms of the winning element share the pair (subscript: number
of sharing atoms) and an atom of another element contests it (100 WFs).
Dashed line: ratio 1. Grey slot (∞): *Z*<sub>eff</sub> of the loser is not
positive (61 WFs). Each row is labelled winner-loser. In 462 of the 466
one-owner WFs the pair goes to the element with the larger Allen
electronegativity.

## Numbers behind the caption (all printed by the script)

### Panel (a)

| Method | Atoms correct | Cases fully correct | Element-mean pass | File |
|---|---|---|---|---|
| In-situ electronegativity (this work) | 758/769 (98.6%) | 36/38 (94.7%) | 38/38 | `results/clean38_r3-T-B.tsv` |
| Min-cost flow | 689/769 (89.6%) | 27/38 (71.1%) | 38/38 | `results/clean38_flow.tsv` |
| Per-WF greedy rule | 589/769 (76.6%) | 27/38 (71.1%) | 31/38 | `results/clean38_greedy.tsv` |
| Min-cost flow, *β* = 0 | 497/769 (64.6%) | 19/38 (50.0%) | 25/38 | `results/clean38_geo.tsv` |
| Lewis-first with Allen table | 769/769 (100%) | 38/38 (100%) | 38/38 | `results/clean38_r2-W-A.tsv` |

### Panel (b), incorrect atoms (per cent of the labelled atoms)

| Case | Labelled | This work | Min-cost flow | Per-WF greedy rule |
|---|---|---|---|---|
| Re2 | 2 | 0 | 2 (100%) | 0 |
| NaN3 | 8 | 0 | 6 (75.0%) | 6 (75.0%) |
| H2S5 | 7 | 0 | 5 (71.4%) | 4 (57.1%) |
| CrO5 | 6 | 0 | 4 (66.7%) | 2 (33.3%) |
| H2S3 | 5 | 0 | 3 (60.0%) | 3 (60.0%) |
| As4S4 | 32 | 0 | 16 (50.0%) | 28 (87.5%) |
| H2S2 | 4 | 0 | 2 (50.0%) | 0 |
| GaSe | 16 | 0 | 8 (50.0%) | 0 |
| K2S2O3 | 56 | 0 | 16 (28.6%) | 16 (28.6%) |
| KMgF3 | 5 | 0 | 0 | 3 (60.0%) |
| SrBe(OH)4 | 160 | 0 | 0 | 66 (41.2%) |
| ReF6 | 7 | 0 | 0 | 2 (28.6%) |
| MoCl2 | 72 | 9 (12.5%) | 16 (22.2%) | 26 (36.1%) |
| Rb2CuCl4 | 28 | 2 (7.1%) | 2 (7.1%) | 24 (85.7%) |
| Total (14 cases) | 408 | 11 | 80 | 180 |

The totals equal 769 minus the atoms correct of panel (a). In the 12 cases
above the dashed line the min-cost flow has 62 incorrect atoms and the greedy
rule 130; this work has none.

### Panel (c), atoms correct / labelled per element

O 216/216, H 153/153, Cl 79/82, S 42/42, C 32/32, F 27/27, Mo 18/24, N 19/19,
K 17/17, Sr 17/17, As 16/16, Be 16/16, Cu 11/13, Cr 9/9, B 8/8, Ga 8/8,
Rb 8/8, Se 8/8, I 7/7, Na 7/7, Ti 7/7, Cs 6/6, Ce 4/4, Re 4/4, La 3/3, Si 3/3,
W 3/3, Au 2/2, Mn 2/2, Sn 2/2, Co 1/1, Fe 1/1, Mg 1/1, V 1/1. Sum 758/769;
31 of 34 elements have every labelled atom correct. The incorrect atoms are
Cl 3 and Mo 6 in MoCl2, Cu 2 in Rb2CuCl4. The 24 atoms without a per-atom
label are the B atoms of B6H10, for which the reference defines the element
mean only.

### Panel (d), in-situ *χ* (e Å<sup>−2</sup>)

| Element | Block | Allen | *χ* median | Interquartile range | Range over atoms | Atoms | Bond ends | Cases |
|---|---|---|---|---|---|---|---|---|
| Mo | d | 1.470 | 5.724 | 4.761-6.925 | 3.238-22.330 | 14 | 49 | MoCl2 |
| As | p | 2.211 | 2.373 | 2.338-2.391 | 2.284-2.395 | 16 | 40 | As4S4 |
| I | p | 2.359 | 3.296 | 3.282-3.304 | 3.242-3.367 | 7 | 9 | [Cu5I7]2- |
| Se | p | 2.424 | 2.956 | 2.954-2.958 | 2.945-2.962 | 8 | 24 | GaSe |
| S | p | 2.589 | 3.358 | 3.284-3.667 | 3.265-4.263 | 30 | 47 | As4S4, H2S2, H2S3, H2S5, K2S2O3 |
| Cl | p | 2.869 | 4.901 | 4.838-4.905 | 1.644-4.975 | 42 | 49 | MoCl2, WCl6 |
| O | p | 3.610 | 7.483 | 7.448-8.819 | 6.093-10.401 | 168 | 260 | CO2, MnO2, NaReO4, SiO2, Cr(CO)6, CrO3, CrO5, HNO3, (CH3)3NO, K2S2O3, SrBe(OH)4 |
| F | p | 4.193 | 13.637 | 13.620-13.637 | 13.620-13.956 | 10 | 10 | CH2F2, NaBF4 |

Order by *χ*: As < Se < I < S < Cl < Mo < O < F; order by Allen: Mo < As < I <
Se < S < Cl < O < F. Spearman rank correlation of the medians with the Allen
values: 0.9643 for the seven p-block elements (one exchange, I and Se), 0.6190
for all eight. The *χ* of an atom is the median over the contested WFs of
panel (f) of which it is the winner or the loser ("bond ends": 488 in all);
the value of an element is the median over its atoms.

### Panel (e), atoms correct of 769 (cases fully correct, element-mean pass)

| Constant (name in the plugin) | Default | Grid | Result unchanged for | Outside |
|---|---|---|---|---|
| Pair screening (`PAIR_SCREEN`) | 5/16 | 0-1, step 0.025, and 5/16 (42 values) | 0.05-0.625: 758 (36, 38) | 0: 749; 0.025: 755; 0.65-0.675: 746; 0.7-1: 734 (35, 37) |
| Dipole-free limit (`PARITY`) | 1/8 | 0-0.5, step 0.0125 (41 values) | 0.075-0.175: 758 (36, 38) | 0: 748; 0.0125-0.0375: 754; 0.05-0.0625: 756; 0.1875: 755; 0.2-0.2375: 751; 0.25-0.3125: 750; 0.325-0.35: 749; 0.3625: 748; 0.375: 745; 0.3875-0.4: 742; 0.4125-0.4375: 740; 0.45-0.5: 734 (31, 33) |
| WF extent (`REACH_SIGMA`) | 1 | 0.5-2, step 0.05 (31 values) | 0.55-1.25: 758 (36, 38) | 0.5: 756; 1.3: 752; 1.35: 746; 1.4: 744; 1.45: 743; 1.5: 741; 1.55: 738; 1.6-1.7: 734; 1.75-1.85: 732; 1.9: 713; 1.95-2: 696 (35, 36) |
| Numerical tolerance (`FLOAT_EPS`) | 10<sup>−6</sup> | 10<sup>−12</sup>-10<sup>−2</sup>, decades (11 values) | 10<sup>−12</sup>-10<sup>−3</sup>: 758 (36, 38) | 10<sup>−2</sup>: 760 (37, 38) |

"Result unchanged" means that every one of the 38 cases has the atoms correct
and the element-mean verdict of the default run; the limits are grid values.
Cases that change first: CH2F2 below 0.05 and B6H10 above 0.625 (pair
screening); LaFeO3 below 0.075 and above 0.175 (dipole-free limit); MoCl2
below 0.55 and above 1.25 (WF extent); Rb2CuCl4 at 10<sup>−2</sup>
(tolerance, 28/28 there). The ordinate runs from 690 to 769.

### Panel (f), 566 contested WFs with a heteronuclear partner

Of the 2768 WFs of clean38, 618 are contested and 566 of these have a partner
atom of another element; 527 are decided by force and 39 by *χ* (32 in As4S4,
7 in MoCl2).

One atom takes the pair (466 WFs, 23 classes; 430 by force, 36 by *χ*):

| Class (winner-loser) | WFs | Force | *χ* | Ratio min | Median | Max | Cases |
|---|---|---|---|---|---|---|---|
| S-As | 32 | 0 | 32 | 1.363 | 1.388 | 1.470 | As4S4 |
| C-H | 11 | 11 | 0 | 1.297 | 1.463 | 1.501 | CH2F2, (CH3)3NO |
| O-N | 37 | 37 | 0 | 1.101 | 1.760 | 5.394 | HNO3, (CH3)3NO |
| H-B | 24 | 24 | 0 | 1.684 | 1.769 | 1.952 | B6H10 |
| N-C | 3 | 3 | 0 | 2.212 | 2.214 | 2.214 | (CH3)3NO |
| S-H | 6 | 6 | 0 | 2.230 | 2.245 | 2.290 | H2S2, H2S3, H2S5 |
| O-Re | 4 | 4 | 0 | 3.043 | 3.125 | 3.279 | NaReO4 |
| O-S | 24 | 24 | 0 | 3.148 | 3.436 | 3.635 | K2S2O3 |
| N-H | 12 | 12 | 0 | 4.069 | 4.070 | 4.226 | BH6N |
| O-Cr | 21 | 21 | 0 | 2.290 | 5.296 | 7.317 | CrO3, CrO5 |
| O-H | 72 | 72 | 0 | 5.968 | 6.211 | 10.510 | HNO3, SrBe(OH)4 |
| O-C | 96 | 96 | 0 | 3.494 | 6.475 | 8.738 | CO2, Cr(CO)6 |
| F-C | 2 | 2 | 0 | 6.851 | 6.853 | 6.856 | CH2F2 |
| C-Cr | 24 | 24 | 0 | 7.003 | 7.060 | 7.103 | Cr(CO)6 |
| N-B | 4 | 4 | 0 | 7.079 | 7.079 | 7.079 | BH6N |
| Cl-W | 18 | 18 | 0 | 7.204 | 7.501 | 7.545 | WCl6 |
| Se-Ga | 24 | 24 | 0 | 7.759 | 7.926 | 7.998 | GaSe |
| Cl-Mo | 25 | 22 | 3 | 1.294 | 8.011 | 10.045 | MoCl2 |
| O-Si | 2 | 2 | 0 | 12.819 | 12.819 | 12.819 | SiO2 |
| O-Mn | 4 | 4 | 0 | 16.822 | 16.852 | 16.879 | MnO2 |
| F-B | 8 | 8 | 0 | 17.886 | 17.893 | 17.900 | NaBF4 |
| I-Cu | 9 | 9 | 0 | 52.256 | 65.774 | 70.672 | [Cu5I7]2- |
| Mo-Cl | 4 | 3 | 1 | 3.769 | ∞ | 3.769 (3 at ∞) | MoCl2 |

Like atoms share the pair (100 WFs, 5 classes; 97 by force, 3 by *χ*):

| Class | Sharing atoms | WFs | Force | *χ* | Ratio min | Median | Max (finite) | At ∞ | Case |
|---|---|---|---|---|---|---|---|---|---|
| As2-S | 2 | 8 | 8 | 0 | 12.074 | 12.427 | 12.779 | 0 | As4S4 |
| B2-3-H | 2-3 | 16 | 16 | 0 | 18.796 | 18.956 | 20.637 | 0 | B6H10 |
| S2-H | 2 | 1 | 1 | 0 | 25.480 | 25.480 | 25.480 | 0 | H2S2 |
| Mo2-6-Cl | 2-6 | 67 | 64 | 3 | 1.229 | ∞ | 90.923 | 50 | MoCl2 |
| S2-O | 2 | 8 | 8 | 0 | - | ∞ | - | 8 | K2S2O3 |

Force rule: 527 WFs, finite ratios 1.101-90.923, median over all 6.617, 61 at
∞. *χ* rule: 39 WFs, ratios 1.262-7.807, median 1.389. Ratios below 1.5: 47
WFs; below 2: 98 WFs. Smallest ratio: 1.101, an O-N pair in (CH3)3NO.

## Provenance

| Panel | Source | What is read or computed |
|---|---|---|
| (a) | `results/clean38_r3-T-B.tsv`, `clean38_flow.tsv`, `clean38_greedy.tsv`, `clean38_geo.tsv`, `clean38_r2-W-A.tsv`; each from `python3 -m pyleos bench clean --split dev --method <plugin> --out <file>` (greedy: `--greedy`), plugins as listed in the script (`METHODS`) | sums of `atoms_ok` and `atoms_labelled`, number of rows with `atoms_ok == atoms_labelled`, number of rows with `status == OK`. The script checks that every table lists exactly the development split of `pyleos/sets/clean_set.json` (the clean set without the 12 cases of `splits.test_regression`) and the same number of labelled atoms per case |
| (b) | `results/clean38_r3-T-B.tsv`, `clean38_flow.tsv`, `clean38_greedy.tsv` | `atoms_labelled - atoms_ok` per case, for the cases where it is not zero for at least one of the three methods |
| (c) | `_cache/r3tb/fig2_analysis.json`, block `elements` | the frozen plugin `pyleos/methods/frozen/r3-T-B-cap-repair.py` re-run in process on the evaluator's own case list, `scrub` and `site_labels` (`pyleos.bench`); per case the atoms correct, atoms labelled and the list of wrong sites must equal `results/clean38_r3-T-B.tsv`, otherwise the script stops |
| (d) | the same cache, block `chi_atoms`; `pyleos.capacity.ENEG_ALLEN` (abscissa only) | for every contested WF with a heteronuclear partner, *Z*<sub>eff</sub> and *ρ*<sup>2</sup> of its two ends from the plugin's own shell model (`_Shells.seen_from`, `PAIR_SCREEN`). The plugin's `_analyse` does not return the atoms of the two ends, so the script repeats its selection of winner and loser and requires every record to equal the plugin's note (elements, *Z*<sub>eff</sub>, distances, rule, *χ* ratio) and, per element pair, the plugin's `pair_margins(case)` |
| (e) | `_cache/r3tb/fig2_scan.json` | 125 runs of the frozen evaluator, `python3 -m pyleos bench clean --split dev --method <copy> --json <copy>.json`, on copies of the frozen plugin in which the one line `NAME = ...` is replaced by `NAME = <value>` (one constant per copy, the other three at their defaults). `--recompute` writes the copies to a temporary directory outside the repository and removes it afterwards; `--scan-dir DIR` keeps them in `DIR`. The frozen file is not edited and no copy enters the repository. The copy with the default value must reproduce `results/clean38_r3-T-B.tsv` case by case, otherwise the script stops |
| (f) | `_cache/r3tb/fig2_analysis.json`, block `wfs` | the `notes` of the plugin's `_analyse(case)` (the records behind `pair_margins`): winner, loser, rule, and the ratio `chi_margin` (rule *χ*) or `force_margin` (rule force); the number of atoms that hold the WF from `_analyse()['geo'][i]['hosts']`. The count of 462 compares winner and loser in `pyleos.capacity.ENEG_ALLEN` |

Both caches record the SHA-256 of the plugin file and are rebuilt when it
does not match. The method itself reads no table; `ENEG_ALLEN` is read by the
figure script only, for the abscissa of (d) and the count of 462 in the caption of (f).

Checks made when the figure was finished (outputs kept outside the
repository): the totals of (a) recomputed with `awk` from the five tables
equal the plotted counts and the team brief (758, 689, 589, 497, 769 of 769;
36, 27, 27, 19, 38 and 38, 38, 31, 25, 38 of 38); the rank correlations of (d)
recomputed with `scipy.stats.spearmanr` (0.9643, 0.6190); the counts of (f)
summed directly from `pair_margins(case)` over the 38 cases (566 WFs, 527 by
force, 39 by *χ*, 26 element pairs); the scan of (e) repeated with a separate
script on another grid (37 shared grid values, no difference).

## What the lower block shows and why

- (d) shows the eight elements for which the method computes *χ* on clean38,
  not all 34: *χ* exists only for an atom that carries lone pairs at a
  contested WF. The seven p-block elements follow the Allen order with one
  exchange (I, Se); Mo (MoCl2, Mo<sub>6</sub> clusters), whose lone pairs
  are own WFs displaced from the nucleus by at least 1/8 of their rms radius,
  lies above Cl. Both rank correlations are printed in the panel.
- The force rule decides 527 of the 566 contested heteronuclear WFs, so the
  agreement of the polarity with the Allen scale for those bonds is shown in
  (f): rows by bond class with the winner first, and 462 of 466 one-owner
  pairs going to the element with the larger Allen electronegativity.
- (e) scans all four constants, including the numerical tolerance, each over
  a range several times wider than the interval in which the result does not
  change.
- (f) is split by bond class and by rule, with the counts of both rules in
  the columns. The like-atom block is separate because there the decision is
  between the atoms of a homonuclear bond or cluster and a neighbour of
  another element.

## Drawing conventions

- (b): marker sizes decrease in the order square, circle, triangle so that
  coincident markers stay visible; the counts in the columns are exact.
- (c): the tile fill is the fraction correct (Cl 96.3%, Cu 84.6%, Mo 75.0%).
  The lanthanide row holds Ce.
- (d): logarithmic ordinate; bars are the interquartile range over atoms (the
  full ranges are in the table above).
- (e): the ordinate starts at 690 and ends at 769, the number of labelled
  atoms; the fourth abscissa is logarithmic.
- (f): logarithmic abscissa from 1; ratios whose loser has
  *Z*<sub>eff</sub> ≤ 0 are infinite and are drawn in the grey slot.
- Typography: page 17.1 cm wide, text 6-9 pt (Arial), formula labels through
  `style.formula`; panel labels on the grid of `fig2_benchmark.py`.
