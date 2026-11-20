# LEOS results

This document collects the test and application results of LEOS. The method and the
design choices are in [PLAN.md](PLAN.md), the terms in [CONTEXT.md](../CONTEXT.md).
Every number names its data file; where a result can be reproduced locally, the
command is given.

## Overview

| Item | Result | Data |
|---|---|---|
| 74-case benchmark (72 with reference oxidation states) | flow 64/72 (89%), greedy 51/72 (71%) | results/bench74_golden.tsv, results/bench74_greedy.tsv |
| The 40 solids of Comas-Vilà 2025 | all 40 receive an assignment | runs with role=comas in cases/INDEX.tsv, results/comas40_summary.tsv |
| Li2O2 delithiation (Li16O16 / Li15O16 / Li14O16) | mean O oxidation state −1.00 / −0.94 / −0.88, in agreement with the charge count | cases/li16o16/, cases/li15o16/, cases/li14o16/, results/li2o2_summary.tsv |
| CeO2 oxygen vacancy (Ce4O7) | the 4 Ce are [3, 3, 4, 4] | cases/ce4o7/, cases/ce4o8/ |
| Li-rich cathodes, 3 systems × 3 delithiation points × 20 frames | in the ordered systems the mean transition-metal oxidation state equals 4 + 2x, with no drift within 50 ps | results/cathode_md/, cases/{li2mno3,li2ruo3,nmc}/md-x* |
| C–H activation of methane on Ce18O36Co4 | one H transfers to the Co4 cluster as H⁻; total Co4 oxidation state 0 → +2 → +2 | cases/ce18o36co4ch4_{is,ts,fs}/v3 |

The cause of each of the 8 failed benchmark cases is known; see the "Failed cases"
table under "74-case benchmark".

## 1. Test design

### Test systems

| Benchmark set | Count | Content | Data (bench_set column of cases/INDEX.tsv) |
|---|---|---|---|
| Legacy regression set | 36 | The project's original data: 10 molecular crystals, 26 oxides and halides, including magnetic, metallic, peroxide and charged systems | cases/*/v1 (v1_mol, v1_oxide) |
| Extension wave 1 | 13 | Mixed valence (Fe3O4, CrO5), the carbonyl Cr(CO)6, the electride Ca3AuN and others | cases/*/v3 (wave1) |
| Extension wave 2 | 25 | Hydrides, Zintl phases, metal clusters, perovskites, the 160-atom SrBe(OH)4 and others | cases/*/v3 (wave2) |

The AuClB complex and CuCl4 have no reference oxidation states; they are only checked
for running to completion and are not scored. This leaves 72 scored cases.

The source check of all structures is in [testset_findings.md](testset_findings.md).
In 3 more scored cases of the benchmark the structure is not the real compound the
case name refers to: CO (actually ethylenedione), SnO (a hypothetical rock-salt
phase), Ti2O3 (a defect model: rutile with one O removed). Their composition matches
the name and the reference also holds for the actual cell content, so all count as
pass.

### Reference oxidation states

The ideal oxidation states of each system are annotated by hand from chemical
knowledge and the IUPAC rules, and stored in pyleos/sets/reference_os.json.
Mixed-valence systems carry the element mean, e.g. Fe +2.667 in Fe3O4 (1 Fe²⁺ plus
2 Fe³⁺). Where the answer itself is disputed (e.g. Ca3AuN), the annotation records
the differing views.

### Pass criterion

For each system, compute the mean oxidation state of each element, compare it with
the reference value, and take the largest absolute deviation over the elements. This
is dev.

| Status | Condition |
|---|---|
| OK (pass) | dev ≤ 0.25 |
| ~ | 0.25 < dev ≤ 1.0 |
| FAIL | dev > 1.0 |
| run | no reference; only checked for running to completion |

dev checks element means only, not single atoms. In As4S4, for example, half of the
As are +1 and half +3; the mean is +2 and dev is 0. Per-atom behaviour is described
case by case.

### The two assignment methods

- **Greedy**: the method built into the VASP patch. Each WF goes independently to the
  nearest atom (electronegativity taken into account). The result is the OS(int)
  column of the per-atom table in EOSCAR.
- **Flow**: pyleos takes the same WF data and assigns all electrons together as a
  minimum-cost flow. The VASP patch writes the same result to the OS(flow) column of
  EOSCAR.

The two methods have identical input; they differ only in how they assign.

## 2. 74-case benchmark

### Pass counts

| Benchmark set (scored cases) | Flow pass | Greedy pass |
|---|---|---|
| Legacy regression set (35) | 33 | 30 |
| Extension wave 1 (12) | 10 | 9 |
| Extension wave 2 (25) | 21 | 12 |
| Total (72) | 64 (89%) | 51 (71%) |

Case by case: 15 cases pass only with the flow, 2 pass only with greedy (the v1 data
of CSe2, and hi3); for the rest the two methods agree. Within wave 2 these two
numbers are 10 and 1.

Data: the flow in results/bench74_golden.tsv, greedy in results/bench74_greedy.tsv.
Both tables come from the same script with the same definition of dev.

**Structure-corrected version.** 12 scored cases of the 74 have a wrong structure, a
broken geometry, or are plainly not the intended phase: CO (actually ethylenedione),
SnO (hypothetical rock-salt phase), Ti2O3 (rutile with oxygen removed), CuO
(hypothetical rock-salt phase), Ce2O3 (hand-built ideal cell), ReCl4, Na2CrO4 and
NaMnO4 (isolated clusters), N2O (an ordered model 14% denser than the experimental
phase), S4N4 (a PBE cell without dispersion correction, 25% too large in volume),
li2o2 and re2 (the POSCAR scaling factor enlarged the Cartesian coordinates). They
are replaced by COD or Materials Project structures; N2O uses a P2₁3 ordered
approximant built on the experimental Pa-3 cell. After recomputation
(`cases/<system>/v31`; li2o2 uses the Föppl crystal of the homonuclear-bond test set,
cases/li2o2/v3; replacement list pyleos/sets/corrected.json) all 12 cases still pass
and the flow stays at 64/72 (results/bench74_corrected.tsv,
`python3 -m pyleos bench full --corrected --out ...`). The golden table still holds
the original runs, and CI compares against the golden table.

Four are not replaced: the correct crystal of RhCl(CO)(PPh3)2 is only in the CSD,
which needs a licence, and the present cell is an SO2 adduct; the cell of the AuClB
complex is a structure corrupted in transcription, and no real compound with this
formula exists; HI3 has no solid phase; neutral CuCl4 does not exist. The last three
also have no reference, or a reference that does not hold for the actual cell.

Reproduce (about 25 s):

```bash
python3 -m pyleos bench full --no-reliability --out /tmp/bench.tsv
python3 -m pyleos bench full --greedy --no-reliability --out /tmp/greedy.tsv
```

Without `--no-reliability` the reliability index R is computed as well, which takes
much longer. `--out` is required; write results/bench74_golden.tsv only when
regenerating the golden table on purpose.

.github/workflows/ci.yml runs the unit tests and this benchmark on every push and
compares the status and dev columns line by line with results/bench74_golden.tsv.
GitHub Actions is not yet enabled for the repository (issue #3).

### Failed cases

| Case | Status / dev | Flow result | Reference | Cause |
|---|---|---|---|---|
| As4S4 (v1 data) | ~ / 1.00 | As +1.00, S −1.00 | As +2, S −2 | 15 of the 88 WFs are unlocalised. The same system recomputed with v3 passes (see "Case studies") |
| CSe2 (two data sets: v1 and wave 1) | FAIL / 4.00 | C 0, Se 0 | C −4, Se +2 | On the Allen scale C (2.544) and Se (2.424) differ by only 0.12, so the assignment is near-degenerate (R_elem about 25) |
| Ca3AuN | FAIL / 2.00 | Au −1, one Ca at 0 | Au −3 | The result corresponds to the electride picture, which the IUPAC technical report considers more accurate (see "Case studies") |
| hi3 | FAIL / 2.00 | H −1 | H +1 | The cell actually holds HI + I2 molecular chains; H (2.300) and I (2.359) are close in electronegativity |
| RhCl(CO)(PPh3)2 | FAIL / 4.00 | Rh +5 | Rh +1 | The structure has an extra SO2 ligand, while the reference is written for the compound without SO2. Projection-assisted arbitration assigns it correctly |
| WAu12 | FAIL / 6.00 | W +6, Au mean −0.5 | W 0 | A difference between two conventions: the electronegativity rule gives W +6, IUPAC assigns W(0) by the 18-electron rule |
| WCl5 | FAIL / 4.67 | W mean +0.33, Cl mean −0.07 | W +5, Cl −1 | Wannierisation failed: 316 of the 598 WFs are unlocalised (see "Case studies") |

Inspect a single case with `python3 pyleos/cli.py case <dir> -v`, e.g.
`cases/hi3/v3`.

### Cases corrected by the flow (wave 2)

| Case | Greedy dev | Flow dev | Note |
|---|---|---|---|
| SrBe(OH)4 | 6.25 | 0.00 | Greedy assigns −1 to 50 of the 64 H, i.e. treats O–H bonds as hydrides, and pushes Sr to between +4 and +10. This is the kind of failure of the "nearest atom" criterion that Gimferrer 2020 points out |
| KMgF3 | 4.00 | 0.00 | Perovskite |
| Rb2CuCl4 | 8.00 | 0.00 | Halide |
| RhH | 5.00 | 0.00 | Metal hydride |
| RbPb | 5.88 | 0.06 | Zintl phase; 48 WFs are stray, yet the flow result stays close to the reference |
| H3PO3 | 1.00 | 0.00 | Two conventions, see "Case studies" |

## 3. Case studies

### As4S4 (realgar)

On the v3 data the flow gives As +2 and S −2 (element means), in agreement with the
reference; greedy gives As +4 and a mean S of −4.

| Data | Flow As / S | Greedy As / S | Unlocalised WFs | Status |
|---|---|---|---|---|
| v3: cases/as4s4/v3 | +2.00 / −2.00 | +4.00 / −4.00 | 0 | OK |
| v1: cases/as4s4/v1 | +1.00 / −1.00 | +4.12 / −4.12 | 15 | ~ |

Per atom: the cell has 8 As–As bonds, and every As sits on exactly one of them. The
electron pair of such a bond can only go to one end as a whole, so 8 As are +1 and
8 are +3. The element mean is unaffected; R_site = 1.8 reflects exactly this
equivalence of the two ends. `--fragment` merges the two bonded As into one unit;
each unit totals +4, i.e. +2 per As:

```bash
python3 pyleos/cli.py case cases/as4s4/v3 -v \
  --fragment 17,19 --fragment 18,20 --fragment 21,23 --fragment 22,24 \
  --fragment 25,27 --fragment 26,28 --fragment 29,31 --fragment 30,32
```

The problem of the v1 data: 15 of the 88 WFs are not localised. Their spread reaches
22.6 Å² (normal is below 1 Å²) and their centres lie in the voids between atoms. With
these 15 excluded (`--drop-stray`), the other 73 WFs are all assigned correctly,
including the electrons of the 27 As–S bonds, which all go to S. Excluding WFs lowers
the electron total, so the element means can no longer be compared directly with the
reference; this check is only possible WF by WF. On the same data greedy gives some S
14 electrons (oxidation state −8), more than the 8 that the S valence shell can hold.

There are two sets of v3 data: cases/as4s4/v3 and cases/as4s4/v31-wp. The two
calculations have the same POSCAR and INCAR and the same assignment atom by atom. The
latter is produced by the v3.1 patch, and the OS(flow) column of its EOSCAR agrees
with pyleos atom by atom.

### H3PO3: two conventions

Textbooks give P +3 in H3PO3, from the habit that "H is always +1". Under the IUPAC
electronegativity rule, H (2.300) is more electronegative than P (2.253) on the Allen
scale, so the electrons of a P–H bond belong to H. The flow follows the latter rule
and gives P +5 and −1 for the H bonded to P. The reference is annotated by the same
rule, so dev is 0.00; greedy has dev 1.00. R is about 50, which reflects how close H
and P are in electronegativity.

### Ca3AuN: electride

The IUPAC technical report on oxidation states (Karen et al.) discusses this system:
the usual treatment of other perovskites gives Au −3, but the report considers Au −1
plus 2 electrons that belong to no atom (an electride) more accurate. The method of
Comas-Vilà 2025 gives Au −3 (R = 60.4).

The flow gives two Ca at +2, one Ca at 0, Au −1 and N −3, consistent with Karen's
picture: the 2 free electrons are booked on one Ca, which therefore shows 0. The WFs
themselves support this reading: 6 WFs have centres near Ca but spreads of 5–10 Å²
(well-localised WFs are below 2 Å²), which is how free electrons look in the Wannier
representation. R_site is 7.1: which Ca carries these 2 electrons is arbitrary.

The reference is annotated as Au −3, so the case counts as failed (dev 2.00).
High-pressure K in the Comas-Vilà comparison set behaves the same way.

### hi3 and WAu12: fragments

A fragment (`--fragment`) makes a specified group of atoms receive electrons as one
unit and reports the oxidation state of the unit. It is meant for delocalised anions
and metal clusters. When electron pairs with almost equal costs appear, pyleos
suggests considering a fragment.

- **hi3**: without fragments, pyleos reports two groups of delocalised electrons
  (margin 0.34). Projection-assisted arbitration returns H to +1; with fragments
  added, the two tri-iodine units are −1 each. The actual content of the cell is
  HI + I2 chains, not I3⁻ (see testset_findings.md); the chemical reading follows the
  actual content.
- **WAu12**: with Au12 as a fragment the result is W +6 and −6 for Au12 as a whole,
  consistent with the electronegativity rule. IUPAC assigns W(0) to this system by
  the 18-electron rule. As with H3PO3, this is a difference of convention.

### Tungsten chlorides

The tungsten chloride structures of the original test set transcribed only the
asymmetric unit; the cell named WCl6, for example, holds only 1 W and 1 Cl. Results
after replacing them with literature or database structures, expanded by the space
group, and recomputing:

| System | Flow result | Note |
|---|---|---|
| WCl6 | OK, W +6 / Cl −1 | No warnings |
| WCl4 | OK | |
| WCl3 | OK, W mean +3 | The W sites within the W6 cluster are equivalent, so the value of a single W is not unique; a fragment can be used |
| WCl5 | FAIL, dev 4.67 | The two calculations, on the experimental lattice and on the self-consistently relaxed lattice, have 316 and 323 unlocalised WFs. The problem lies in the wannierisation, not in the structure transcription |

Data: cases/wcl3/v3, cases/wcl4/v3, cases/wcl5/v3, cases/wcl5/v3-relaxed (relaxed
lattice), cases/wcl6/v3; corrected structures and their sources:
cases/wcl*/structure-fixed/.

On a wrong structure that holds only the asymmetric unit, the program returns the
result for the actual composition of that incomplete cell: the cell named WCl6 that
really holds one W and one Cl gives W +1.

## 4. Comas-Vilà comparison set (40 solids)

Comas-Vilà 2025 (JCTC) tests 40 solids with a method based on QTAIM atomic overlap
matrices. We assign all 40:

- 30 are the runs with role=comas in cases/INDEX.tsv, all full MLWF calculations. The
  default cost and projection-assisted arbitration give the same element means.
- 10 overlap with the 74-case benchmark (CrO3, KMgF3, NaBF4, SrBe(OH)4, Cs3AuO,
  Rb3AuO, Ca3AuN, GaSe, WCl4, Rb2CuCl4); their results are taken from
  results/bench74_golden.tsv.

The structures come from transcription of the SI tables of the paper, experimental
CIFs and databases, 29 new structures in all; the exact structures of MgXe and the
MH4 hydrides are in `cases/<system>/structure-fixed/`. A further 6 v3-scdm-only runs
use the SCDM-only fast mode and are not part of the comparison. The system-by-system
comparison of the two methods is in
[comas_vila_2025_notes.md](comas_vila_2025_notes.md).

| Class | Systems and assignment | Note |
|---|---|---|
| Zintl phase | Ba3Si4: Ba +2, Si −1.5 | Si4⁶⁻ butterfly cluster |
| Platinide | Cs2Pt: Cs +1, Pt −2 | The IUPAC example of a metal in a negative oxidation state |
| Electride | High-pressure K: K mean 0 (−1 to +1), with stray warnings | The paper has K +1, with the electrons booked on non-nuclear attractors. LEOS has no such bookkeeping unit and books the free electrons on K, as for Ca3AuN |
| Mixed valence | YBaFe2O5, two phases: Y +3, Ba +2, Fe +2.5, O −2 | Element means correct; site charge order not checked |
| MH4 hydrides | Mg, Ca, Sr: M +2 + 2H⁻ + H2; Sc, Y: M +3 + 3H⁻ + ½H2; Zr, Ce, Th: M +4 + 4H⁻ | Metal oxidation states and H element means agree with the paper. The site distribution differs for Sc and Y: the paper has 2 H⁻ plus 2 H at −½. Greedy gives metal +10 / +9 / +10 for Sr, Ce, Th |
| Xenon fluorides | XeF2 +2, XeF3 +3, XeF4 +4 | |
| High-pressure xenon oxides | XeO, two phases: +2; XeO2, two phases: +4; XeO3, 130 GPa phase: Xe mean +5 (+4 to +6); XeO3, 150 GPa phase: +6 | The two XeO3 phases differ from the paper. There the 130 GPa phase is Xe +4 (+5 with O2 as a fragment, but with a lower R), and the 150 GPa phase is a mixture of +6 and +4 whose assignment the paper considers intrinsically ambiguous |
| Mercury(IV) | HgF4: hgf4_a and hgf4_c are +4; hgf4_b has Hg mean +3 (+2 to +4) | All three structures are +4 in the paper. How the hgf4_b structure corresponds to the paper's is not verified |
| Polyatomic anions | CsF3: F −1/3; CsF5: F −0.2; CsH3: H −1/3; LiN5: N −0.2 | The charge is shared among several atoms; fragments can be used. The H in the CsH3 cell form two isolated equilateral H3 triangles (H–H 1.24 Å), whereas the paper has linear H3; the LiN5 cell holds linear N3 plus an elongated N2, whereas the paper has a pentazole ring (see testset_findings.md). These two are not compared on the same structure |
| Noble-gas anion | MgXe (200 GPa): Mg +2, Xe −2 | Agrees with the Xe anion reported in the literature; the lattice is inferred from an equation of state, with an error of about ±1.5% |

12 of the 30 calculations of the Comas comparison set carry stray warnings (number of
strays in brackets): ba3si4 (36), xeo3_b (9), ceh4 (8), thh4 (8), k_hp (6), mgxe (5),
ybafe2o5_b (4), csh3 (3), mgh4 (3), sch4 (3), lin5 (1), xeo3_a (1).

Data: the runs with role=comas. The status column of results/comas40_summary.tsv is
the batch-run status (whether EOSCAR was produced normally); it does not say whether
the assignment agrees with the paper.

## 5. Applications

### Li2O2 delithiation: the oxidation state of oxygen

When a lithium-air battery is charged, the Li atoms of Li2O2 are removed one by one.
Each removed Li leaves a hole, which turns one peroxide (O2²⁻, −1 per O) into a
superoxide (O2⁻, the two O totalling −1). The element that changes oxidation state
here is O, whereas the DFT+U+V occupation-matrix method can only follow d/f metal
sites.

Calculations: 32-atom supercells, namely the complete Li16O16, Li15O16 with 1 Li
removed, and Li14O16 with 2 Li removed.

| System | Mean O from the charge count | Flow | Greedy |
|---|---|---|---|
| Li16O16 | −1.000 | −1.00 | −1.00 |
| Li15O16 | −0.9375 | −0.94 | −0.94 |
| Li14O16 | −0.875 | −0.88 | −0.88 |

Removing 1 Li removes 1 electron, shared among 16 O, so the mean oxidation state goes
from −1 to −15/16 = −0.9375.

Both methods agree with the charge count in the element mean; they differ in the
single atoms and the single O2 units. The 8 O2 units of the flow have exactly these
charges: 8 at −2 in Li16O16; 7 at −2 plus 1 at −1 in Li15O16; 6 at −2 plus 2 at −1 in
Li14O16. Every O lies between −2 and 0. Greedy assigns −3 or +1 to individual O (in
Li16O16 the O range from −3 to +1).

The calculations use a 2×2×1 supercell of Föppl-type Li2O2 (COD 1514084), lattice
a = 6.366, c = 7.726 Å; single points, PBE, ISPIN = 2. The earlier v3 runs had the
two columns of O coordinates swapped, which put 8 Li at the midpoints of O–O bonds
(Li–O 0.78 Å); on that geometry greedy gave −0.25 and 0.00. The table above is the
recomputation on the corrected geometry (v31); the element means of the flow are the
same as before. See "Geometry re-check" in docs/testset_findings.md.

Molecular dynamics: Li15O16 at room temperature, 40 steps, analysed every 2 steps
(20 frames in all). The flow gives a mean O of −0.94 in all 20 frames. This
trajectory starts from the relaxed, correct geometry and is not affected by the
problem above.

Data: the table above is cases/li16o16/v31, cases/li15o16/v31, cases/li14o16/v31
(SAI: LEOS/recompute_20261021); the earlier runs with the wrong geometry are kept in
the respective v3 directories. The MD is cases/li15o16/v3-md2 (20 frames);
cases/li15o16/v3-md is the first MD attempt, which crashed after the first step and
has only 1 frame.

**Position of the hole.** Physically there is only 1 hole, but in the pure PBE
calculation several O2 units appear oxidised in turn from frame to frame, and no O–O
bond shortens (the bond of a localised superoxide should shorten from 1.55 Å to about
1.33 Å). PBE thus spreads the hole over all O2 units, a known symptom of the PBE
self-interaction error.

To obtain a localised hole, Li15O16 is started from a structure with "the target
peroxide pre-shortened to 1.354 Å + asymmetric magnetic moments" and scanned with
PBE+U (U on O 2p); PBE and HSE06 single points are then run on the U = 8 eV geometry:

| Calculation | Target O2 unit (O16–O20) | Conclusion |
|---|---|---|
| U = 0 (pure PBE) | Bond length returns to 1.530 Å, −1.98; the hole is spread over 6 units (−1.7 to −2.0) | Symmetry breaking in the initial structure alone is not enough; U is needed |
| U = 4 | The hole does not land on the target unit; two O2 pairs dissociate | Localisation incomplete |
| U = 6 | The hole lands on O25–O29 (1.377 Å, −1.13); 5 O2 pairs are stretched | Localised, but with large side effects |
| U = 8 | 1.401 Å, −1.11; the other peroxides −1.93 to −2.00 | Best localisation |
| PBE on the U = 8 geometry | −1.22, plus a diffuse part of about 0.3 electrons, ghost = 7 | Partly localised |
| HSE06 on the U = 8 geometry | −1.11, the same as PBE+U | Confirmed by the hybrid functional |

The unit charges in the table are counted from the greedy occupations. A side effect
of U on O 2p is that the non-target peroxides are stretched (at U = 8, 3 pairs to
about 2.2 Å), so the conclusion needs a cross-check by a single point without U or
with a hybrid functional: the last two rows of the table.

Data: cases/li15o16/polaron-{u0,u4,u6,u8,sp-pbe,sp-hse}.

The flow cannot say reliably which O2 unit holds the hole. In all calculations the
mean O of both the flow and greedy equals −0.9375, but moving the hole to a
neighbouring unit costs only about 0.005 cost units more, and R_site is only 0.3–0.5.
Between such equivalent sites the unit on which the flow solution lands has no
physical meaning and may contradict the geometry; greedy keeps the original
assignment of each WF and so does agree with the geometry. The sound procedure:
confirm the total with R_elem (about 91), detect site degeneracy with R_site, and
leave the specific site to the bond lengths (1.40 against 1.55 Å) and the spin
density.

### CeO2 oxygen vacancy

Removing one O atom from CeO2 leaves 2 electrons, which are usually taken to sit on
two Ce (Ce⁴⁺ → Ce³⁺). The calculations use DFT+U to localise the Ce 4f electrons.

| System | Ideal result | Flow | Quality of the calculation |
|---|---|---|---|
| Ce4O8 (complete) | 4 Ce⁴⁺ | [4, 4, 4, 4] | Electron count conserved to three decimal places, no warnings |
| Ce4O7 (1 oxygen vacancy) | 2 Ce³⁺ + 2 Ce⁴⁺ | [3, 3, 4, 4] | Same |

In the legacy data, the Ce2O3 calculation without U loses 0.30 electrons (dNe
warning).

Data: cases/ce4o7/v3, cases/ce4o8/v3.

### Li-rich cathodes: oxidation-state time series

Each of the three systems is taken at three delithiation points, x = 0, 0.25 and
0.50. Workflow: 50 ps of molecular dynamics for each with MACE-MPA-0 (pretrained
model, not fine-tuned) → one frame every 2.5 ps (20 frames per point, 180 frames in
all) → for each frame a DFT single point (PBE+U, U on the transition metals) and full
MLWF → flow assignment.

| System | Result |
|---|---|
| Li2MnO3 (ordered, 3d) | At x = 0 / 0.25 / 0.50 the mean Mn oxidation state is +4.0 / +4.5 / +5.0, equal to the charge count 4 + 2x; the largest drift in each 50 ps series is 0.000 |
| Li2RuO3 (ordered, 4d) | Same as Li2MnO3: +4.0 / +4.5 / +5.0, largest drift 0.000 |
| NMC (Li1.2Ni0.13Mn0.54Co0.13O2, SQS disorder) | Li and O are +1 / −2 throughout, with no drift; the transition metals show charge fluctuations between frames; the Ni fluctuation grows from ±0.9 at x = 0 to ±1.7 at x = 0.50, and x = 0.25 does not lie between the two |

On the same data the greedy Mn mean is 3.984, with a frame-to-frame fluctuation of
±0.047.

Li2RuO3 is as stable as Li2MnO3 in the full MLWF workflow, so the problem of 4d
systems in the SCDM-only fast mode (see "Efficiency") is specific to that mode.

Under PBE+U (U on the metals only) the holes created by delithiation are all carried
by the transition metals (Mn and Ru go from +4 to +5), and O stays at −2. This agrees
with the Li2O2 result: an O 2p hole localises only with U on O or with a hybrid
functional. Which element is oxidised depends on the functional; LEOS reports the
result for the functional used.

Other notes: the sampled frames typically have 48–56 ghosts, which come from the
margin of the automatic frozen window and the metallic character of the delithiated
states, and are handled by the cap on the electron total. No O–O dimerisation occurs
in the MACE trajectories within 50 ps. The SQS initial structure of NMC relaxes
during roughly the first 1 ps (the instantaneous temperature reaches 930 K and then
falls back), which does not affect the sampled frames.

Data: results/cathode_md/timeseries_*.tsv (frame-by-frame series),
results/cathode_md/acceptance_*.json (summaries). The per-frame assignment is done by
tools/timeseries.py, which enables projection-assisted arbitration by default
(`--tier2 arbiter`).

### Fe3O4: no charge order appears

The cubic conventional cell is relaxed with the correct magnetic order, Verwey-type
alternating initial values and U(Fe d) = 4 eV. The element means are correct
(Fe +2.667, O −2.000), but no site charge order appears, and there is a ghost = 48
warning. OUTCAR shows that the system converges to a half-metallic state (bands cross
the Fermi level; magnetic moment 39.7 μB/cell, ideal value 32). Under these
conditions the cubic Fd-3m cell is a charge-disordered half-metal anyway; Verwey
charge order needs the low-temperature monoclinic supercell (P2/c or Cc) and a
stronger treatment of correlation, which are not computed here.

Data: cases/fe3o4/v3-verwey/.

## 6. Recomputation of systems the old program failed on

These 7 directories were originally computed with the v1/v2 programs and hand-written
wannier90 input, and all had problems. Results of the v3 recomputation:

| System | Symptom in the old version | v3 result |
|---|---|---|
| ReO3 | Aborted with NELECT mismatch (a metal with an odd electron count) | Re +6, O −2, R = 84; 5 ghosts excluded automatically |
| MUMNOX (K complex) | Lost 5.7 electrons; skipped | dNe = 8e-11, no warnings; K, O, N, H all correct, C in its normal mixed valence |
| LASWEH (Ag complex) | Lost 3.4 electrons; when the calculation was forced, Ag was assigned +2 or +3 | dNe = 5e-11, no warnings; Ag is +1, O, N, H all correct. S is +6 whereas the ideal value for a sulfonamide is +4; the difference is one S–C bond with almost tied costs |
| VIRCEF (Ce complex) | N and C scrambled | N and H correct; Ce is +3 (the surplus WFs vanish once the WF number is set per spin); the residual errors in C and H come from strays |
| Ce18O36Co4CH4 initial state | Disentanglement did not converge over many rounds (the Co4 cluster is metallic) | 18 Ce +4, 36 O −2, CH4 correct, Co4 fragment total 0; stray 2, ghost 6 |
| Ce18O36Co4CH4 transition state, final state | Same, and out of memory on rerun | Completed on 8 GPUs; Ce and O unchanged; the C–H bond breaks, one H transfers to the cluster as H⁻, and the Co4 total becomes +2 |

The two causes of the old failures:

1. NELECT mismatch came from the hand-written wannier90.win lacking the
   disentanglement frozen window `dis_froz_max`; it had nothing to do with the number
   of disentanglement iterations. In that batch only VIRCEF, whose input had the
   window, lost no electrons. The automatic frozen window (AUTO_WINDOW) of v3 solves
   this.
2. `kmesh_tol = 1e-5` in the old inputs is what let their 4×4×2 k-mesh pass the
   wannier90 b-vector check (usage in patch/README.md).

Data:
cases/{ce18o36co4ch4_is,ce18o36co4ch4_ts,ce18o36co4ch4_fs,lasweh,mumnox,vircef}/v3;
per-atom results in results/legacy_rerun_pyleos.txt.

### The three images of methane activation

The Co4 fragment (`--fragment 55,56,57,58`) has 60 → 58 → 58 electrons and a total
oxidation state of 0 → +2 → +2 (initial state → transition state → final state). One
H⁻ appears in each of the transition and final states, and the charge of the cell is
conserved: CH3(−1) + H⁻(−1) + Co4(+2) = 0, matching CH4(0) + Co4(0) of the initial
state.

The oxidation states only say that "one of the H is −1", not which atom. A check
against the POSCAR bond lengths: only H60 breaks its bond along the path (C–H
distance 1.10 → 2.35 → 2.92 Å) and bonds to the cluster (Co–H 1.60 / 1.66 Å), so the
geometry and the oxidation states point to the same atom.

On the same data greedy gives physically impossible values: C −6 in the final state,
H between −4 and +1 in the transition state, and a Ce mean that varies over 3.83–3.89
from image to image.

The assignment among the Co atoms is degenerate: R is 80–92 for Ce and O, but only
4–22 for the four Co. For the same total of +2, the transition-state solution books
it all on Co57, and the final-state solution gives +1 each to Co57 and Co58; both are
optimal. In this situation read the total of the Co4 fragment, not the single Co.

Ground state: each of the three images is recomputed from two initial magnetic
moments, 6 and 10. The initial state converges to moment 6, the transition state to
4, the final state to 6; the two initial values give the same energy for the final
state (−490.23423 eV). Data: results/methane_groundstate_verdict.tsv.

The md5 sums of the raw cluster data and of the EOSCAR files in the local repository
agree; the jobs are 602731, 604671 and 604672.

## 7. Parameter choices

### Electronegativity weight β

β is the weight of the electronegativity term in the cost function; the default
is 10. Full scan over β (results/scan_beta.tsv, tools/beta_scan.py):

| β | Pass count | Change relative to β = 10 |
|---|---|---|
| 3 | 55 | 7 X–H and sulfur-chain cases go from pass to fail |
| 5 | 58 | |
| 8 | 63 | RhCl(CO)(PPh3)2 becomes pass |
| 9–12 | 62 | None |
| 15 and above | | Cr(CO)6 becomes fail |
| 20 and above | | TMAO also becomes fail |

Results are stable between 8 and 12; the default is 10. When β is too small, X–H
bonds are misassigned, and without any warning. This scan predates the corrected
tungsten chloride structures, which is why β = 10 gives 62 here and not the present
64.

### Projection-assisted arbitration (TAU, PMARGIN)

The default cost uses only the geometric position of the WF centre and the
electronegativity. Projection-assisted arbitration also uses the #WP block of EOSCAR:
the weight of each WF on the PAW projectors of each atom. It is not the default and
must be switched on explicitly.

Rule: a WF is reassigned by projection only when the cost difference between its two
best candidate atoms is below TAU and the projection weight leads by more than
PMARGIN. Assignments that are geometrically clear are never changed.

This rule is chosen over "subtracting the projection weight linearly from the cost".
On the 38 cases with #WP data both reach the same best score, but the linear form
works only near one parameter value (γ = 0.5 and 3.0 both fall back to the baseline).
The projection itself also favours heavy atoms, so the linear form pulls electrons
towards heavy atoms (this happened on CSe2).

Parameter scan: a 6 × 5 grid of TAU × PMARGIN, run on the 38 #WP benchmark cases and
the 6 old-program systems (tools/arbiter_scan.py,
results/scan_arbiter_20260712.tsv).

- Switched off, the baseline is 31/37 (the denominator excludes WCl5, whose
  wannierisation failed).
- The 8 grid points with PMARGIN in 0.20–0.30 and TAU in 0.25–1.5 are all safe: they
  add no new errors.
- When PMARGIN is too low, artefacts of the H 1s projection misassign C–H and N–H
  bonds: at 0.05, 37 H in the old-program systems are misassigned (all 24 H of VIRCEF
  become −1, and 13 in MUMNOX). For TAU ≥ 0.5, every PMARGIN ≤ 0.10 is unsafe.
- The default is (TAU, PMARGIN) = (1.0, 0.20), and all neighbouring grid points are
  safe. The baseline rises from 31/37 to 33/37 (hi3 and RhCl(CO)(PPh3)2 are
  corrected) with no new errors.

Scope: on organic ligands with many C–H bonds (e.g. VIRCEF) the H 1s projection has
too much weight near the bond midpoint and may wrongly give the electrons of a C–H
bond to H. The hydrogen-containing systems among the 38 benchmark cases show no such
error. Projection-assisted arbitration therefore suits inorganic systems and
near-degenerate situations; use it with caution on molecular systems with many C–H
bonds.

## 8. Quality warnings and reliability

### Behaviour of the reliability index R

R is defined in PLAN.md, "Reliability index R". Over the 74 cases
(results/bench74_golden.tsv), the median R_elem is 83.0 for the 64 passing cases and
25.55 for the 8 failing ones. R_elem is very low for hi3 (25.2), RhCl(CO)(PPh3)2
(21.7), CSe2 (25.0 and 25.9) and WCl5 (3.6). R_site = 7.1 for Ca3AuN reflects the
arbitrary placement of the free electrons; R = 49.7 for H3PO3 reflects how close H
and P are in electronegativity.

A high R does not mean agreement with the reference. Ca3AuN (75.6) and WAu12 (74.0)
fail with a fairly high R_elem; in both, the convention of the reference differs from
the electronegativity rule. R measures how stable the solution is for the given
input. If the input itself is systematically wrong (wrong structure, failed
wannierisation), the solution can be both stable and wrong. To judge whether a result
can be trusted, read R together with the three warnings stray, ghost and dNe.

### Re-check of magnetic systems with ghost warnings

13 spin-polarised systems with ghost warnings were recomputed with the v3.1 patch
(WF number set per spin); the element means do not change. In the spin-asymmetric
systems the ghost count drops, e.g. from 42 to 20 for Fe3O4.

Data: results/ghost_ab_summary.tsv.

### Systematic bias of ghosts on metal oxidation states

When the cap on the electron total removes surplus WFs, it removes the costliest
first. The cost contains the electronegativity term, so electrons lost upstream are
booked mostly against the least electronegative atoms (usually metals): Ag in LASWEH
was once assigned +2 or +3 (it should be +1), and the same happened to Ce in VIRCEF
and C in MUMNOX.

So with ghost > 0 the oxidation state of an electropositive metal may be
overestimated by 1–2; with ghost = 0 the problem does not arise. The remedy is to fix
the upstream calculation (add a frozen window, or recompute with the automatic frozen
window).

A VIRCEF comparison: offline wannier90 with a shared WF number gives Ce +4, library
mode with the WF number set per spin gives Ce +3; the difference comes from this
bias. Data: cases/li2mno3/w90-library, cases/li2mno3/w90-offline,
cases/vircef/w90-offline.

## 9. Efficiency

### SCDM-only fast mode

SCDM-only skips the iterative optimisation of the Wannier functions and assigns
directly from the centres of the SCDM initial guess.

| Test | Agreement with the full calculation | wannier90 time |
|---|---|---|
| 38 cases | element means 32/38 (84%), per-atom 19/38 | about 86 times lower (CO2: 2 s → 0.023 s) |
| 96-atom Li2MnO3 delithiation frames, 12 pairs | element means 0/12 | 67 times lower (554 → 8.3 s/frame) |

The 6 cases among the 38 that disagree fall into two classes: halides of 4d/5d
transition metals (WCl4/5/6, MoCl2; d orbitals are sensitive to the initial guess)
and delocalised bonds (GaSe). CSe2 happens to give the reference in the fast mode.

Li2MnO3 is a mixed-valence system with metallic character, and the fast mode fails on
every frame, but every failed frame carries a warning of about 288 strays (0–4 in the
full calculation). The tiered mode of the batch script
(`batch.py stage --mode tiered`) switches to the full calculation automatically on a
quality warning or when 4d/5d metals are present; all 12 frames were redirected to
the full calculation.

All other results in this document use full MLWF. The fast mode is only for
high-throughput screening of small systems or insulators and for frame-by-frame
analysis of trajectories.

The raw data of the 38-case comparison are not in the repository.

### Offline wannier90 and threaded MKL

wannier90 is split off from the GPU job and run offline on CPU nodes. The test frame
is the 96-atom Li2MnO3, with the two spins run at the same time
(results/w90_e2e_pipeline.tsv):

| Mode | wannier90 wall time | GPU part |
|---|---|---|
| wannier90 library mode inside the GPU job | 506.5 s | whole frame 816.9 s, all 4 cards occupied throughout |
| GPU only exports, wannier90 runs offline on CPU | 70.0 s | 307.1 s (GPU occupancy reduced by 62%) |
| Same, plus threaded MKL (2 threads per process) | 62.0 s | 307.1 s |
| Same, with the WF number of the second spin reduced from 288 to 240 | 48.1 s | 307.1 s |

Threaded MKL: the reference maths library does not respond to the thread count. With
threaded MKL, at 8 processes × 4 threads, disentanglement takes 0.055–0.065 s per
iteration and localisation 0.074 s, about 16 and 13 times faster than the single-core
reference build (scan of 18 settings, results/w90_threads_scan.tsv). Channels that
converge normally give results bit-identical to the reference build, and Ω_I is the
same in all 18 runs.

Two limits. This system has few k-points (nk = 8), and the fixed IO of reading
.mmn/.amn, about 30 s, takes most of the wall time; the larger the system, the more
threading gains. With the two spins running concurrently, memory bandwidth is the
bottleneck: 2 threads per process is best, and 4 threads is slower (72.2 s).

On converged channels the offline mode and the library mode give the same result WF
by WF. Build and run instructions are in patch/README.md.

### gpuw90

gpuw90 is the PyTorch rewrite, in this repository, of the disentanglement and
localisation of wannier90 (one code for CPU and CUDA). On the 4 entangled-band
acceptance channels Ω_I differs from wannier90 on the order of 10⁻⁹ Å², and the
number of disentanglement iterations is the same as in wannier90.

On a 96-atom Li-rich cathode frame (288 WFs), a V100 takes 0.079 s per
disentanglement iteration and 0.041 s per localisation iteration, about 14 and 22
times faster than serial wannier90 on the cluster; the torch version on pure CPU is
about 4.6 times faster. Against the baseline of 554 s/frame, the Wannier part for the
two spins would come down to about half a minute.

Details in [gpuw90/README.md](../gpuw90/README.md).

### Warm start

A warm start writes the Wannier centres of the previous frame into the wannier90
input of the next frame as guiding centres. On a 40-step AIMD of the 96-atom
delithiated Li2MnO3 supercell (analysed every 5 steps, 8 analysis steps in all) it
gives no speed-up:

| Setting | Mean analysis step | Whole trajectory |
|---|---|---|
| Warm start on | 572 s | 6742.6 s |
| Warm start off | 546 s | 6532.4 s |

The reason is that the initial guess here comes from SCDM and is already good enough,
so the guiding centres only add overhead. The option is kept for cases without an
SCDM initial guess. Also, an analysis step costs about 525 s more than an ordinary MD
step (wannier90 plus LEOS); this is the upper limit of what SCDM-only can save in a
trajectory.

Data: results/cathode_md/warmstart_*.

### Analysis time

On a laptop pyleos assigns the 74 cases by the flow in about 25 s (pure Python
standard library, without R); computing R takes a few minutes.

## 10. Known limitations

- **Occupations are not updated along a trajectory.** EOS_PREPARE_MLWF is called only
  at the first step; from step 2 on, the trajectory still uses the WF number of the
  first step. If the occupation increases, the electron-conservation check
  (EOS_NELECT_ABORT) aborts the calculation and reports it; no wrong result is
  produced.
- **Length limit of the wannier90 input.** The .win is processed up to 20000 lines of
  255 columns each; anything beyond is truncated. Real inputs are far below this
  limit; if it is exceeded, wannier90 reports an error.
- **WCl5**: wannierisation fails on both lattices (see "Case studies").
- **The residual C and H errors of VIRCEF** come from strays, i.e. from the quality
  of the wannierisation, not from ghosts.
- **#WP data differ between platforms.** The per-atom aggregated #WP values from the
  GPU (nvhpc) and CPU (Intel) builds of the patch differ by about a factor of two,
  for an unknown reason. #WP data from different platforms must not be mixed (see
  patch/README.md).
- **Site degeneracy.** The assignment of electrons among equivalent atoms (the two
  ends of a homonuclear bond, the atoms of a metal cluster, equivalent O2 units) is
  not unique. Read the element mean or the fragment total; single-atom values are
  unreliable.
