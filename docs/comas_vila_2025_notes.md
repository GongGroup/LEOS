# Comas-Vilà et al. (JCTC 2025): reading notes

These notes record how Comas-Vilà et al. extend effective oxidation state (EOS) analysis
to solids, their reliability index R(%), how the method compares with LEOS, and the LEOS
results on their 40-solid test set. For LEOS numbers, docs/RESULTS.md and data/ prevail.

## Reference

| Item | Content |
|---|---|
| Authors | Gerard Comas-Vilà, Leila Pujal, Alberto Otero-de-la-Roza, Davide Tiana, Julia Contreras-García, Pedro Salvador |
| Title | Effective Oxidation State Analysis for Solids |
| Source | *J. Chem. Theory Comput.* 2025, 21, 7075–7086 |
| DOI | 10.1021/acs.jctc.5c00482 |
| PDF | docs/refs/Comas-Vilà 等 - 2025 - Effective Oxidation State Analysis for Solids.pdf |
| Molecular predecessor | Ramos-Cordoba, Postils, Salvador, *J. Chem. Theory Comput.* 2015, 11, 1501–1508; reference 20 of the paper |

Page numbers below are journal pages (7075–7086).

How the paper relates to LEOS:

- Same question: per-atom integer oxidation states in a solid from DFT wavefunctions.
- Different route. They use QTAIM atomic domains and atomic overlap matrices; LEOS
  uses Wannier centre positions and a capacity-constrained minimum-cost flow.
- The name "EOS" comes from the molecular method.
- Their 40-solid test set is the source of the LEOS extension benchmark. LEOS has
  assignments for all 40; see "Comparison on the 40 solids" below.

## Their method

They build effective atomic orbitals and their occupations in each atomic domain, then
sort all orbitals of the system by occupation and assign integer electrons in one pass.

### Idea

The Theory section (p.7076–7077) gives the steps below for a single-determinant
wavefunction and QTAIM partition; the solid version keeps the same algorithm (p.7078):

1. Partition space. QTAIM (Bader's topological partition of the electron density)
   cuts space into non-overlapping atomic domains Ω_A, one basin per atom.
2. Compute the atomic overlap matrix (AOM). For each atomic domain, compute the
   overlap integrals of the occupied orbitals inside the domain:
   S_ij^(σ,A) = ∫_{Ω_A} φ_i^σ*(r) φ_j^σ(r) dr (eq. 1, p.7077).
3. Diagonalise the AOM. The eigenvectors are effective atomic orbitals (eff-AOs), defined
   only inside Ω_A; each eigenvalue λ, between 0 and 1, is the orbital's occupation
   (eqs. 2–4, p.7077): near 1 for core shells and lone pairs, intermediate for orbitals
   involved in bonding. For fragments they are called effective fragment orbitals (EFOs).
4. Global assignment. Sort the eff-AOs of all atoms or fragments by decreasing
   occupation; the first n^σ are occupied, the rest empty; n^σ is the electron count
   of spin channel σ. An atom's oxidation state is its atomic number minus the α and
   β electrons it receives. Spin-polarised systems treat the two channels separately.

Step 4 does not round the occupations one by one; all occupations of the system compete
and the list is cut (p.7077; stressed again for Rb2CuCl4 and XeO2, p.7081, p.7082). So
in an electron-deficient system the last occupied eff-AO can have an occupation below
0.5, and in an electron-rich one an eff-AO above 0.5 can be left empty (p.7077).

### Program chain

Plane-wave pseudopotential DFT does not give the all-electron density QTAIM needs, so
they chain four programs (start of Results, p.7078; Computational Details, p.7084).

| Step | Program | Input | Output | Purpose |
|---|---|---|---|---|
| 1a. Self-consistent DFT | Quantum ESPRESSO, norm-conserving pseudopotentials, PBEsol | crystal structure | Bloch states, electron density | provides the wavefunctions |
| 1b. All-electron density reconstruction | Quantum ESPRESSO, a second, PAW calculation with the same geometry and grid | crystal structure | all-electron density on a grid | the QTAIM topological analysis needs the all-electron density |
| 2. Wannierisation | wannier90 | Bloch states from 1a | MLWFs | turns Bloch states into localised orbitals, which replace the molecular orbitals of the molecular case |
| 3. Atomic-domain integration | critic2 | all-electron density, MLWFs | one AOM per atom | QTAIM topological analysis finds the basins; the Henkelman grid method then integrates products of MLWFs in each basin |
| 4. EOS assignment | APOST-3D | all AOMs | eff-AO occupations, integer oxidation states, R(%) | diagonalises the AOMs and does the global sort of step 4 above |

Parameters (p.7084): PBEsol functional; wavefunction cutoff 80 Ry, density cutoff
320 Ry; k-point grid 4×4×4.

### Computational cost

All-electron density reconstruction and basin grid integration dominate; no timings given.

- Two DFT runs. The norm-conserving run gives the Bloch states, the PAW run the
  all-electron density; geometry and grid must be identical (p.7084).
- Basin integration. Each atomic domain needs an N_wann × N_wann AOM, each element a
  numerical 3D grid integral. The 320 Ry density cutoff means a dense real-space grid.
- Metals. Partly filled bands cannot be wannierised, so the AOM comes straight from Bloch
  states (p.7084); the authors call this less efficient and costlier (Rb3AuO, p.7079).
- Electrides. The free electrons are booked on non-nuclear attractors (NNAs) of the
  density, used as extra basins. In Ca2N the authors find the presence and number of
  NNAs very sensitive to the k-point count and the topological-analysis grid (p.7079).

### Input and output

- Input: the crystal structure. The paper's 40 structures come from COD, CCDC or the
  original literature (Table 1, p.7079).
- Output: an integer oxidation state per atom or fragment, the eff-AO occupations of
  each spin channel, and one R(%) per cell. The frontier eff-AO occupations and the
  assignments of the 40 systems are in SI Tables S1–S16 (p.7078).
- Fragments: QTAIM atomic domains do not overlap, so a fragment's AOM is the sum of
  its member atoms' AOMs (p.7078).

## Reliability index R(%)

R(%) rates how clear the assignment is by the occupation gap of the two frontier eff-AOs.

### Definition

After sorting spin channel σ, let λ_LO^σ be the occupation of the last occupied
eff-AO and λ_FU^σ that of the first unoccupied one. Eq. 5 (p.7077):

    R^σ(%) = 100 · min(1, max(0, λ_LO^σ − λ_FU^σ + 1/2))

The cell takes the smaller of the two spin channels: R(%) = min(R^α, R^β). The solid
version treats the frontier orbitals differently from the molecular one (p.7078):

- The molecular version requires the two frontier eff-AOs to be on different fragments.
- The solid version lets them belong to the same atom, by a supercell argument. The
  paper's example is a binary solid AB: under the molecular rule the primitive cell gets
  R = 90, but in a doubled supercell the frontier orbitals fall on two equivalent atoms
  and R changes. Without the restriction, primitive cell and supercell give the same R.

### Meaning of the values

| R(%) | Meaning | Source |
|---|---|---|
| 100 | λ_LO − λ_FU ≥ 0.5; the authors consider the assignment indisputable | p.7077 |
| 50–100 | the higher the value, the closer the formal ionic picture is to the actual electron distribution | p.7077 |
| 50 | the two frontier eff-AOs have equal occupations, i.e. they are degenerate | p.7078 |
| below 50 | the assignment did not go from high to low occupation. This never happens in a normal EOS run, only when a set of oxidation states is imposed by hand to test how well the system fits it | p.7077–7078 |

When symmetry causes the degeneracy, the authors share the last electrons equally among
the degenerate eff-AOs, which gives fractional oxidation states. General rule: p.7078;
examples: GaSe (p.7080), MH4 hydrides (p.7081–7082), high-temperature YBaFe2O5 (p.7081).

The authors make two remarks on the meaning of R(%) (p.7078):

- R is not a probability. R = 65% does not mean the assignment is 65% likely to be right;
  it measures how close the formal ionic picture is to the actual electron distribution.
- The functional form of eq. 5 is somewhat arbitrary. It matters more to identify the
  frontier eff-AOs and analyse their shapes and occupations.

### Uses in the paper

The paper uses R(%) in four ways, each with examples.

| Use | Example | Source |
|---|---|---|
| Flag whether an assignment is clear | CrO3: Cr is +6 while its QTAIM partial charge is only +1.76; R = 71. The authors use this to show that partial charges cannot be used to infer oxidation states | p.7078–7079 |
| Choose between fragment partitions | GaSe: R = 50 by atoms, R = 77.7 with Ga–Ga as a fragment, hence the (Ga–Ga)⁴⁺ picture. MH4: Ha–Ha fragments are clearer for Mg/Ca/Sr; the partition by atoms has the higher R for Zr/Ce/Th | p.7079–7080; p.7081–7082 |
| Follow trends | HgF4: R rises monotonically with pressure over the three structures at 50, 100 and 200 GPa, so the Hg(+4) picture becomes more credible | p.7083 |
| Show where the method does not apply | LiN5: the five N have almost equal occupations, R = 50, no atom-level assignment. XeO3 150 GPa phase: R = 51.2; the authors consider the assignment inherently ambiguous | p.7084; p.7083 |

The Conclusions section calls the last use the method's ability to "self-assess its
limits of applicability" (p.7084).

## Comparison with LEOS

Both methods get integer oxidation states from one global assignment; neither rounds
partial charges. They differ in bookkeeping unit, data needed and computational chain.

### Principle

| | Comas-Vilà | LEOS |
|---|---|---|
| Bookkeeping unit | eff-AO; its occupation is continuous between 0 and 1 | Wannier function (WF), each holding an integer number of electrons: 2 if closed-shell, 1 per channel if spin-polarised |
| Basis of attribution | all eff-AOs sorted by occupation | minimum-cost flow with cost = d² + β·(1 − χ/χ_max). d is the distance from the WF centre to the nearest periodic image of the atom (Å), χ the Allen electronegativity, β = 10 |
| Integrality | the first n^σ after sorting | supplies and capacities are integers, so the optimum is automatically integer |
| Physical constraints | no explicit capacity constraint | each atom has an upper and a lower capacity bound: the upper is the physical capacity of the valence shell, the lower the semicore electrons |
| Data needed | one AOM per atom, N_wann × N_wann, whose elements are basin integrals | WF centre coordinates and spreads (EOSCAR), crystal structure (POSCAR) |
| Space partition | QTAIM atomic domains | none; only the distances from the centres to the atoms |

They put the chemical information in the overlap integrals over the domains. LEOS
uses only the WF centre positions and the element electronegativities, and capacity
constraints exclude physically impossible assignments, such as S taking 14 electrons.
For the flow formulation see docs/PLAN.md, "Assignment as a minimum-cost flow".

### Computational cost

| | Comas-Vilà | LEOS |
|---|---|---|
| Programs | QE (two runs) → wannier90 → critic2 → APOST-3D | one run of LEOS-patched VASP, which calls wannier90 in library mode and writes EOSCAR; pyleos can recompute offline |
| Main cost | all-electron density reconstruction, QTAIM topological analysis, basin grid integration | the DFT and the wannierisation themselves; the assignment problem has only tens to hundreds of variables and solves in milliseconds |
| Metals | AOM from Bloch states, slower | still WF centres; the surplus WFs are excluded as ghosts, which is an approximation |
| Timing | not reported in the paper | pyleos recomputes the 74 cases in about 25 s on a laptop, without the R index |

A ghost is a surplus WF with no electron to hold; LEOS excludes ghosts with the cap on the
electron total (the true electron count). Known bias: the cap drops the highest-cost edges
first and the cost has an electronegativity term, so electrons lost upstream all land on
the least electronegative atoms. With ghost > 0 the oxidation state of electropositive
metals can be 1–2 too high (docs/RESULTS.md, "Systematic bias of ghosts on metal oxidation
states"). Their Bloch-state route for metals costs more but avoids this approximation.

Their cost comes mainly from the current implementation, not from the EOS framework:

- The molecular EOS suits Hilbert-space and real-space partitions alike (p.7076–7077).
- The Conclusions say the scheme extends to any real-space partition, including Voronoi
  polyhedra and Hirshfeld-type partitions, provided its AOM can be computed (p.7084).
- For electrides the authors suggest that a partition combining ELF and Voronoi
  (Weaver et al.) may suit better than QTAIM (p.7079).

A cheaper partition still needs an N_wann × N_wann matrix of domain integrals per atom;
LEOS does no spatial integration. Cost is also why docs/PLAN.md, "The default cost uses
only geometry and electronegativity", keeps the full AOM out of the LEOS cost function.

### Features: trajectories

All 40 systems of the paper are single-point calculations, with no molecular dynamics.
In their chain, every frame would need its own PAW reconstruction and basin integration.

The LEOS EOSCAR carries the in situ assignment column OS(flow), so MD can be assigned
frame by frame. Trajectory results so far:

| System | Size | Result | Source |
|---|---|---|---|
| Li15O16 | 40 steps of room-temperature MD, analysed every 2 steps, 20 frames | the element-mean oxidation state of O is −0.94 in all 20 frames | docs/RESULTS.md, "Li2O2 delithiation: the oxidation state of oxygen" |
| Li-rich cathodes Li2MnO3, Li2RuO3, NMC | 3 systems × 3 delithiation points × 50 ps, 180 frames | oxidation-state time series | docs/RESULTS.md, "Li-rich cathodes: oxidation-state time series"; results/cathode_md/ |

Per-frame wannierisation time, for a 96-atom delithiated Li2MnO3 frame: the wannier90
stage of a full MLWF run takes 554 s/frame; the SCDM-only fast mode takes 8.3 s/frame,
but on this system the element means agree in 0/12, and the stray warning catches every
failed frame; gpuw90, scaled from its measured iteration speed, needs about half a minute
for the Wannier stage of both spins. Data and scope: docs/RESULTS.md, "Efficiency".

### Features: fragments

Both sides can merge several atoms into one fragment, for cases where the picture
"each electron pair belongs to one or two atoms" does not hold.

| | Comas-Vilà | LEOS |
|---|---|---|
| Implementation | fragment AOM = sum of the member atoms' AOMs (p.7078) | the member atoms merge into one sink: capacities add, the cost is that of the nearest member. CLI option `--fragment` |
| Use cases | OH in SrBe(OH)4, Si4 in Ba3Si4, Ga–Ga in GaSe, Ha–Ha in MH4, O2 in XeO3, H3/F3/F5 in CsH3/CsF3/CsF5 | the I3 unit of hi3, Au12 of WAu12, Co4 of Ce18O36Co4CH4 |
| When to use a fragment | compare R(%) between partitions | when the two candidate atoms of an electron pair have close costs, the CLI suggests `--fragment` |
| Symmetry degeneracy | share the last electrons equally among equivalent eff-AOs, giving fractional oxidation states | conclude at the level of the element mean or the fragment total |

For the validation on the LEOS side see docs/RESULTS.md, "hi3 and WAu12: fragments".

### Features: reliability index

The LEOS counterparts of R(%) are R_site and R_elem. Definition: after solving the
minimum-cost flow, force the electron count of one atom (R_site) or of one element
total (R_elem) to change by ±1 and re-solve globally; the extra cost m is normalised
to 0–100 as R = 100·m/(m+1). The R_site of the cell is the minimum over all atoms.

| | R(%) | R_site / R_elem |
|---|---|---|
| What is measured | occupation gap between the two frontier eff-AOs | extra cost of the global re-solve after a forced ±1; chain rearrangements are fully costed |
| Normalisation | 100·min(1, max(0, λ_LO − λ_FU + 1/2)) | 100·m/(m+1) |
| Granularity | one number per cell, the smaller of the two spin channels | two levels: R_site counts any ±1 rearrangement; R_elem only those that change an element total |
| Under site degeneracy | R = 50, e.g. GaSe by atoms (p.7080) | R_site near 0, e.g. As4S4 1.8, GaSe 0.2 |
| Known limits | not a probability (p.7078) | with systematically wrong input, R can be high and the assignment wrong; read it together with stray, ghost and dNe |

A stray is an unlocalised WF: spread above 3 Å², or centre more than 1.6 Å from every
atom. dNe is the difference between the total electron count in the Wannier
representation and the DFT value. R_site isolates site degeneracy, such as a
homonuclear bond split half and half or equivalent sites swapping; R_elem looks only
at element totals, which matches the element-mean scoring of the 74-case benchmark.
Median R_elem of two groups of cases among the 74:

| Group (status column) | Cases | Median R_elem |
|---|---|---|
| Correct (OK) | 64 | 83.0 |
| Deviating (FAIL or ~) | 8 | 25.55 |

Data: results/bench74_golden.tsv, grouped by the status column; R_site of As4S4 and
GaSe is from the As4S4 and sega rows of that file. Definition and reading:
docs/RESULTS.md, "Behaviour of the reliability index R".

### Features: electrides

The two methods detect free electrons differently.

| | Comas-Vilà | LEOS |
|---|---|---|
| Basis | NNAs of the density, which take part in the assignment as extra basins | stray: the WF centre lies near an atom but the spread is abnormally large |
| High-pressure K | 2 NNAs found, each assigned 2 electrons; K is +1, R is low (p.7079) | K averages 0; the free electrons are booked on K, with a stray warning |
| Ca3AuN | no NNA found; Au(−3), R = 60.4 (p.7079) | Ca₂(+2), Ca(0), Au(−1), N(−3); 6 WFs have spreads of 5–10 Å², against normal values below 2 Å² |
| Known problem | presence and number of NNAs are sensitive to numerical parameters, see Ca2N (p.7079) | which atom the free electrons are booked on is arbitrary; R_site of Ca3AuN = 7.1 |

The two pictures of Ca3AuN: as the paper relates it, the IUPAC technical report
(Karen et al.) notes that following the practice for other perovskites gives Au(−3),
which violates the 12−N rule, so the report considers the electride picture
Au(−1) + 2e⁻ more accurate (p.7079). The EOS analysis of the paper assigns Au(−3).
LEOS gives Au(−1) and books the 2 free electrons on one Ca, which therefore shows 0.
Both reliability readings are low: R = 60.4 in the paper, R_site = 7.1 in LEOS.

Data: ca3aun row of results/bench74_golden.tsv; see docs/RESULTS.md, "Ca3AuN: electride".

## Comparison on the 40 solids

On the 40 solids most element-mean oxidation states agree with the paper; the differences
are the electrides (Ca3AuN, high-pressure K), high-pressure XeO3 and one HgF4 structure.

### Test set and basis of the references

The paper divides the 40 systems into these groups (p.7078):

| No. | Systems | What the paper compares with |
|---|---|---|
| 1–15 | examples from the IUPAC technical report (Karen et al.): metal oxides, perovskites, antiperovskites, transition-metal systems | the IUPAC iBOS or BVS algorithm. The authors report full agreement (p.7081, p.7084) |
| 16–23 | the MH4 metal hydrides of Bi and Zurek | assignments of the original literature |
| 24–33 | Xe compounds: high-pressure oxides 24–29, fluorides 30–32, high-pressure MgXe 33 | assignments of the original literature |
| 34–36 | high-pressure HgF4 | assignments of the original literature |
| 37–40 | systems with homonuclear bonds: CsH3, CsF3, CsF5, LiN5 | assignments of the original literature |

The original assignments for 16–40 rest on partial charges, ELF shapes or geometric
criteria (p.7078): for example QTAIM charges and H–H bond lengths for MH4 (p.7081),
ELF shapes for the Xe oxides (p.7082), PDOS for HgF4 (p.7083). The last 25 systems
have no IUPAC reference oxidation states; where LEOS and the paper differ on them, the
two methods simply give different results, and this cannot show which one is right.

### Location of the LEOS data

- 30 systems are in `cases/<case>/v3` (role=comas in INDEX.tsv), all full MLWF runs.
- 10 systems are also in the 74-case benchmark; their assignments come from
  results/bench74_golden.tsv: cro3, kmgf3, nabf4, srbe_oh_4, cs3auo, rb3auo, ca3aun, sega
  (GaSe), wcl4, rb2cucl4. 9 are OK; ca3aun is FAIL because the benchmark reference takes
  Au(−3) (pyleos/sets/reference_os.json). WCl4 uses the corrected structure cases/wcl4/v3.
- 6 further v3-scdm-only runs use the SCDM-only fast mode (EOS_SCDM_ONLY = .TRUE. in
  INCAR) and are not comparison results. cse2, mocl2, wcl5 and wcl6 are not in the
  paper's test set; the v3-scdm-only runs of sega and wcl4 share the names of GaSe and
  WCl4, but the cases/wcl4/v3-scdm-only cell is in fact WCl3 (docs/testset_registry.tsv).
- results/comas40_summary.tsv has 36 rows (30 + 6), all with status OK. This column
  is the batch-run status: EOSCAR was produced normally (the triage of
  pyleos/batch.py). It does not mean the assignment agrees with the reference.
- Structure sources: 29 new structures from table-by-table transcription of the SI,
  experimental CIFs and same-phase database structures, plus the exact high-pressure
  phases in `cases/<system>/structure-fixed/` (MgXe and the four hydrides of
  Mg/Sr/Ce/Th). The source and the audit verdict of each structure are recorded in
  docs/testset_registry.tsv; cs2pt, csh3 and lin5 are marked SUSPECT.
- Settings: β = 10, Allen electronegativity. The 30 runs of the comparison set (role=comas
  in INDEX.tsv) give the same element means with the default cost
  (geometry + electronegativity) and with projection-assisted arbitration (default TAU and
  PMARGIN). The 10 shared systems use the default-cost results of the 74-case benchmark.

### Comparison by class

R in the paper column always means R(%). The LEOS column is written as in
docs/RESULTS.md, "Comas-Vilà comparison set (40 solids)". Rows marked OK come from
results/bench74_golden.tsv: OK means dev = 0.00, i.e. each element-mean oxidation state
equals the benchmark reference; the values shown are from pyleos/sets/reference_os.json.

| No. | System | Paper (page) | LEOS | LEOS data |
|---|---|---|---| --- |
| 1 | CrO3 | Cr +6, R = 71 (p.7078–7079) | Cr+6 / O−2, OK | bench74_golden.tsv |
| 2–4 | KMgF3, NaBF4, SrBe(OH)4 | agree with IUPAC, R at its maximum for all (p.7079) | K+1/Mg+2/F−1, Na+1/B+3/F−1, Sr+2/Be+2/O−2/H+1, all OK | bench74_golden.tsv |
| 5–6 | Cs3AuO, Rb3AuO | agree with IUPAC, R = 100; Rb3AuO is a metal, AOM from Bloch states (p.7079) | Cs+1/Au−1/O−2, Rb+1/Au−1/O−2, both OK | bench74_golden.tsv |
| 7 | Ca3AuN | Au −3, R = 60.4, no NNA found (p.7079) | Ca₂(+2), Ca(0), Au(−1), N(−3); FAIL against the benchmark reference | bench74_golden.tsv |
| 8 | K (high pressure) | 2 NNAs with 2 electrons each, K +1, R low (p.7079) | K averages 0; free electrons booked, with stray | cases/k_hp/v3 |
| 9 | GaSe | R = 50 by atoms, Ga +2 after equal sharing; Ga–Ga fragment R = 77.7 (p.7080) | Ga+2 / Se−2, OK | sega row of bench74_golden.tsv |
| 10 | Ba3Si4 | Si4 as fragment: Ba +2, Si4 −6, R = 100 (p.7080) | Ba+2 / Si−1.5, Si4⁶⁻ butterfly cluster | cases/ba3si4/v3 |
| 11 | Cs2Pt | Pt −2, assignment undisputed (p.7080) | Cs+1 / Pt−2 | cases/cs2pt/v3 (recomputed with the experimental lattice: cases/cs2pt/v31) |
| 12 | WCl4 | W +4, R = 92.6 (p.7080) | W+4 / Cl−1, OK | bench74_golden.tsv; cases/wcl4/v3 |
| 13 | Rb2CuCl4 | Cu +2, β channel R = 75.4 (p.7080–7081) | Rb+1/Cu+2/Cl−1, OK | bench74_golden.tsv |
| 14–15 | YBaFe2O5, two phases | low-temperature phase: Fe1 +3, Fe2 +2, R = 55.6; high-temperature phase: both Fe shared equally at +2.5, R = 61.3 (p.7081) | both phases Y+3 / Ba+2 / Fe+2.5 / O−2; agrees at element level, site charge order is a separate matter | cases/ybafe2o5_a/v3, cases/ybafe2o5_b/v3 |
| 16–18 | MgH4, CaH4, SrH4 | M +2, Hb −1; Ha–Ha as a fragment is neutral H2, clearer than the partition by atoms (p.7081) | M+2 + 2H⁻ + H2; site degeneracy inside the H2 unit | cases/mgh4/v3, cases/cah4/v3, cases/srh4/v3 |
| 19–20 | ScH4, YH4 | M +3, Hb −1; the last electron is shared equally between the equivalent Ha, −½ each (p.7081–7082) | M+3 + 3H⁻ + ½H2 | cases/sch4/v3, cases/yh4/v3 |
| 21–23 | ZrH4, CeH4, ThH4 | M +4, all four H −1; R is higher by atoms, 75.4 for ZrH4 (p.7082) | M+4 + 4H⁻ | cases/zrh4/v3, cases/ceh4/v3, cases/thh4/v3 |
| 24–25 | XeO | 100 GPa: Xe +2, R = 74.6; 200 GPa: Xe1 +4, Xe2 0, R = 73.3 (p.7082) | +2 in both phases | cases/xeo_a/v3, cases/xeo_b/v3 |
| 26–27 | XeO2 | 150 GPa: Xe +4, R = 77.6; 200 GPa: Xe +4, R = 58.7 (p.7082) | +4 in both phases | cases/xeo2_a/v3, cases/xeo2_b/v3 |
| 28 | XeO3 (130 GPa) | by atoms: Xe +4, O1 −2, O2 0, R = 60.9; with O2 as a fragment, mixed-valence Xe +5 but a lower R (p.7082–7083) | +5, with stray | cases/xeo3_a/v3 |
| 29 | XeO3 (150 GPa) | Xe1 +6, Xe2 +4, R = 51.2; the authors consider it inherently ambiguous (p.7083) | +6 | cases/xeo3_b/v3 |
| 30–32 | XeF2, XeF3, XeF4 | +2; XeF3 is a mixture of XeF2 and XeF4 molecules, Xe1 +2, Xe2 +4; +4. R is 89–98 (p.7083) | XeF2 +2 / XeF3 +3 / XeF4 +4. The +3 of XeF3 is the element mean; per site it is +2 and +4 | cases/xef2/v3, cases/xef3/v3, cases/xef4/v3 |
| 33 | MgXe | Mg +2, Xe −2, R = 59.6 (p.7083) | Mg+2 / Xe−2 (200 GPa); the lattice is an equation-of-state estimate, error about ±1.5% | cases/mgxe/v3 |
| 34–36 | HgF4 | all three structures are I4/m, all Hg +4; R rises with pressure, 56.3 for No. 34 (p.7083) | hgf4_a and hgf4_c are +4; in hgf4_b Hg averages +3 (+2 to +4) | cases/hgf4_a/v3, cases/hgf4_b/v3, cases/hgf4_c/v3 |
| 37 | CsH3 | H3 as a fragment is −1; by atoms the central H is +1, R = 68 (p.7083) | H−1/3; structure in doubt, see the note after the table | cases/csh3/v3 |
| 38 | CsF3 | F3 as a fragment is −1; by atoms R = 66 (p.7083) | F−1/3, i.e. F3⁻ | cases/csf3/v3 |
| 39 | CsF5 | F5 as a fragment is −1; by atoms Cs +1, F −1, F2 0 (p.7084) | F−0.2, i.e. F5⁻ | cases/csf5/v3 |
| 40 | LiN5 | the five N have almost equal occupations, R = 50, no assignment possible (p.7084) | N−0.2; structure in doubt, see the note after the table | cases/lin5/v3 |

### Agreements and differences

Systems whose element-mean oxidation states agree:

- All of the first 15 except Ca3AuN and high-pressure K. YBaFe2O5 agrees only at the
  element-mean level; site charge order is a separate matter.
- MH4 hydrides: the metal oxidation state and the element mean of H agree with the
  paper. The site distribution differs for ScH4 and YH4: the paper has 2 Hb(−1) plus
  2 Ha(−½), LEOS has 3 H⁻ plus 1 neutral H.
- XeO, XeO2, XeF2/XeF3/XeF4, MgXe and two of the HgF4 structures. XeF3 also agrees at
  site level: the two kinds of Xe are +2 and +4 (cases/xef3/v3).
- The greedy baseline gives the metal +10, +9, +10 in SrH4, CeH4, ThH4 (the greedy
  column of EOSCAR); the flow gives +2, +4, +4.

Systems that differ:

| System | Paper | LEOS | Note |
|---|---|---|---|
| Ca3AuN | Au −3, R = 60.4 | Au −1 plus 2 free electrons | see "Features: electrides" above |
| High-pressure K | K +1, electrons booked on the NNAs | K averages 0 | LEOS has no bookkeeping unit like the NNA; the free electrons are booked on K, with a stray warning |
| XeO3 130 GPa | Xe +4, R = 60.9 | +5 | the alternative the paper gets with O2 as a fragment is also Xe +5, but with a lower R |
| XeO3 150 GPa | Xe1 +6, Xe2 +4, R = 51.2 | +6 | the paper considers the assignment of this phase inherently ambiguous |
| HgF4 (hgf4_b) | +4 | Hg averages +3 | the paper's three structures are all I4/m; how the three structures of this project map onto them is not verified, and the registry records hgf4_c as a Pm-3m-type framework |

Two more points affect the comparison:

- The LEOS cells of CsH3 and LiN5 differ from what the paper describes: linear H3
  units (p.7083) and the pentazolate anion N5⁻ (p.7084). In this project the H atoms
  of the CsH3 cell form two isolated equilateral H3 triangles, and the LiN5 cell
  holds linear azide N3 units plus a stretched N2, not a pentazole ring
  (docs/testset_registry.tsv, docs/testset_findings.md). The numbers in these two
  rows do not compare the same structure.
- When pyleos recomputes the 30 full MLWF runs of the Comas-Vilà comparison set, 12
  carry stray warnings. Stray counts: ba3si4 36, xeo3_b 9, ceh4 8, thh4 8, k_hp 6,
  mgxe 5, ybafe2o5_b 4, csh3 3, mgh4 3, sch4 3, lin5 1, xeo3_a 1.

To reproduce:

```bash
# one Comas-Vilà system: per-element assignment, stray/ghost warnings, R_elem/R_site
python3 pyleos/cli.py case cases/xeo3_a/v3
# 74-case benchmark (includes the 10 shared systems); --out must be given explicitly, otherwise the golden table is overwritten
python3 -m pyleos bench full --out /tmp/bench.tsv
```

## Index to the paper

| Content | Location |
|---|---|
| Molecular EOS (Ramos-Cordoba et al. 2015); applies to Hilbert-space and real-space partitions | end of Introduction and start of Theory, p.7076 |
| Definitions of eff-AO, EFO and AOM, eqs. 1–4 | Theory, p.7077 |
| Global sorting algorithm; no one-by-one rounding | Theory, p.7077 |
| Definition of R(%), eq. 5; meaning of R below 50 | Theory, p.7077–7078 |
| Equal sharing of electrons under degeneracy; R is not a probability; the form of eq. 5 is somewhat arbitrary | p.7078 |
| Program chain: PAW reconstruction, norm-conserving calculation, MLWFs, AOM under QTAIM; fragment AOMs add | start of Results, p.7078 |
| The solid-version R lets the frontier orbitals belong to one atom (supercell argument) | p.7078 |
| Test-set groups and basis of the references; SI Tables S1–S16 | p.7078 |
| Formulas and structure sources of the 40 systems | Table 1, p.7079 |
| CrO3: +6, QTAIM charge +1.76, R = 71 | p.7078–7079 |
| Perovskites 2–4; antiperovskites 5–6; Rb3AuO is a metal, Bloch states used instead | p.7079 |
| Ca3AuN: Au(−3), R = 60.4, no NNA found; the Au(−1) + 2e⁻ of Karen et al. | p.7079 |
| NNAs of high-pressure K; NNA sensitivity in Ca2N; the suggested ELF + Voronoi partition | p.7079 |
| GaSe: R = 50, Ga–Ga fragment R = 77.7 | p.7079–7080 |
| Ba3Si4, Cs2Pt, WCl4 (R = 92.6) | p.7080 |
| Rb2CuCl4 (β channel R = 75.4) | p.7080–7081 |
| YBaFe2O5 (R = 55.6, 61.3); the first 15 agree with iBOS/BVS | p.7081 |
| MH4 hydrides and the choice of fragment partition | p.7081–7082 |
| XeO, XeO2 (R = 74.6, 73.3, 77.6, 58.7) | p.7082 |
| XeO3, two phases (No. 28 R = 60.9; No. 29 R = 51.2) | p.7082–7083 |
| XeF2/XeF3/XeF4 (R 89–98); MgXe (R = 59.6); HgF4; CsH3, CsF3 | p.7083 |
| CsF5; LiN5 (R = 50) | p.7084 |
| Conclusions: extends to other real-space partitions; R used to self-assess the range of applicability; agreement with iBOS/BVS | Conclusions, p.7084 |
| Computational parameters and programs used | Computational Details, p.7084 |
