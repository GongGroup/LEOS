# Structure audit of the test set

**Status (audited 2026-07-12, issue #19; re-checked 2026-10-21 against the geometry
VASP actually reads):** the four tungsten chlorides have been recomputed with the
correct structures, and MgXe and the four MH4 hydrides now use the exact structures
from the literature. The re-check also found that the geometry of two v1 cases, li2o2
and re2, and of the delithiated Li2O2 supercells is broken; see "Geometry re-check"
below. The other problem cases keep their original structures, and their assignments
are read according to what the cell actually contains. Case-by-case records are in
[testset_registry.tsv](testset_registry.tsv); the source, verdict and evidence for
every structure are registered in [cases/STRUCTURES.tsv](../cases/STRUCTURES.tsv).

## Audit method

For all 113 test units (the 74-case benchmark, the Comas-Vilà comparison set and the
SCDM comparison copies), 12 parallel automated audit agents read each structure file,
compared the chemical formula, checked the reference oxidation states and searched
the literature for known crystal phases.

| Verdict | Now | At audit | Meaning |
|---|---|---|---|
| OK | 89 | 86 | composition matches and the source is credible, or a routine idealised approximation |
| SUSPECT | 14 | 15 | source unknown or chemically doubtful: hand-built cluster, hypothetical polymorph, deviation from the known phase |
| WRONG | 10 | 12 | the structure clearly conflicts with the formula or the known phase, or the geometry itself is broken |

The "now" counts are the current verdicts of the registry: WCl3–6 went from WRONG to
OK after switching to the corrected structures; li2o2 went from SUSPECT, and re2 from
OK, to WRONG (see "Geometry re-check").

Registry columns: case, batch, dir, formula_expected (the composition the case name
implies), formula_cell (the actual composition of the cell), stoich_match,
known_phase (phase known in the literature), source_guess (inferred source of the
structure), ref_os_entry (entry in the reference oxidation states), ref_os_ok
(whether the reference is charge-balanced), verdict, note.

## Problem types

### 1. CIF transcribed as the asymmetric unit only (6 cases)

The lattice parameters match experiment, but the cell holds only the atoms of the
asymmetric unit; the space group was not applied.

| Case | Actual cell | Should be |
|---|---|---|
| wcl3 | 1 W + 3 Cl (R-3 phase) | 24 atoms (W6Cl18). The stoichiometry happens to be right; the volume of 441 Å³/atom gives the problem away |
| wcl4 | 1 W + 3 Cl (C2/m) | 4 W + 16 Cl |
| wcl5 | 2 W + 9 Cl (C2/m) | 48 atoms |
| wcl6 | 1 W + 1 Cl (one W at 3a and one Cl at 18f) | 21 atoms |
| scdm_wcl4/5/6 | the same structures as above | |

### 2. Species in the cell do not match the case name (5 cases)

| Case | Actual cell |
|---|---|
| CO (legacy molecule set) | two ethylenedione molecules O=C–C=O (C–C 1.53 Å), at about twice the density of solid CO |
| hi3 | van der Waals chains of HI and I2 molecules (H–I 1.70 Å, I–I 2.78 Å), not I3⁻ |
| lin5 | linear N3 units plus one stretched N2, not a pentazole ring |
| rhcl_co__p_c6h5_3_2 | one extra coordinated SO2, so the cell is RhCl(CO)(SO2)(PPh3)2; the reference is written for the compound without SO2 |
| AuClB complex | the cell holds C42H32, which does not match the name; one benzene ring is bonded to gold through a C–Au σ bond |

### 3. Hand-built idealised cells or hypothetical polymorphs (8 cases)

| Case | Problem | Verdict |
|---|---|---|
| sno | hypothetical rock-salt phase; SnO really exists only as the layered litharge phase | WRONG |
| ti2o3 | a 5-atom defect model, rutile TiO2 with one O removed, not corundum-type Ti2O3 | WRONG |
| cuo | hypothetical rock-salt polymorph; CuO is really monoclinic tenorite, and the Jahn–Teller coordination is wiped out | SUSPECT |
| ce2o3 | idealised cell; the Ce coordination number does not match the A-type phase, and the c axis is 7% too long | SUSPECT |
| reo2 | a rutile template with rounded parameters, matching neither the ground-state monoclinic phase nor the experimental phase; reference missing | SUSPECT |
| n4s4 | density 25% below that of the known monoclinic phase | SUSPECT |
| cs2pt | c/a 16.6% shorter than experiment, volume 12.7% smaller | SUSPECT |
| csh3 | the H atoms form isolated H3 triangles, unlike the linear picture in the literature | SUSPECT |

### 4. Unrelaxed isolated clusters (5 cases, legacy oxyanion salts)

na2cro4, na3vo4, namno4 and recl4 are hand-built, unrelaxed clusters whose bond
lengths are off by 26–35% (for example Cr–O 2.08 Å against a real 1.66–1.68 Å). An
isolated multiply charged anion cannot bind its electrons in the gas phase, so the
electronic structure of the high oxidation states is unreliable. li2o2 (legacy oxide
set) is an isolated molecule in a box, not a crystal. nareo4 is also a gas-phase ion
pair, but its bond lengths are right and it is judged OK.

### Missing reference oxidation states

na3vo4 and reo2 have no entry in the reference table, so they are not among the
scored cases.

## Actions taken

| Case | Verdict | Action | Current benchmark result |
|---|---|---|---|
| wcl3, wcl4, wcl5, wcl6 | WRONG | recomputed after applying the space group or switching to the literature structure; cases/wcl*/v3 now hold the new structures | WCl3, WCl4, WCl6 pass; WCl5 fails (wannierisation failed) |
| scdm_wcl4/5/6 | WRONG | not recomputed (copies used for the SCDM-only comparison) | not counted in the 74-case benchmark |
| CO | WRONG | recomputed with α-CO (mp-11875), cases/co/v31 | original and new runs both pass |
| sno, ti2o3 | WRONG | recomputed with romarchite (COD 9011217) and corundum-type Ti2O3 (COD 9008082), `cases/<system>/v31` | original and new runs both pass |
| rhcl_co__p_c6h5_3_2 | WRONG | kept | fails (dev 4.00) |
| AuClB complex | WRONG | kept | no reference, not scored |
| hi3, lin5, csh3 | SUSPECT | kept: HI3 has no solid phase; LiN5 and CsH3 can be compared only with the predicted structures in the SI of the Comas-Vilà paper, which neither COD nor MP has (the current CsH3 cell is mp-1006815) | hi3 fails; lin5 and csh3: see the Comas-Vilà comparison set |
| cuo, ce2o3, recl4, na2cro4 | SUSPECT | recomputed with tenorite (COD 9014580), A-type Ce2O3 (mp-2721), β-ReCl4 (COD 4123951) and Cmcm Na2CrO4 (COD 2107101), `cases/<system>/v31` | original and new runs both pass |
| namno4, n2o, n4s4 | SUSPECT | recomputed with AgMnO4-type NaMnO4 (mp-765753, isotypic with experimental AgMnO4, COD 2105345), P2₁3-ordered N2O built on the experimental Pa-3 cell, and experimental S4N4 (COD 9008994), `cases/<system>/v31` | original and new runs both pass |
| cs2pt | SUSPECT | recomputed with the experimental lattice (Karpov 2003, a = 5.676, c = 9.471 Å), cases/cs2pt/v31; the original cell is mp-13548 | both give Cs +1, Pt −2 |
| li2o2, re2 (geometry re-check) | WRONG | li2o2 now uses the Föppl crystal, cases/li2o2/v3; re2 recomputed as the intended 2.40 Å dimer, cases/re2/v31 | original and new runs both pass |
| other SUSPECT | SUSPECT | kept, flagged in the data sources | |

In addition, MgXe and the MH4 hydrides (mgh4, srh4, ceh4, thh4) of the Comas-Vilà
comparison set first used approximate 0 GPa structures, which the audit judged OK.
They were later switched to the exact literature structures (SI of Bi & Zurek 2021;
for MgXe a lattice estimated from Miao 2015) and recomputed; the data are in
`cases/<system>/v3`.

The source, space group and self-check of every replacement structure are recorded in
`cases/<system>/structure-fixed/SOURCE.md`; the summary table is
results/structures_fixed_manifest.tsv.

**Cases in the golden table (original runs of the 74 cases) that still have wrong
structures:** CO, sno, ti2o3, li2o2 and re2, plus the doubtful cuo, ce2o3, recl4,
na2cro4, namno4, n2o and n4s4. All have been recomputed on the correct structures and
pass; see RESULTS.md, "Structure-corrected version", and
results/bench74_corrected.tsv. When benchmark results are quoted, the two tables can
be given side by side.

## Geometry re-check (2026-10-21)

A scan of the shortest interatomic distance in every structure, using the geometry
VASP actually reads (with Cartesian coordinates in POSCAR, the coordinates too are
multiplied by the scale factor), found three broken geometries:

| Structure | Problem | Evidence | Effect |
|---|---|---|---|
| li2o2/v1 | scale factor 10 applied to the Cartesian coordinates: the two Li coincide at (0, 5, 0) and the two O are 5 Å apart, not a peroxide | ionic positions in cases/li2o2/v1/OUTCAR | counted as a pass in the 74 cases (v1_oxide) |
| re2/v1 | scale factor 12 applied to the Cartesian coordinates: the two Re are 4.8 Å apart, not a 2.4 Å dimer | cases/re2/v1/OUTCAR | counted as a pass in the 74 cases (v1_oxide) |
| v3 of li16o16, li15o16 and li14o16, and li15o16/v3-md | the second set of Li in the 2×2×1 supercell sits on the columns of the peroxide ions, at the O–O bond midpoints (Li–O 0.78 Å) | POSCAR; single-point calculations, NSW = 0 | the three-row table of RESULTS.md, "Li2O2 delithiation", uses this geometry; v3-md is an MD attempt that crashed after the first step |

The original audit read li2o2 and re2 with unscaled coordinates, so the bond lengths
it derived do not hold; both rows are corrected in the registry. All three have been
recomputed on correct geometries (see "Actions taken" above and RESULTS.md). The
v3-u, v3-md2 and polaron series of li15o16 come from a structure obtained by 31
relaxation steps starting from the broken geometry
(SAI: LEOS/inputs/m5a/li15o16_relax); it matches Föppl-type Li2O2 and its geometry is
correct.

The legacy batch (v1) consists mostly of hand-built structures; the new batch (the
Comas-Vilà comparison set) mostly of literature and database structures, of clearly
better quality.
