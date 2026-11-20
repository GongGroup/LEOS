# WCL3: trigonal R-3 W6Cl18 cluster phase, rhombohedral primitive cell (24 atoms)

- Phase: "WCl3" is in fact the molecular cluster compound W6Cl18 ([W6Cl12] inner
  cluster + 6 terminal Cl), space group R-3 (No. 148). Experimental structure:
  ICSD 23129 (hexagonal cell a=14.91, c=8.47 Å; obtained through the AFLOW entry
  Cl18W6_ICSD_23129, whose original cell agrees with ICSD).
- Lattice: the existing hexagonal cell of the test set is kept, a=15.155, c=8.873 Å
  (the registry judges it consistent with the R-3 lattice of W6Cl18). As the task
  requires, the output is the **complete 24-atom cell**, i.e. the rhombohedral primitive
  cell of the R-centred hexagonal cell (a_rh = 9.2363 Å, α = 110.254°, converted exactly
  from the existing hexagonal cell), which holds 1 complete W6Cl18 cluster. The
  corresponding complete hexagonal conventional cell has 72 atoms (Z=3); the two are
  equivalent descriptions of the same crystal.
- Atomic coordinates: the relaxed AFLOW structure of ICSD 23129 (24-atom rhombohedral
  primitive cell) is expanded to the hexagonal cell (72 atoms) and folded back into the
  rhombohedral primitive cell of the existing lattice. Every step is a numerical
  transformation; nothing is adjusted by hand.
- Lineage check: a rigid-translation fit of the 4 atoms of the existing POSCAR (the
  1W+3Cl asymmetric unit) to the 72-atom template gives the optimal translation (0,0,0)
  and a mean deviation of 0.85 Å. The existing 4 atoms are in the same coordinate frame
  but have drifted under relaxation, so they should not be expanded directly; the
  database coordinates are used instead.
- Self-checks: stoichiometry W6Cl18 (= 6×WCl3) ✓; the W6 octahedron has
  W-W = 2.803 / 2.926 Å in two alternating sets (trigonal-antiprismatic distortion;
  experiment is of the order of 2.85-2.95 Å) ✓; each W is coordinated by
  4 bridging Cl + 1 terminal Cl, W-Cl = 2.314-2.333 Å ✓; shortest interatomic distance
  2.314 Å > 1.5 Å ✓ (the coordinates are transplanted from a PBE-relaxed cell about
  3.8% larger, so the bond lengths are correspondingly slightly shorter than
  experiment).
- Applies to case: wcl3 (m4_wave2). Pressure: ambient-pressure experimental phase.
