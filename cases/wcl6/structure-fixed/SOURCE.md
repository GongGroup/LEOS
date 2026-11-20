# WCL6: trigonal R-3 α-WCl6, complete hexagonal cell (21 atoms)

- Phase: α-WCl6 molecular crystal, space group R-3 (No. 148), hexagonal setting, Z=3
  (same phase as Materials Project mp-571518; consistent with known_phase in the
  registry).
- Lattice: the existing cell of the test set is kept (data/v3/comas_gap/scdm_wcl6 and
  data/v3/m4_wave2/wcl6, byte-identical): a = 6.2356 Å, c = 17.0110 Å (experiment
  6.088/16.68; the existing cell is about +2%; the registry judges the lattice
  correct).
- Atomic coordinates: W 3a (0,0,0); Cl 18f (0.295, 0.295, 0.08), taken from the
  experimental CIF COD 1537510 (Ketelaar & van Oosterhout, Recl. Trav. Chim. Pays-Bas
  **62**, 197 (1943), = ICSD 26081) and expanded to the complete hexagonal cell with
  all 18 symmetry operations of R-3.
- Why the Cl coordinates of the existing POSCAR are not used: the Cl of the existing
  2-atom cell has drifted about 0.98 Å from the standard 18f orbit (a trace of the
  2-atom incomplete cell having been relaxed), and direct expansion would give distorted
  octahedra; Cl therefore takes the experimental coordinates. The offset of the existing
  W position from 3a is likewise drift (about 1.6 Å); the expansion uses the standard
  origin.
- Self-checks: stoichiometry W3Cl18 = 3×WCl6 ✓; each W is coordinated by exactly 6 Cl,
  W-Cl = 2.288 Å (experiment 2.24-2.29 Å; 2.29 Å in the existing, slightly expanded
  cell) ✓; any two atoms are ≥ 2.288 Å > 1.5 Å apart ✓.
- Applies to cases: wcl6 (m4_wave2) and scdm_wcl6 (comas). Pressure: ambient-pressure
  experimental phase.
