# WCL5: monoclinic C2/m W2Cl10 dimer molecular crystal, complete cell (72 atoms)

- Phase: solid WCl5 = W2Cl10 dimers of two edge-sharing octahedra, monoclinic C2/m
  (No. 12). Experimental determination: Cotton & Rice, *Acta Cryst.* **B34**,
  2833 (1978): a=17.438, b=17.706, c=6.063 Å, β=95.51° (= ICSD 2398; the original cell
  of the AFLOW entry Cl10W2_ICSD_2398 agrees with it exactly).
- **Atom-count correction**: this cell has Z=12 WCl5 (6 W2Cl10 dimers) = **72 atoms**
  (W12Cl60), not the 48 written in the registry and the task sheet. Independent check:
  the volume of the existing cell gives 1981 Å^3 / 72 = 27.5 Å^3/atom, in agreement
  with the 27.3 Å^3/atom of the existing WCl6 cell; with 48 atoms it would be
  41 Å^3/atom, which is not a reasonable density. The asymmetric unit = 2 W + 9 Cl =
  11 sites, exactly the 11 lines of the existing POSCAR (the original transcription
  took the asymmetric unit for the full cell).
- Lattice: the existing cell of the test set is kept (a=17.779, b=17.980, c=6.231 Å,
  β=96.08°; the registry judges it correct).
- Atomic coordinates: the relaxed structure of the AFLOW ICSD_WEB entry
  Cl10W2_ICSD_2398 (36-atom monoclinic C-centred primitive cell) is transformed
  numerically to the C2/m conventional cell (72 atoms), and the fractional coordinates
  are then placed in the existing lattice.
  Orbits: W 4g + W 8j; Cl 4i×3 + Cl 8j×6 (4+8=12 W; 12+48=60 Cl) ✓.
- Lineage check: the existing 11 atoms are compared one by one with the 72-atom
  template (a C2/m origin shift is allowed); 7/11 match within 0.8 Å, with a mean
  deviation of 0.67 Å. This agrees with the judgement that the incomplete cell (the
  asymmetric unit) had also been relaxed.
- Self-checks: stoichiometry W12Cl60 = 12×WCl5 ✓; all 12 W are six-coordinate ✓;
  W-Cl = 2.139-2.437 Å; shortest interatomic distance 2.139 Å > 1.5 Å ✓. Remark: the
  AFLOW relaxed cell (PBE without dispersion correction) is 3-6% larger than the
  existing cell. After the fractional coordinates are transplanted, the intramolecular
  W-Cl is systematically about 4% shorter than experiment
  (terminal 2.24-2.28/bridging ~2.4-2.5 Å), and the shortest terminal bond, 2.14 Å, is
  slightly below the 2.2 Å guideline. A geometry difference of this size does not
  affect the oxidation-state analysis.
- Applies to cases: wcl5 (m4_wave2) and scdm_wcl5 (comas). Pressure: ambient-pressure
  experimental phase.
