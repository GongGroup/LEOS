# THH4: I4/mmm (ThCr2Si2-type) ThH4, cell optimised at 0 GPa

- Source: Bi & Zurek, *Chem. Eur. J.* **27**, 14858 (2021),
  DOI 10.1002/chem.202102679, Supporting Information, Section S2
  "Structural Information for Metal Tetrahydrides at 0 GPa"
  (SI PDF obtained from the Wayback Machine snapshot of 2025-07-10).
- Space group: I4/mmm (No. 139). Conventional body-centred tetragonal cell, Z=2,
  10 atoms.
- Lattice parameters (0 GPa, PBE-optimised): a = b = 3.588 Å, c = 6.348 Å.
- Wyckoff positions: Th 2a (0,0,0); H 4d (0,1/2,1/4); H 4e (0,0,z), z = 0.3689.
- Pressure point: **0 GPa**. The phase is thermodynamically stable at high pressure;
  Bi & Zurek optimise it statically at 0 GPa and analyse the bonding. The Comas-Vilà
  JCTC 2025 comparison set (cases 16-23) uses these 0 GPa structures: the H-H distance
  in its SI, 1.664 Å, equals |1-2z|·c = 1.664 Å of this structure exactly. The
  structure is therefore taken at 0 GPa, so that the comparison is like for like.
- Self-checks: stoichiometry Th2H8 = 2×ThH4 ✓; H-H of the H2 molecular unit = 1.664 Å
  (= literature value) ✓; nearest Th-H = 1.664 Å; the element order follows the
  existing order of data/v3/comas_gap/thh4/POSCAR (H Th), so that the POTCAR assembly
  logic can be reused.
- Replaces: data/v3/comas_gap/thh4/POSCAR (the placeholder structure, registered as an
  approximate phase at 0 GPa).
