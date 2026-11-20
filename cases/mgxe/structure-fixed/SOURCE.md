# MGXE: CsCl-type (B2) MgXe, 200 GPa

- Phase and pressure: Miao et al., *J. Am. Chem. Soc.* **137**, 14122 (2015),
  DOI 10.1021/jacs.5b08162: MgXe enters the convex hull at ≥125 GPa; at >100 GPa the
  stable structure is CsCl-type (Pm-3m). Comas-Vilà JCTC 2025 (case 33) does the EOS
  analysis at **200 GPa** in Pm-3m; this structure matches that setting.
- Space group: Pm-3m (No. 221). Mg 1a (0,0,0), Xe 1b (1/2,1/2,1/2), Z=1, 2 atoms.
- Estimate of the lattice constant (the literature gives no table of values; the ACS SI
  could not be obtained):
  V(MgXe, 200 GPa) = V_Mg + V_Xe + ΔV
  where V_Mg(200 GPa) = 9.07 Å^3 (Vinet EOS of bcc Mg: V0=23.1398 Å^3, K0=28.9 GPa,
  K'=4.80; Stinton et al., PRB **90**, 134105 (2014));
  V_Xe(200 GPa) = 15.66 Å^3 (Vinet EOS of hcp Xe: V0=57.41 Å^3/atom, K0=4.887 GPa,
  K'=6.2955, valid for 10-260 GPa; fit of Dewaele et al. 2012);
  ΔV(200 GPa) ≈ -1.9±0.3 Å^3, read from Figure 3b of the Miao paper (the curve
  ΔV=V_MgXe-V_Mg-V_Xe).
  => V = 22.83 Å^3, a = 2.837 Å (uncertainty about ±0.04 Å, i.e. ~1.5%).
- Consistency check: the text of the paper says that at 125 GPa the volume of MgXe is
  about 10% smaller than the sum Mg+Xe
  => ΔV(125) ≈ -2.88 Å^3, consistent with the Figure 3b reading
  (about -2.6 to -2.9 Å^3) ✓.
- Self-checks: stoichiometry Mg1Xe1 ✓; Mg-Xe = 2.457 Å (= a·√3/2, eight-coordinate
  CsCl geometry); Xe-Xe = a = 2.837 Å, comparable to the nearest-neighbour distance in
  elemental Xe at 200 GPa (~2.8 Å) and consistent with the paper's description of a
  strong Xe-Xe interaction.
- Note: a is an estimate, not a value transcribed from the literature. It is the only
  lattice constant in this package that is not transcribed directly. If the Miao SI
  (ja5b08162_si_001.pdf) is obtained later, replace a with its tabulated value.
