# Fig. 1 - caption draft and provenance

Files: `fig1_method.py` -> `fig1_method.pdf` (vector, 17.2 x 7.1 cm, Arial
embedded) and `fig1_method.png` (300 dpi preview).
Regenerate with `python3 fig1_method.py` from this directory (about 30 s; the
script prints the provenance table that the numbers below are copied from).

## Caption

**Fig. 1** Assignment of Wannier functions (WFs) to atoms. **(a)** Atoms
(labelled circles) and WF centres of an H<sub>2</sub>S<sub>2</sub> molecule,
orthographic projection. Filled circles: bond-pair WFs, i.e. centres that lie
between their two nearest atoms and within 0.10 Å of the line joining them;
open circles: all other centres (lone pairs). Each WF holds two electrons.
Lines join the atoms bridged by a bond-pair WF. Scale bar: 1 Å in the plane of
projection. **(b)** Position *f* of the bond-pair WF centre along each terminal
X–H bond (*f* = 0 at X, *f* = 1 at H) in the 12 benchmark compounds that
contain such bonds, grouped by the heavy atom X and ordered by the Allen
electronegativity difference χ<sub>X</sub> − χ<sub>H</sub> (in parentheses);
*n* is the number of bonds. Every dot is one bond; within a row the compounds
are offset vertically. Dashed line: *f* = 0.5, where a nearest-atom criterion
transfers the electron pair from X to H. The horizontal line separates
χ<sub>X</sub> > χ<sub>H</sub> (the pair belongs to X in the ionic
approximation) from χ<sub>X</sub> < χ<sub>H</sub> (the pair belongs to H);
in the shaded regions the nearest-atom criterion contradicts this rule.
**(c)** The assignment for the molecule in (a) as a flow network. Each WF *i*
receives its electrons from the source *s* (2 *e* per WF; 1 *e* per WF and spin
channel in spin-polarised calculations), sends them to atoms *j* at a cost
*c<sub>ij</sub>* per electron (*d<sub>ij</sub>*: distance between WF centre and
atom; χ<sub>j</sub>: Allen electronegativity; χ<sub>max</sub>: largest
electronegativity in the cell; β = 10 Å<sup>2</sup>), and each atom passes at
least *n<sub>j</sub>*<sup>min</sup> and at most *n<sub>j</sub>*<sup>max</sup>
electrons to the sink *t* (bounds in brackets); the total flow does not exceed
the number of valence electrons, *N*<sub>elect</sub> = 14. WFs and atoms are
listed from top to bottom as they occur in (a), with the markers of (a). Blue
arrows: edges that carry flow in the minimum-cost solution; grey lines: edges
without flow; dashed arrow: alternative of equal cost for the S–S bond pair.
Both S–H bond pairs go to S (*c* = 0.80 Å<sup>2</sup>) although H is the
nearer atom (*c* = 1.33 Å<sup>2</sup>).

## Numbers behind the caption (all printed by the script)

Panel (b), top to bottom; lanes inside a row in the order listed.

| Bond | χ<sub>X</sub> − χ<sub>H</sub> (Allen) | *n* | Compound (*n*, range of *f*) | *f* > 0.5 |
|---|---|---|---|---|
| O–H | +1.310 | 88 | HNO<sub>3</sub> (8, 0.470); H<sub>3</sub>PO<sub>3</sub> (16, 0.445–0.456); SrBe(OH)<sub>4</sub> (64, 0.536–0.541) | 64 |
| N–H | +0.766 | 12 | NH<sub>3</sub>BH<sub>3</sub> (12, 0.560–0.564) | 12 |
| S–H | +0.289 | 6 | H<sub>2</sub>S<sub>2</sub> (2, 0.658); H<sub>2</sub>S<sub>3</sub> (2, 0.657–0.658); H<sub>2</sub>S<sub>5</sub> (2, 0.655–0.656) | 6 |
| C–H | +0.244 | 41 | CH<sub>2</sub>F<sub>2</sub> (2, 0.670); (CH<sub>3</sub>)<sub>3</sub>NO (9, 0.654–0.669); RhCl(CO)(SO<sub>2</sub>)(PPh<sub>3</sub>)<sub>2</sub> (30, 0.649–0.668) | 41 |
| I–H | +0.059 | 2 | HI·I<sub>2</sub> chain structure, case `hi3` (2, 0.602) | 2 |
| P–H | −0.047 | 8 | H<sub>3</sub>PO<sub>3</sub> (8, 0.719–0.725) | 8 |
| B–H | −0.249 | 36 | B<sub>6</sub>H<sub>10</sub>, terminal H only (24, 0.720–0.734); NH<sub>3</sub>BH<sub>3</sub> (12, 0.773–0.785) | 36 |

- χ<sub>X</sub> > χ<sub>H</sub>: 125 of 149 bond pairs lie nearer to H than to X
  (all N–H, S–H, C–H, I–H; 64 of 88 O–H). χ<sub>X</sub> < χ<sub>H</sub>: 44 of 44.
- Panel (c): costs of an S–H pair are 0.797 (to S) and 1.331 (to H); of the
  S–S pair 1.056 to either S (difference < 10<sup>−4</sup>); of a lone pair
  0.226–0.227 to its S. Capacities H [0, 2], S [0, 8]. Minimum total cost 7.1108.

## Provenance

- **Panels (a), (c)**: `cases/h2s2/v3/POSCAR` and `EOSCAR` (#WF block:
  7 centres; header: NELECT(Wannier) = 14). Read with `pyleos.parsers.load_case`;
  capacities from `pyleos.capacity.capacity_bounds`, electronegativities from
  `pyleos.capacity.ENEG_ALLEN`, optimal edges from
  `pyleos.transport.assign_case(case, beta=10)`. The script recomputes
  *c<sub>ij</sub>* from the formula in the figure and asserts that it
  reproduces the solver's objective.
- **Panel (b)**: the cases of the 74-case benchmark
  (`pyleos/sets/reference_os.json`, directories resolved through
  `cases/INDEX.tsv` by `pyleos.bench.collect`, as in `pyleos bench full`) that
  have a reference answer and contain H:
  `cases/{b6h10,bh6n,ch2f2}/v1` (centres from `wannier90.wout`),
  `cases/{hno3,tmao}/v3` (benchmark set `wave1`) and
  `cases/{h2s2,h2s3,h2s5,h3po3,hi3,rhcl_co_p_c6h5_3_2,rhh,srbe_oh_4}/v3`
  (set `wave2`; centres from the EOSCAR #WF block). All are
  non-spin-polarised.

## Selection rules (constants at the top of the script)

1. *Terminal X–H bond*: X is the non-H atom nearest to the H atom (all periodic
   images count); the bond is used if the next non-H neighbour is farther than
   1.25 × *d*(X–H). Accepted bonds have ratios ≥ 1.43; rejected H atoms have
   ratios ≤ 1.05, so the result does not depend on the threshold.
2. *Bond-pair WF of that bond*: the WF centre nearest to the H atom, accepted
   if it lies within 0.10 Å of the X–H axis and 0 < *f* < 1. All terminal bonds
   pass; the largest offset is 0.055 Å (O–H in H<sub>3</sub>PO<sub>3</sub>).
3. *f* = (**r**<sub>WF</sub> − **r**<sub>X</sub>)·(**r**<sub>H</sub> − **r**<sub>X</sub>) / |**r**<sub>H</sub> − **r**<sub>X</sub>|<sup>2</sup>, minimum-image convention.
4. *Panel (a) classes*: same axis rule applied to the two atoms nearest to each
   centre. Offsets: 0.006 Å (S–H), 0.035 Å (S–S); lone pairs 0.45–0.46 Å.
5. *View of panel (a)*: chosen by the script among 20 000 directions: all drawn
   objects kept apart, then the least foreshortened bonds (the shortest-looking
   bond is drawn at 0.79 of its length). Atom circles are symbols, not radii.
6. *Atom labels* H1, H2, S1, S2 number the atoms of each element in POSCAR order.

## Not drawn, and why

- **Rh–H** (case `rhh`): H occupies an interstitial site with six equidistant
  Rh neighbours (2.03 Å); there is no two-centre bond, and the WF centre
  nearest to H is 1.35 Å away from it and off every Rh–H axis (≥ 0.37 Å).
- **Bridging H of B<sub>6</sub>H<sub>10</sub>** (16 atoms): three-centre B–H–B
  bonds, excluded by rule 1.
- `cases/auclb_c6h5_3_p_c6h5_3_2/v1` (32 C–H bonds, *f* = 0.652–0.669):
  no reference answer, hence not a scored case.

## Caveats for the text

- **O–H is split.** Only the hydroxide SrBe(OH)<sub>4</sub> (O–H 0.97 Å) has
  its centre on the H side (*f* = 0.54). In the two hydrogen-bonded acids
  (O–H 1.02–1.03 Å, H···O 1.53–1.62 Å) the centre is on the O side
  (*f* = 0.445–0.470), and a nearest-atom criterion assigns these 24 bonds
  correctly. The statement "O–H centres lie nearer to H" holds for 64 of 88 bonds.
- **"Nearest-atom criterion" here is the pure distance rule** (each WF to its
  nearest atom). It is not the in-situ baseline stored in the EOSCAR column
  OS(int) (`results/bench74_greedy.tsv`), which gives an off-site pair to the more
  electronegative of its two nearest atoms and then caps O at 8 electrons.
  H atoms assigned −1 in the 13 surveyed cases: pure nearest-atom rule
  185 of 210; EOSCAR baseline 102 (B<sub>6</sub>H<sub>10</sub> 40,
  NH<sub>3</sub>BH<sub>3</sub> 12, SrBe(OH)<sub>4</sub> 50); min-cost flow 63
  (B<sub>6</sub>H<sub>10</sub> 40, NH<sub>3</sub>BH<sub>3</sub> 12,
  H<sub>3</sub>PO<sub>3</sub> 8, `hi3` 2, RhH 1).
- **I–H and P–H are convention-dependent** (|χ<sub>X</sub> − χ<sub>H</sub>| < 0.06
  on the Allen scale; on the Pauling scale I–H is +0.46 and P–H −0.01).
  `clean_set.json` excludes `h3po3` and `hi3` for this reason; `hi3` is in
  addition flagged SUSPECT (the cell is HI + I<sub>2</sub> chains), and the
  Rh complex is the SO<sub>2</sub> adduct, not the compound it is named after.
  These flags concern the reference answers, not the WF geometry shown here.
  The min-cost flow itself gives H(−1) in `hi3` (a known FAIL in
  `results/bench74_golden.tsv`).
- **The example in (c) is degenerate in the S–S pair**: the solver returns
  S1(−2), S2(0) (element mean −1, the reference); sending the pair to S2 or
  one electron to each S costs the same. The dashed arrow shows this.
- A suitable citation for the failure of nearest-atom counting for X–H bonds
  is the one listed in `docs/PLAN.md` (Gimferrer, Comas-Vilà, Salvador,
  *Molecules* 2020).
