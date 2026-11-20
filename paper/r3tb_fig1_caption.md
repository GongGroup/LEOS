# Method figure of the in-situ electronegativity (r3-T-B) - caption draft and provenance

Files: `r3tb_fig1_method.py` -> `r3tb_fig1_method.pdf` (vector, 17.2 x 7.1 cm,
Arial embedded) and `r3tb_fig1_method.png` (300 dpi preview). Counterpart of
`fig1_method.py` for the plugin `pyleos/methods/frozen/r3-T-B-cap-repair.py`.
Regenerate with `python3 r3tb_fig1_method.py` from this directory (about 8 s;
the script prints the provenance table that the numbers below are copied from
and keeps it in `_cache/r3tb/fig1_method.run.log`). Each panel is a function
(`draw_molecule(ax)`, `draw_dividing_point(ax)`, `draw_rule(ax)`) whose text
and line weights scale with the axes.

## Caption

**Fig.** Assignment of bond electrons by the in-situ electronegativity.
**(a)** Atoms (labelled circles) and Wannier-function (WF) centres of an
H<sub>2</sub>S<sub>2</sub> molecule, orthographic projection, in the two
classes of the method. Filled circles: contested WFs (bond pairs), whose
centre is nearer the midpoint of two atoms than either nucleus; open circles:
own WFs (lone pairs), which stay on their atom. Each WF holds two electrons.
Lines join the two atoms nearest to each bond pair. Scale bar: 1 Å in the
plane of projection. **(b)** Position *f* of the WF centre along each terminal
X–H bond (*f* = 0 at X, *f* = 1 at H) in the nine compounds of the clean38
benchmark that contain hydrogen, grouped by the heavy atom X; in parentheses
the effective nuclear charge *Z*<sub>eff</sub> of X (*Z*<sub>eff</sub> of H is
0.69); *n* is the number of bonds. Every symbol is one bond; within a row the
compounds are offset vertically. Dashed line: bond midpoint, *f* = 0.5, the
dividing point of a nearest-atom rule. Solid line: force-balance point of the
method, *f*\* = √*Z*<sub>eff,X</sub> / (√*Z*<sub>eff,X</sub> +
√*Z*<sub>eff,H</sub>), where *Z*<sub>eff,X</sub>/*d*<sub>X</sub><sup>2</sup> =
*Z*<sub>eff,H</sub>/*d*<sub>H</sub><sup>2</sup>. A centre left of *f*\* goes to
X and H is +1; a centre in the shaded region right of *f*\* goes to H and H is
−1. Filled circles: contested WFs, decided by *Z*<sub>eff</sub>/*d*<sup>2</sup>;
open circles: one-centre WFs of H in NH<sub>3</sub>BH<sub>3</sub>
(*f* > 0.75), which stay on H. 129 of the 137 centres lie nearer H than X. The
method gives 101 pairs to X (O–H, S–H, N–H, C–H) and 36 pairs to H (B–H), and
all 153 H atoms of clean38 receive their reference oxidation state.
**(c)** The decision rule. Top: the quantities, all read from the WFs of the
calculation. An atom keeps its own WFs: dipole-free shells (grey ring; centre
within 1/8 of the rms radius from the nucleus) and lone pairs (open circles),
whose second moment about the nucleus is the squared valence radius
*ρ*<sup>2</sup>. A contested WF (filled circle) lies at distances
*d*<sub>A</sub> and *d*<sub>B</sub> from its two ends. *Z*<sub>val</sub>:
valence charge of the pseudopotential; *n*<sub>inner</sub>: electrons of the
atom's own WFs inside the sphere of radius *d*; 5/16: screening by the partner
electron of the pair. Bottom: the pair goes to the end with the larger in-situ
electronegativity *χ* = *Z*<sub>eff</sub>/*ρ*<sup>2</sup> when both ends carry
lone pairs, and to the end with the larger force
*Z*<sub>eff</sub>/*d*<sup>2</sup> otherwise. Two examples with the numbers of
the method (*Z*<sub>eff</sub> in *e*; *χ* and
*Z*<sub>eff</sub>/*d*<sup>2</sup> in *e* Å<sup>−2</sup>): an As–S bond pair of
As<sub>4</sub>S<sub>4</sub> and an S–H bond pair of H<sub>2</sub>S<sub>2</sub>.
The shaded row takes the pair.

## Numbers behind the caption (all printed by the script)

clean38 (`python3 -m pyleos bench clean --split dev`): 38 cases, 793 atoms,
2768 WFs. Classes of the plugin: 2150 own WFs (1931 lone pairs, 219
dipole-free) and 618 contested WFs. Of the contested WFs, 566 have ends of two
elements: 527 decided by *Z*<sub>eff</sub>/*d*<sup>2</sup>, 39 by *χ* (32 As–S
in As<sub>4</sub>S<sub>4</sub>, 7 Mo–Cl in MoCl<sub>2</sub>); 52 have ends of
one element only. Atoms correct 758/769, cases fully correct 36/38
(`results/clean38_r3-T-B.tsv`, reproduced case by case by the script).

Panel (a): 7 WFs. Contested: two S–H pairs (centre 0.464 Å from H, 0.893 Å
from S; both go to S) and the S–S pair (1.028 Å from each S; shared 0.500 /
0.500). Own: four lone pairs, 0.475–0.477 Å from their S, rms radius
1.052–1.053 Å, offset/rms 0.452–0.453 (dipole-free below 0.125). Oxidation
states H +1, S −1.

Panel (b), top to bottom; lanes inside a row in the order listed.

| Bond | *Z*<sub>eff,X</sub> | *f*\* | *n* | Compound (*n*, range of *f*) | nearer H than X | pair to X | pair to H | smallest force ratio |
|---|---|---|---|---|---|---|---|---|
| O–H | 5.6875 | 0.7420 | 72 | HNO<sub>3</sub> (8, 0.470); SrBe(OH)<sub>4</sub> (64, 0.536–0.541) | 64 | 72 | 0 | 5.968 |
| S–H | 5.6875 | 0.7420 | 6 | H<sub>2</sub>S<sub>2</sub> (2, 0.658); H<sub>2</sub>S<sub>3</sub> (2, 0.657–0.658); H<sub>2</sub>S<sub>5</sub> (2, 0.655–0.656) | 6 | 6 | 0 | 2.230 |
| N–H | 4.6875 | 0.7231 | 12 | NH<sub>3</sub>BH<sub>3</sub> (12, 0.560–0.564) | 12 | 12 | 0 | 4.069 |
| C–H | 3.6875 | 0.6984 | 11 | CH<sub>2</sub>F<sub>2</sub> (2, 0.670); (CH<sub>3</sub>)<sub>3</sub>NO (9, 0.654–0.669) | 11 | 11 | 0 | 1.297 |
| B–H | 2.6875 | 0.6641 | 36 | B<sub>6</sub>H<sub>10</sub>, terminal H (24, 0.720–0.734, contested); NH<sub>3</sub>BH<sub>3</sub> (12, 0.773–0.785, one-centre WFs of H) | 36 | 0 | 36 | 1.684 |
| all | | | 137 | | 129 | 101 | 36 | |

- *Z*<sub>eff,H</sub> = 0.6875. For every contested X–H WF *n*<sub>inner</sub> = 0
  at both ends, so *Z*<sub>eff</sub> = ZVAL − 5/16 (ZVAL: O 6, S 6, N 5, C 4,
  B 3, H 1).
- Force ratio: (*Z*<sub>eff</sub>/*d*<sup>2</sup> of the end that takes the
  pair) / (*Z*<sub>eff</sub>/*d*<sup>2</sup> of the other end).
- 125 of the 137 WFs are contested; for each of them *f* < *f*\* is equivalent
  to the plugin's verdict "pair to X" (asserted). 12 are one-centre WFs of H.
- H atoms of clean38: 153; reference −1 for 52, +1 for 101; the plugin gives
  the reference value to 153/153. A pure nearest-atom rule gives −1 to 145.
- The nine cases are fully correct in `results/clean38_r3-T-B.tsv`:
  B<sub>6</sub>H<sub>10</sub> 40/40, NH<sub>3</sub>BH<sub>3</sub> 32/32,
  CH<sub>2</sub>F<sub>2</sub> 5/5, HNO<sub>3</sub> 40/40,
  (CH<sub>3</sub>)<sub>3</sub>NO 14/14, H<sub>2</sub>S<sub>2</sub> 4/4,
  H<sub>2</sub>S<sub>3</sub> 5/5, H<sub>2</sub>S<sub>5</sub> 7/7,
  SrBe(OH)<sub>4</sub> 160/160.

Panel (c), examples (the figure shows two or three decimals).

| Example | End | ZVAL | *n*<sub>inner</sub> | *Z*<sub>eff</sub> | *d* (Å) | *Z*<sub>eff</sub>/*d*<sup>2</sup> | *ρ*<sup>2</sup> (Å<sup>2</sup>) | *χ* | OS |
|---|---|---|---|---|---|---|---|---|---|
| As–S pair, As<sub>4</sub>S<sub>4</sub> (WF 63), decided by *χ* | S | 6 | 0 | 5.6875 | 0.9491 | 6.3141 | 1.7420 | 3.2649 | −2 |
| | As | 5 | 0 | 4.6875 | 1.3248 | 2.6710 | 1.9571 | 2.3952 | +2 |
| S–H pair, H<sub>2</sub>S<sub>2</sub> (WF 3), decided by *Z*<sub>eff</sub>/*d*<sup>2</sup> | S | 6 | 0 | 5.6875 | 0.8927 | 7.1373 | 1.3340 | – | −1 |
| | H | 1 | 0 | 0.6875 | 0.4635 | 3.2001 | – | – | +1 |

- Ratio of the two ends: *χ* 1.3631 for the As–S pair (the smallest of the 32
  As–S pairs of the case); *Z*<sub>eff</sub>/*d*<sup>2</sup> 2.2304 for the S–H
  pair (the smaller of the two S–H pairs).
- The four constants of the method: 5/16, REACH_SIGMA = 1, PARITY = 1/8, 10<sup>−6</sup>.

## Provenance

- **Method**: every class, verdict, *Z*<sub>eff</sub>, *d* and *ρ*<sup>2</sup>
  comes from the plugin itself, run on the case as the clean evaluator shows it
  (`pyleos.bench.load_plugin`, `pyleos.bench.scrub`): `_analyse(case)` (WF
  classes, hosts and the decision records `notes`), `_Shells.seen_from(d)`
  (*n*<sub>inner</sub>, *ρ*<sup>2</sup>) and `assign(case)` (oxidation states).
  The script adds only *f* and *f*\*.
- **Cases**: the 38 cases of `pyleos/sets/clean_set.json` outside
  `splits.test_regression`, directories resolved by
  `pyleos.bench.collect_clean` through `cases/INDEX.tsv`. References:
  `pyleos/sets/reference_os.json`, per-atom labels by
  `pyleos.bench.site_labels`.
- **Panel (a)**: `cases/h2s2/v3/POSCAR` and `EOSCAR` (#WF block, 7 centres).
- **Panel (b)**: `cases/{b6h10,bh6n,ch2f2}/v1` (centres from `wannier90.wout`)
  and `cases/{hno3,tmao,h2s2,h2s3,h2s5,srbe_oh_4}/v3` (centres from the EOSCAR
  #WF block). All are non-spin-polarised.
- **Panel (c)**: `cases/as4s4/v3` and `cases/h2s2/v3`.

## Selection rules (constants at the top of the script)

1. *Terminal X–H bond*, as in Fig. 1(b): X is the non-H atom nearest to the H
   atom (all periodic images count); the bond is used if the next non-H
   neighbour is farther than 1.25 × *d*(X–H). Accepted bonds have ratios
   ≥ 1.57 and lengths 0.970–1.356 Å. The 16 bridging H atoms of
   B<sub>6</sub>H<sub>10</sub> (three-centre B–H–B) are not terminal; their WF
   is a one-centre WF of H and they are −1, as in the reference.
2. *WF of that bond*: the WF centre nearest to the H atom, accepted if it lies
   within 0.10 Å of the X–H axis and 0 < *f* < 1. All terminal bonds pass; the
   largest offset is 0.030 Å.
3. *f* = (**r**<sub>WF</sub> − **r**<sub>X</sub>)·(**r**<sub>H</sub> − **r**<sub>X</sub>) / |**r**<sub>H</sub> − **r**<sub>X</sub>|<sup>2</sup>, minimum-image convention.
   On the bond axis a WF is contested for 0.25 < *f* < 0.75.
4. *Rows of panel (b)* are ordered by *f*\* (largest first), then by bond length.
5. *Examples of panel (c)*: in each case the contested WF with the smallest
   ratio between its two ends.
6. *View of panel (a)*: the projection of Fig. 1(a) (same routine): all drawn
   objects kept apart, then the least foreshortened bonds (the
   shortest-looking bond is drawn at 0.79 of its length). Atom circles are
   symbols, not radii. Atom labels number the atoms of each element in POSCAR
   order.

## Checks made by the script (it stops if one fails)

- Each decision record of the plugin is matched to its WF and its two atoms,
  and its *Z*<sub>eff</sub> equals ZVAL − *n*<sub>inner</sub> − 5/16 with
  *n*<sub>inner</sub> from the plugin's shells.
- For all 125 contested X–H WFs: distances recomputed with numpy from POSCAR
  and the WF centres agree with the plugin's to 10<sup>−9</sup> Å; the ratio of
  the two forces equals the plugin's margin; *f* < *f*\* agrees with the verdict.
- For the two examples: *d* and *ρ*<sup>2</sup> recomputed with numpy from the
  WF centres and spreads agree with the plugin's to 10<sup>−9</sup>; the ratio
  of the two *χ* equals the plugin's margin.
- Per-case numbers of correct atoms, recomputed from `assign(case)` and the
  per-atom reference, equal those of `results/clean38_r3-T-B.tsv` for all 38
  cases.
