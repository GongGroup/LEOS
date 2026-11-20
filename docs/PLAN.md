# LEOS project plan

This document states the problem LEOS solves, the method, the main design decisions
and the engineering conventions. Test and application results are in
[RESULTS.md](RESULTS.md); terms are defined in [CONTEXT.md](../CONTEXT.md).

## 1. Problem

An oxidation state is an integer that keeps the electron account of an atom: assign
all electrons of every chemical bond to the more electronegative end, then count how
many electrons each atom has gained or lost relative to the neutral atom. In CO2, O
is more electronegative than C, so the electrons of the C–O bonds all count for O:
C is +4 and O is −2. The oxidation state is a convention, not a measurable physical
quantity, but it is the common language for redox and bonding.

The electron density from DFT is continuous and belongs to no particular atom. To
get the oxidation state of every atom automatically from a DFT result, one must
decide which atom each portion of the electrons counts for.

LEOS (Lattice Effective Oxidation States) does this in two steps:

1. **Wannierisation**: inside the VASP calculation, transform the occupied states
   into a set of localised Wannier functions (WFs). Each WF is a localised packet of
   electrons with centre coordinates and a spread (spatial extent). In a
   closed-shell calculation each WF holds 2 electrons; in a spin-polarised
   calculation each WF of each spin channel holds 1.
2. **Assignment**: decide which atom receives the electrons of each WF, then compute
   the oxidation states from the electron counts. Intuitively, a WF whose centre
   sits on an atom is usually a lone pair of that atom; a WF between two atoms is a
   covalent bond and goes to one end according to electronegativity.

The code has two parts:

- **VASP patch** (patch/leos_v3.patch, Fortran): calls the wannier90 library inside
  VASP 6.4.3 to wannierise (SCDM initial guess), writes the centre, spread,
  occupation and projection weights of every WF to the EOSCAR file, and does the
  assignment inside VASP.
- **pyleos** (Python standard library, no third-party dependency): reads EOSCAR,
  redoes the assignment, and gives diagnostics and the reliability index.

Both parts use the same assignment algorithm and their results agree atom by atom.
The algorithm is settled in pyleos first, then ported to Fortran.

A separate tool, gpuw90 (PyTorch), reimplements the disentanglement and localisation
of wannier90 on the GPU to speed up the wannierisation; see
[gpuw90/README.md](../gpuw90/README.md).

## 2. Method

### Assignment as a minimum-cost flow

The most direct assignment is greedy: each WF independently picks the nearest atom
(or also weighs electronegativity). Its flaw is that the WFs know nothing of each
other: when many WFs are all "nearest" to the same atom, they all pile onto it and
the result is physically impossible (see the As4S4 example under "Key decisions").

LEOS assigns all WFs together, as a minimum-cost flow:

- each WF has a fixed number of electrons to hand out (2 closed-shell, 1
  spin-polarised);
- each atom has a lower and an upper limit on the electrons it can receive;
- giving the electrons of a WF to an atom has a cost;
- the assignment with the lowest total cost is sought.

Electron counts and capacities are integers, so the optimum of the minimum-cost flow
is automatically an integer assignment and the oxidation states need no rounding.
The implementation is in pyleos/transport.py.

### Cost

The default cost uses only geometry and electronegativity:

    cost(WF, atom) = d² + β · (1 − χ / χmax)

d is the distance from the WF centre to the atom (Å), χ the electronegativity of the
atom, χmax the largest electronegativity in the system; β is 10 by default. The
default electronegativity scale is Allen (pyleos/capacity.py); `--eneg pauling`
switches to Pauling. For the choice of β see "Parameter choices" in RESULTS.md.

### Capacity

The number of electrons an atom receives is limited to [n_min, n_max], set by the
physical limits of the valence shell (OS_BOUNDS in pyleos/capacity.py). For example,
O holds at most 8 valence electrons, so its oxidation state cannot go below −2.
Semicore electrons (such as Na 2s2p) must be counted in full; n_min ensures this.

### Cap on the electron total

The total number of electrons assigned does not exceed NELECT, the true electron
count of the DFT calculation. This handles two accounting problems at once:

- **ghost**: when there are more WFs than occupied electron pairs (metals or
  fractional occupations), the surplus WFs have no electrons to hold. The cap makes
  the solver drop the most expensive part, so no electrons are over-counted.
- **electron loss**: when some electrons did not enter the Wannier subspace, the
  same cap makes the assignment use the number of electrons actually present.

### Projection-assisted arbitration (optional)

The #WP block of EOSCAR records the weight of each WF on the PAW projectors of each
atom. Projection-assisted arbitration reassigns a WF according to the projection
only when the costs of its two best candidate atoms differ by less than TAU and the
projection weight leads by more than PMARGIN. (TAU, PMARGIN) = (1.0, 0.20), defined
as ARBITER_TAU and ARBITER_PMARGIN in pyleos/transport.py. It is off by default; the
parameter scan and the range of validity are under "Parameter choices" in
RESULTS.md.

Use: the library call `assign_case(case, tier2="arbiter")`, or
`python3 -m pyleos bench full --tier2 arbiter`. **`--gamma` of cli.py uses a
different mode, a linear mix (`tier2="linear"`), not the projection-assisted
arbitration described here.**

### Fragments

`--fragment 2,3,4` merges a group of atoms into one unit that receives electrons
(capacities add; the cost is that of the member nearest to the WF) and reports the
oxidation state of the unit. Use it for delocalised anions (such as I3⁻) and metal
clusters (such as Au12, Co4), where the picture "each WF belongs to one or two
atoms" does not hold. When the two candidate atoms of an electron pair have almost
the same cost, pyleos suggests a fragment.

### Reliability index R

After the optimal assignment is found, ask for each atom: if it is forced to take 1
electron more or fewer, how much extra cost m does a global re-solve pay? This is a
full re-solve; the other WFs may rearrange. m is normalised to 0–100 by
R = 100·m/(m+1); the higher R, the more stable the current assignment.

- **R_site**: the minimum over all atoms. Exchanges of electrons across a
  homonuclear bond or between equivalent sites count too, so with site degeneracy
  R_site is close to 0.
- **R_elem**: counts only perturbations that change the total electron count of an
  element; exchanges between atoms of the same element do not count. It corresponds
  to the element-mean oxidation state that dev checks.

R is the counterpart of the R(%) index of Comas-Vilà 2025. The lower limit for
semicore electrons is a hard constraint and takes no part in the perturbations of R.

## 3. Reading the output

### pyleos output

```bash
python3 pyleos/cli.py case cases/as4s4/v3 -v
```

```
=== as4s4  (ISPIN=1, WFs=88, electrons=176) ===
greedy header: NELECT(Wannier)=176.0  |dNe|=0.0000
reliability R_elem = 69/100 (element ledger)   R_site = 2/100 (...)
note: 8 electron pair(s) sit at a near-tie between hosts (...) - ... consider --fragment
 ion  el   ne  OS_flow  OS_greedy     R
   1   S    8       -2         -4    69
  ...
  17  As    4       +1         +4    17
  el     flow min..max (mean)     greedy min..max (mean)
  As           +1..+3 (+2.00)             +4..+4 (+4.00)
   S           -2..-2 (-2.00)             -6..-2 (-4.00)
```

- `ne`: valence electrons assigned to the atom; `OS_flow`: result of the flow;
  `OS_greedy`: greedy result; `R`: reliability index of the atom.
- The last table gives the minimum, maximum and mean of each element.
- The `note` line flags near-tied electron pairs, the places where a fragment may be
  considered. Here they are 8 As–As bonds.

### The two columns of EOSCAR

In the per-atom table of EOSCAR, the `OS(int)` column is the greedy result and the
`OS(flow)` column is the result of the flow (present from patch v3.1 on). When
reading EOSCAR directly, use `OS(flow)`; `OS(int)` can show values outside the
physical range. Analyse an old EOSCAR that has no `OS(flow)` column with pyleos.

### Quality warnings

| Warning | Meaning | Criterion | What to do |
|---|---|---|---|
| stray | The WF is not well localised; its centre carries no chemical information | spread > 3 Å², or centre more than 1.6 Å from every atom | Element means are usually still usable, single atoms are not; redo the wannierisation (more iterations or a better initial guess) |
| ghost | More WFs than occupied electron pairs | The part of the total WF capacity that exceeds NELECT | Handled automatically by the cap on the electron total; but with ghost > 0 the oxidation state of electropositive metals may be overestimated by 1–2 (see "Systematic bias of ghosts on metal oxidation states" in RESULTS.md) |
| dNe | The electron count of the Wannier representation differs from DFT | Below 1e-6 in a healthy calculation | Some electrons did not enter the Wannier subspace; the diagnostic block of EOSCAR shows in which spin, k-point and band they are lost |

### Using R

- **Low R means the solution is not unique, not that it is wrong.** The typical case
  is degeneracy between equivalent atoms, such as the As–As bonds of As4S4 (R_site
  1.8) or the Co atoms inside the Co4 cluster (4–22). Then read the element mean or
  the total oxidation state of the fragment.
- **High R does not guarantee a correct result.** When the input itself is wrong
  (wrong structure, failed wannierisation), the solution can be both stable and
  wrong. Judge reliability from R together with the three quality warnings.

## 4. Key decisions

### Global assignment replaces greedy

The implementation at the start of the project (it produced the data in cases/*/v1)
used greedy assignment. A review found three problems:

1. **Greedy can give impossible results.** In As4S4 one S received 14 electrons,
   while the valence shell of S holds at most 8. The old code used EOS_O_MAX to keep
   O from exceeding 8 electrons, but this protected only O and no other element.
   Adding limits element by element never ends.
2. **The other assignment branch (Hungarian algorithm) was unusable.** It gave each
   atom as many slots as the neutral atom has electron pairs, so an anion could
   never receive more electrons than the neutral atom (O²⁻ needs 4 pairs; there
   were only 3 slots). None of the inputs of the 198 legacy calculations ever
   enabled it; it has been deleted.
3. **Advertised features were not wired in.** Frame-by-frame analysis of molecular
   dynamics (EOS_EVERY_N) was hard-coded in the main program to run only at the
   first step; the code comments claimed to follow the IUPAC 2014 rules but used
   Pauling electronegativities, whereas the IUPAC report goes with the Allen scale.

The root of all three is that each WF decides on its own, so all WFs are now
assigned together by a minimum-cost flow. The capacity upper limits solve the excess
problem for all elements at once.

### The main contribution is the minimum-cost-flow formulation

Assigning oxidation states from Wannier centres has precedents: Jiang, Levchenko
and Rappe 2012 proposed it first; Gimferrer and Salvador 2020 analysed the failures
of the "nearest atom" criterion systematically but gave no fix; Comas-Vilà 2025 made
a version for solids with QTAIM atomic overlap matrices. No precedent was found for
writing the electron assignment as a global minimum-cost flow with capacity
constraints. The main contribution of the paper is therefore this formulation and
its correction of the known failure modes.

### Two implementation paths

The VASP patch does the assignment internally (for frame-by-frame analysis of
molecular dynamics); pyleos does offline analysis, diagnostics and statistics. The
assignment problem has only tens to hundreds of variables, and Python solves it in
milliseconds; development and testing are much faster in Python than in Fortran.
Hence pyleos first, then the patch.

### The default cost uses only geometry and electronegativity

The cost can be built at three tiers: geometry plus electronegativity; the same plus
projection weights; the full atomic overlap matrix. The third tier is the method of
Comas-Vilà; it is expensive and gives up the cost advantage of LEOS, so it is not
adopted. Projection weights serve as an optional arbitration and are not in the
default; all default results depend only on geometry and electronegativity, with no
adjustable empirical parameter (except β).

### Capacities follow physical limits

The capacity limits are set by the physical limits of the valence shell, not by "the
chemically common range of oxidation states". Constraining with the common range
would constrain the answer with the answer, and would exclude results such as
high-pressure phases and rare oxidation states. Physical limits exclude only the
impossible.

### Allen electronegativity scale

The IUPAC 2014/2016 technical report on oxidation state goes with the Allen scale.

### β = 10

Results do not change for β in 8–12; the middle value 10 is taken. With β too small
(such as 3), X–H bonds are assigned wrongly without any warning. For the scan data
see "Parameter choices" in RESULTS.md.

### Projection-assisted arbitration uses a threshold rule

The projection may change an assignment only when the geometric decision is nearly
tied; assignments that geometry makes clear are never changed. Mixing the projection
weights in linearly is sensitive to the parameters, and the projection itself is
biased towards heavy atoms. For details see "Parameter choices" in RESULTS.md.

### Warm start is kept as an option

Adjacent frames of a trajectory have similar structures, so the Wannier centres of
the previous frame can be written into the wannier90 input of the next frame as
guiding centres. The feature is implemented but gives no speed-up in the workflow
with the SCDM initial guess ("Efficiency" in RESULTS.md); it is kept for workflows
without the SCDM initial guess.

### Choice of application systems

- **Li2O2 delithiation**: on charging, O changes oxidation state, not the metal. The
  DFT+U+V occupation-matrix method can only track d/f metal sites and cannot handle
  this anionic redox.
- **CeO2 oxygen vacancy**: the 2 electrons left by the vacancy turn two Ce⁴⁺ into
  Ce³⁺; this also tests the wannierisation of 4f electrons.
- **Li-rich cathodes + machine-learning potential**: Li2MnO3 (ordered, 3d), Li2RuO3
  (ordered, 4d), NMC (SQS disordered). The pretrained MACE model is used only to
  sample configurations; all oxidation states are assigned from DFT single points.

## 5. Paper and literature

The paper has two main lines:

1. the minimum-cost-flow formulation of oxidation-state assignment (global optimum,
   integer solution, capacity constraints);
2. a fast, large-scale Wannier workflow (in situ assignment, SCDM-only screening,
   offline wannier90, gpuw90).

Greedy serves only as the baseline; the main comparison is with Comas-Vilà 2025.

| Reference | Content | Relation to LEOS |
|---|---|---|
| Comas-Vilà et al., JCTC 2025 | QTAIM atomic domains + atomic overlap matrices, four programs in a chain; the R(%) reliability index; a test set of 40 solids; the origin of the name "EOS" | LEOS does no spatial integration, finishes within one VASP run and can do trajectories (the paper reports no timings, so speed cannot be compared directly); R and fragments both have counterparts; all 40 solids have been computed. See [comas_vila_2025_notes.md](comas_vila_2025_notes.md) |
| Gimferrer, Comas-Vilà, Salvador, Molecules 2020 | Shows systematically that the criterion "centre of a localised orbital + nearest atom" assigns O–H and N–H bonds as hydrides | The corrections of cases such as SrBe(OH)4 address this kind of failure |
| Marzari group, npj Comput. Mater. 2025 | Tracks oxidation states along molecular dynamics with DFT+U+V occupation matrices | Overlaps with the trajectory feature of LEOS, but handles only Hubbard sites (d/f metals), not O |

Reference PDFs are in docs/refs/; they are copyrighted and not in git.

**A note on terminology.** To condensed-matter readers "transport" first suggests
electron transport (conductivity, mobility), which has nothing to do with the
transportation problem of operations research. At first mention in the paper write
"minimum-cost-flow / optimal-assignment formulation", or state "optimal transport
(in the operations-research sense)".

## 6. To do

To-do items are kept only as GitHub issues: https://github.com/GongGroup/LEOS/issues

## 7. Engineering conventions

- **Data layout.** The cluster keeps the full run directories (including large files
  such as wavefunctions). The local copy holds only small text files: calculations
  are stored by system in `cases/<system>/<run>/` (EOSCAR, input files, wout), and
  every run has one row in cases/INDEX.tsv; summary tables lie flat in results/. No
  directories per batch. .gitignore blocks large files and pseudopotentials.
- **Cluster location.** Calculations run on the SAI cluster (ssh alias sai). Since
  the home-directory migration of 2026-07-25 the project directory is
  `/org/sjtu-gongxueqing/zhangxianyang/temp/org/LEOS/` (its README.txt describes
  the directory), and the patched VASP is in
  `/org/sjtu-gongxueqing/zhangxianyang/temp/org/apps/vasp.6.4.3_leos_v31_20260711/`.
  Old job scripts, logs and symlinks under that directory still point to the
  pre-migration `/home/sjtu-gongxueqing/zhangxianyang/`; replace the paths before
  reusing them. In old documents `~/<name>` means `<name>` under the LEOS directory.
- **Current state only.** Documents state only the currently correct results, not
  how they evolved. Plans and review reports of finished stages are not kept in the
  working tree; use git log / git show when needed.
- **To-do items only in issues.** Documents contain no to-do items.
- **Numbers have a source.** Every result number in a document must map to a file
  under results/ or cases/, with the command that produced it or how it is counted
  (for example, whether a timing includes the computation of R). To change a number,
  first change or regenerate the data file, then the document.
- **Golden table.** results/bench74_golden.tsv is what CI compares case by case. A
  change that alters assignment results must first regenerate this table and explain
  the differences case by case, then update the numbers in the documents.
- **Snapshot documents state their status.** A document that records the conclusions
  of one moment, such as a health-check report, states its current status at the top.
