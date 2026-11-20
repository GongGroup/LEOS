# LEOS glossary

The documents of this repository use the terms below. Method: docs/PLAN.md.
Results: docs/RESULTS.md.

## Basic concepts

**Oxidation state (OS)**
The number of electrons an atom has lost or gained relative to the neutral
atom once every bond pair is given whole to the more electronegative end. An
integer. It is a bookkeeping convention, not a measurable quantity. Do not mix
it with valence (bonding capacity) or charge (the continuous values of Bader
and similar methods).

**Wannier function (WF)**
A localised packet of electrons obtained by transforming the occupied DFT
states; it has a centre and a spread. In a closed-shell calculation
(ISPIN = 1) each WF holds 2 electrons; in a spin-polarised one (ISPIN = 2)
each WF of each spin channel holds 1.

**Wannier centre**
The centre of a WF (Cartesian, Å), the geometric basis of the assignment. A
centre close to one atom is usually a lone pair; a centre between two atoms is
a covalent bond.

**Spread**
The spatial extent of a WF (Å²), a measure of the quality of the Wannier step.
A well-localised WF is generally below 1–2 Å².

**EOSCAR**
The file written by the LEOS VASP patch and read by pyleos. It holds the
per-atom table (greedy column OS(int), flow column OS(flow)), the #WF block
(centre, spread and occupation of every WF), the #WP block (weight of every WF
on the PAW projectors of each atom) and the electron-loss diagnostics.

## Assignment

**Greedy**
The baseline: each WF goes, on its own, to the nearest atom (with an
electronegativity correction), and the WFs know nothing of each other. It can
give physically impossible answers, such as one S receiving 14 electrons.
Column OS(int) of EOSCAR.

**Minimum-cost flow ("the flow")**
The LEOS assignment: the electrons of all WFs are distributed over the atoms
together, each atom has a lower and an upper limit on what it may receive, and
the plan of least total cost is taken. Electron counts and capacities are
integers, so the optimum is integral. The code calls it transport
(`transport.py`); in the paper the first mention reads "minimum-cost-flow /
optimal-assignment formulation", so that it is not read as electron transport.

**OS(flow)**
The last column of the per-atom table of EOSCAR: the flow assignment computed
inside the VASP patch (INCAR tag EOS_FLOW, on by default). It equals the
pyleos default atom by atom.

**Cost**
The cost of giving the electrons of a WF to an atom: d² + β·(1 − χ/χmax). d is
the distance from the WF centre to the atom, χ the electronegativity (Allen
scale by default), β = 10 by default.

**β**
The weight of the electronegativity term, 10 by default. Results are stable
for 8–12; below that X–H bonds are assigned wrongly, with no warning.

**Capacity**
The lower and upper limits [n_min, n_max] on the electrons an atom receives,
set by the physical limits of the valence shell. They exclude the impossible,
not the rare; they are not set from "chemically common oxidation states".

**Semicore**
An inner shell that must be counted full (2s2p of Na, for example): the part
of n_min that has to be met.

**Projection-assisted arbitration**
An optional rule: when the costs of the two best candidate atoms of a WF
differ by less than TAU (1.0) and one of them leads in #WP projector weight by
more than PMARGIN (0.20), the WF follows the projection. Off by default. Not
suited to organic systems with many C–H bonds. `tier2="arbiter"` in the code;
"tier 2" and "arbiter" in older text mean this, and "tier 1" means the default
cost.

**Fragment**
A user-defined group of atoms that receives electrons as one unit (capacities
add; the cost uses the member nearest to the WF) and reports one oxidation
state. For delocalised anions and metal clusters (I3⁻, Au12, Co4 and the
like).

**Site degeneracy**
Equivalent atoms can exchange electrons at almost no change in total cost: the
two ends of a homonuclear bond, equivalent sites of a cluster, equivalent O2
units. The element mean is then reliable and the value of a single atom is not
unique; R_site is close to 0. Read the element mean or the fragment total.

**Reliability index R (R_site, R_elem)**
After the optimum is found, one atom is forced to take one electron more or
less and the problem is solved again; m is the extra cost. R = 100·m/(m+1),
from 0 to 100. R_site counts every perturbation, so degeneracies such as the
ends of a homonuclear bond bring it close to 0; R_elem counts only the
perturbations that change the electron total of an element. R is the
counterpart of the R(%) of Comas-Vilà. A low R means the solution is not
unique; a high R does not guarantee that it is right.

## Quality warnings

**stray**
An unlocalised WF: spread > 3 Å², or centre more than 1.6 Å from every atom.
The Wannier step itself has failed; the assignment cannot repair it. Redo the
Wannier step.

**ghost**
When there are more WFs than occupied electron pairs, the extra WFs that have
no electrons to hold (common in metals and with fractional occupations). The
cap on the electron total (NELECT) makes the solver drop the costliest part.
With ghost > 0 the oxidation state of electropositive metals may be
overestimated.

**dNe**
The electron total of the Wannier representation minus the true DFT total.
Below 1e-6 in a healthy calculation; a clear non-zero value means that
electrons did not enter the Wannier subspace (electron loss).

**Diagnostic EOSCAR**
The EOSCAR the patch writes before it aborts when the electron loss exceeds
EOS_NELECT_ABORT. It carries the mark `# LEOS_ABORTED` and blocks that locate
the loss (#NL for the window type, #NLD for the diffuse type). It is not
written when wannier90 itself fails.

## Benchmarks and data

**74-case benchmark**
Three batches of test systems: the legacy regression set (36), extension
wave 1 (13) and extension wave 2 (25). 72 of them have reference oxidation
states.

**Reference oxidation states**
Oxidation states labelled by hand from chemical knowledge and the IUPAC rules,
in pyleos/sets/reference_os.json. Mixed-valence systems carry the element mean
(Fe in Fe3O4: +2.667).

**dev**
The largest absolute difference between the mean oxidation state of an element
and its reference. dev ≤ 0.25 is a pass. Only element means are checked, not
single atoms.

**Golden table**
results/bench74_golden.tsv. CI compares a fresh benchmark run with it line by
line (columns status and dev). Its counterpart in gpuw90 is the wannier90
reference output (.wout).

**Comas-Vilà**
The reference method: Comas-Vilà et al., JCTC 2025, effective oxidation states
from QTAIM atomic domains and atomic overlap matrices, tested on 40 solids.
The name "EOS" comes from that paper.

**Case / run**
A case is one chemical system, `cases/<case>/`; a run is one calculation of it,
`cases/<case>/<run>/`. The run name gives the program version and what differs:
v1 is the early program of the project (greedy), v3 the current patch, v31-*
a recomputation or diagnostic with v3.1 (v31-wp carries #WP, v31-ghost is the
ghost A/B test); structure/ and structure-fixed/ hold structures only. The
purpose (role), benchmark set (bench_set) and former path (origin) of each run
are in cases/INDEX.tsv.

**Benchmark set**
The 74-case benchmark has four groups: v1_mol (10 molecular crystals),
v1_oxide (26 oxides and halides), wave1 (13) and wave2 (25). They are the
bench_set column of INDEX.tsv and the batch column of the golden table.
wave1 and wave2 are also called extension wave 1 and 2 (formerly M4). The
application systems were formerly M5: M5-A is Li2O2, M5-B the three Li-rich
cathode systems, M5-C CeO2.

## Wannier workflow

**In situ / offline**
In situ: Wannier step and assignment inside VASP (the greedy and OS(flow)
columns), usable frame by frame in molecular dynamics. Offline: outside VASP,
with pyleos, gpuw90 or a stand-alone wannier90.x.

**Localisation (wannierise)**
The Marzari–Vanderbilt method: a unitary transformation within a given
subspace that minimises the total spread and yields the Wannier centres and
spreads.

**Disentanglement**
When there are more bands than WFs (num_bands > num_wann; metals, or
conduction bands included), the subspace to be localised is selected first
(Souza–Marzari–Vanderbilt). Localising without this step gives wrong results
and no error.

**Frozen window**
The energy range that disentanglement must keep unchanged. It usually covers
the occupied states, so that no valence electron is lost in the choice of
subspace.

**Automatic frozen window (AUTO_WINDOW)**
The patch sets the frozen window from the occupations (INCAR tag
LWANNIER90_AUTO_WINDOW); no hand-written dis_froz_max is needed. A missing
window in hand-written input is the main cause of the NELECT mismatch failure.

**SCDM**
The initial guess of the Wannier step (selected columns of the density
matrix); it needs no hand-set projections. SCDM-only is the fast mode that
skips the iterative localisation and uses the centres of the initial guess.

**Guiding centres**
How the warm start is implemented: the Wannier centres of one frame are
written into the wannier90 input of the next.

**gpuw90**
The disentanglement and localisation tool of this repository, written in
PyTorch. The same code runs on CPU and CUDA; wannier90 and VASP are not
modified.

**Failure triage (FAIL_DFT / FAIL_W90 / FAIL_LEOS)**
How batch.py classifies failed cases. FAIL_DFT: the SCF crashed or OUTCAR is
incomplete. FAIL_W90: DFT is fine but wannier90 did not finish. FAIL_LEOS: the
Wannier step is fine but there is no usable EOSCAR (missing, or marked
LEOS_ABORTED). Each class has its own automatic retry settings.
