# pyleos

pyleos reads the EOSCAR written by the LEOS VASP patch, assigns the electrons
of every Wannier function (WF) to atoms, and reports integer oxidation states,
quality warnings and a reliability index. The assignment is a
capacity-constrained minimum-cost flow: globally optimal and integral. Python
standard library only, no third-party dependency.

Method: [docs/PLAN.md](../docs/PLAN.md); how to read the output: its section
"Reading the output".

## Usage

```bash
# One case. The directory needs POSCAR, EOSCAR and wannier90.wout
# (wannier90.1.wout and wannier90.2.wout when spin-polarised)
python3 -m pyleos case cases/as4s4/v3 -v

# Legacy regression set: flow against greedy on the 36 v1 runs, scored
# against the reference
python3 -m pyleos regress

# Benchmarks (see "Benchmarks" below)
python3 -m pyleos bench full --no-reliability --out /tmp/bench.tsv
```

`python3 pyleos/cli.py ...` is the same as `python3 -m pyleos ...`.

Options of `case`:

| Option | Default | Effect |
|---|---|---|
| `--beta B` | 10.0 | Weight of the electronegativity term. Below about 8, X–H bonds are assigned wrongly with no warning; see RESULTS.md, "Parameter choices" |
| `--eneg allen\|pauling` | allen | Electronegativity scale |
| `-v` | | Print the per-atom result |
| `--fragment 2,3,4` | | Merge these atoms (numbered from 1) into one unit that receives electrons; repeatable |
| `--drop-stray` | | Exclude unlocalised WFs (stray) before assigning, to tell whether an error comes from the input or from the assignment. The electron total is then smaller, and the element means can no longer be compared with the reference |
| `--gamma G` | 0.0 | Linear mixing of the projector weights: the cost is lowered by G times the normalised weight. Sensitive to the parameter and biased towards heavy atoms; not recommended. The recommended threshold rule is under "Library use" |

Options of `regress`: `--beta`, `--eneg`, `--reliability` (also compute R).

## Benchmarks

`python3 -m pyleos bench` has three subcommands. All of them find the case
directories through cases/INDEX.tsv.

| Subcommand | What it runs |
|---|---|
| `full --out FILE` | The 74-case benchmark, in the format of results/bench74_golden.tsv (`--out` is required, so that the golden table is not overwritten by accident). `--greedy` scores the greedy column instead; `--no-reliability` skips R |
| `clean` | The clean set: element-mean dev and per-atom accuracy. `--method` selects the assignment plugin (default methods/baseline.py) |
| `test --set regression\|homonuclear\|both --method PLUGIN` | The pre-registered held-out test sets. `--selftest` only checks the scoring path on the dev set |

An assignment plugin is a .py file in methods/ that defines `assign(case)`.
methods/frozen/ holds the frozen candidates, which are never edited again
(their sha256 are in sets/test_set.json). A plugin sees neither the results of
other methods nor the file system; the rules are at the top of bench.py.

The benchmarks are defined in sets/: reference_os.json (reference oxidation
states of the 74 cases, grouped as dataset, oxides, wave1, wave2),
clean_set.json (per-atom labels) and test_set.json (held-out test sets).

## Library use

```python
from pathlib import Path
from pyleos.parsers import load_case
from pyleos.transport import assign_case, reliability_margins

case = load_case(Path("cases/as4s4/v3"))
res = assign_case(case)                      # default: beta=10, Allen, geometry + electronegativity
res = assign_case(case, tier2="arbiter")     # projection-assisted arbitration; needs the #WP block of EOSCAR
res = reliability_margins(case, res=res)     # compute R_site / R_elem
```

The thresholds of the projection-assisted arbitration are ARBITER_TAU (1.0)
and ARBITER_PMARGIN (0.20) in transport.py.

## Modules

| File | Content |
|---|---|
| parsers.py | Reads POSCAR, EOSCAR (with the #WF and #WP blocks) and the wannier90 .wout into a CaseData |
| capacity.py | Capacity limits OS_BOUNDS of each element; Allen and Pauling electronegativity tables |
| transport.py | Minimum-cost-flow assignment, stray and ghost diagnostics, reliability index R |
| cases.py | Reads cases/INDEX.tsv (the case registry) |
| bench.py | The three benchmarks: 74 cases, clean set, held-out test sets |
| cli.py | Command line (`python3 -m pyleos`) |
| batch.py | Cluster batch driver (below) |
| xsd2poscar.py | Converts a Materials Studio .xsd to POSCAR |

One compatibility rule in the EOSCAR parser: the VRHFIN field of the official
Xe pseudopotential reads "X", which truncates the element name in EOSCAR;
parsers.py matches element names by prefix.

## Model

- Graph: source → WF nodes → atom nodes → sink. Closed-shell (ISPIN = 1): each
  WF supplies 2 electrons. Spin-polarised (ISPIN = 2): each WF of each channel
  supplies 1, and both channels share the atom capacities in one graph.
- Cost: `d²(WF centre, atom) + β·(1 − χ/χmax)`, χ the Allen electronegativity
  by default.
- Atom capacities [n_min, n_max] come from OS_BOUNDS in capacity.py and follow
  the physical limits of the valence shell. Semicore shells must be counted
  full; n_min is enforced with negative-cost edges.
- The total flow does not exceed NELECT. Of the surplus WFs (ghost) the solver
  drops the costliest part.
- Electron counts and capacities are integers, so the optimal flow is integral
  and the oxidation states need no rounding.
- Stray criterion: spread > 3 Å², or centre more than 1.6 Å from every atom.

Three points of the solver need care:

1. The shortest-path potential update must cover the nodes that were not
   reached; otherwise reduced costs turn negative and the algorithm becomes
   exponential.
2. The large reward and the floating-point tolerance must be far apart in
   magnitude (1e3 and 1e-9); otherwise rounding residue counts as an
   improvement and the loop never ends.
3. The number of iterations is capped, and an anomaly raises a clear error
   instead of looping for ever.

## batch.py: cluster batch runs

batch.py depends on the standard library only and can be copied to a cluster
as a single file; inside the repository it is also
`python3 -m pyleos batch ...`.

```bash
python3 batch.py stage   --dest ~/leos_batch1 --cases cases.json
python3 batch.py submit  --dest ~/leos_batch1
python3 batch.py status  --dest ~/leos_batch1
python3 batch.py retry   --dest ~/leos_batch1
python3 batch.py harvest --dest ~/leos_batch1 --out ~/leos_batch1_results
```

Format of cases.json: `[{"name": "co2", "src": "<directory with POSCAR, POTCAR,
KPOINTS>", "incar": {"ENCUT": 500}}, ...]`. Each case directory receives the
LEOS INCAR template (merged with the settings of the case) and a 4-GPU sbatch
script. The state is kept in `<dest>/batch_state.json`.

Failure triage (printed by status and harvest):

| Status | Meaning | What retry does |
|---|---|---|
| FAIL_DFT | The SCF crashed (EDDDAV, ZHEGV, "I REFUSE") or OUTCAR is incomplete | ALGO = Normal first, then ALGO = All with EDIFF = 1e-5 |
| FAIL_W90 | DFT is fine, but wannier90 printed no Final State | Switches EOS_SCDM_ONLY off and raises NBANDS to 1.3 times |
| FAIL_LEOS | The Wannier step is fine, but there is no usable EOSCAR (missing, or marked LEOS_ABORTED) | Reruns with EOS_NELECT_ABORT relaxed to 0.5, for diagnosis |
| OK | EOSCAR exists and carries no LEOS_ABORTED mark | |

"Not enough bvectors" from wannier90 counts as FAIL_W90 (a k-mesh problem, not
an SCF problem); the remedy is in patch/README.md, "Common failures".

Two further modes:

- **Tiered screening** (`stage --mode tiered`, then `escalate`): cases without
  4d, 5d or f elements start in the SCDM-only fast mode; cases with them go
  straight to the full MLWF. After the run, `escalate` re-queues with the full
  MLWF every fast-mode case whose EOSCAR carries a quality warning (stray, dNe
  deviation, too many ghosts).
- **Offline wannier90** (`stage --offline-w90`, then `w90` repeatedly): the
  GPU job only exports .mmn/.amn (`LWRITE_MMN_AMN = .TRUE.`,
  `LWANNIER90_RUN = .FALSE.`) and wannier90 runs on the DSPRHBM CPU partition.
  The `w90` command advances each case: export finished → write the .win of
  each spin (tools/offline_w90_prep.py, kept next to batch.py in a single-file
  deployment) and submit one CPU job (both spins concurrently) → once both
  .wout have a Final State, copy them back to the case directory for pyleos.
  An export-only run yields no usable EOSCAR; the offline mode reads the state
  from the exported files and the .wout.

**The three constants VASP_BIN, OFFLINE_VASP_BIN and W90_BIN at the top of
batch.py are hard-coded to paths from before the home-directory migration of
the SAI cluster (/home/sjtu-gongxueqing/zhangxianyang/...), and VASP_BIN still
points at the old vasp.6.4.3_leos_v2. Change them to the current paths before
use; see patch/README.md, "Cluster deployment".**

## Tests

```bash
python3 -m unittest discover pyleos/tests -v   # 23 tests
```

| File | Content |
|---|---|
| test_mcf.py | Minimum-cost-flow solver and assignment |
| test_batch.py, test_batch_tier.py, test_batch_offline.py | Failure triage, retry, tiered screening and offline mode of batch.py |
| test_registry.py | Consistency of cases/INDEX.tsv, cases/STRUCTURES.tsv and cases/ (the same check as `python3 -m pyleos check`) |

## Known limitations

- The assignment between equivalent sites (the two ends of a homonuclear bond,
  equivalent atoms of a cluster) is not unique: the element mean is right, and
  the value of a single atom depends on which optimum the solver lands on. Use
  `--fragment`.
- Capacities follow physical limits only. A diagnostic mode that tightens them
  to chemically common oxidation states is not implemented.
