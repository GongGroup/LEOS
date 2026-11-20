# LEOS VASP patch

`leos_v3.patch` is a single diff against stock VASP 6.4.3. It changes 4 files:
`src/.objects`, `src/eos.F` (new), `src/main.F`, `src/mlwf.F`. With the patch, VASP
calls the wannier90 library at the end of the calculation to wannierise, and writes
the EOSCAR file. The current version is v3.1 plus the SCDM blocked assembly (see
"Version history").

## Installation

Prerequisite: the **serial** static library of wannier90 3.1 (see the first item of
"Running notes").

```bash
tar xzf vasp.6.4.3.tgz && cd vasp.6.4.3/..   # top-level directory must be named vasp.6.4.3
patch -p1 < leos_v3.patch
cd vasp.6.4.3
# add to makefile.include:
#   CPP_OPTIONS += -DVASP2WANNIER90
#   WANNIER90_ROOT ?= /path/to/wannier90-3.1.0
#   LLIBS += -L$(WANNIER90_ROOT)/lib -lwannier
#   INCS  += -I$(WANNIER90_ROOT)/include
make DEPS=1 -j8 std gam
```

Build environment on the SAI cluster: nvhpc 24.3 + openmpi 5.0.8 + fftw 3.3.10
(modules `nvhpc/24.3-pgi-tuned`, `openmpi/5.0.8-nvhpc24.3-pgi-auto`, `fftw/3.3.10`);
a full build takes about 5–11 minutes.

## INCAR tags

The minimal setting is `LEOS = .TRUE.`. With it on, the following four VASP tags are
set to `.TRUE.` automatically unless INCAR gives them: `LWANNIER90`,
`LWANNIER90_RUN`, `LSCDM`, `LWANNIER90_AUTO_WINDOW`. When `NUM_WANN` is not given,
it is set automatically from the number of occupied bands of each spin.

| Tag | Default | Effect |
|---|---|---|
| LEOS | .FALSE. | Master switch |
| EOS_FLOW | .TRUE. | Do the minimum-cost flow assignment inside VASP and write the result to the OS(flow) column of EOSCAR (β = 10, Allen electronegativity, same as the pyleos defaults) |
| EOS_WF_DETAIL | .TRUE. | Write the #WF block: centre, spread and occupation of every WF. pyleos depends on it |
| EOS_PROJ_WEIGHT | .TRUE. | Write the #WP block: weight of every WF on the PAW projectors of each atom, used by projection-assisted arbitration |
| EOS_EVERY_N | 0 | 0: analyse only at the last step (single-point calculations included); 1: analyse every step; N > 1: analyse every N steps and at the last step |
| EOS_NELECT_WARN | 0.01 | Warn when the Wannier electron count differs from NELECT by more than this |
| EOS_NELECT_ABORT | 0.35 | When the difference exceeds this, write a diagnostic EOSCAR and abort |
| EOS_SCDM_ONLY | .FALSE. | Skip the iterative localisation of wannier90 and use only the centres of the SCDM initial guess (for fast screening) |
| EOS_WARM_START | .TRUE. | When EOS_EVERY_N > 0, pass the WF centres of the previous analysed step to wannier90 as guiding centres |
| EOS_ENEG_TOL | 0.05 | Greedy assignment: electronegativities that differ by less than this count as equal |
| EOS_ONSITE_CUT | 0.35 | Greedy assignment: a WF centre closer to an atom than this distance (Å) is a lone pair of that atom |
| EOS_O_MAX | 8.0 | Greedy assignment: the most electrons O may receive |
| EOS_ENEG_PAULING | .FALSE. | Greedy assignment uses Pauling electronegativities instead |

The last four tags affect only the greedy column OS(int) of EOSCAR, not OS(flow).

A working single-point INCAR (cases/as4s4/v3/INCAR):

```
PREC = Accurate
ENCUT = 500
EDIFF = 1e-06
ISMEAR = 0
SIGMA = 0.05
ISYM = 0
NSW = 0
IBRION = -1
LWAVE = .FALSE.
LCHARG = .FALSE.
LEOS = .TRUE.
LSCDM = .TRUE.
LWANNIER90_AUTO_WINDOW = .TRUE.
EOS_WF_DETAIL = .TRUE.
EOS_PROJ_WEIGHT = .TRUE.
```

## Running notes

### The wannier90 library and offline wannier90

- **Library mode must link the serial libwannier.a.** The patch lets only the master
  process enter the wannier90 library; the other processes wait outside. With an MPI
  build of the library, the library starts collective communication of its own,
  which waits on the processes outside while they wait on it: a deadlock at the
  first disentanglement step, with every process spinning at 100% CPU and no .werr
  file. Many sites ship precompiled MPI builds (for example every version under
  /data/app on Chinahpc). To check: `nm libwannier.a | grep -ci mpi_`; a non-zero
  result means an MPI build, and the serial library must be rebuilt with `COMMS=`
  (left empty).
- **Build the offline wannier90.x with MPI + threaded MKL.** The restriction above
  applies to library mode only. Key make.inc settings: `F90=ifort`, `COMMS=mpi`,
  `MPIF90=mpiifort`, and LIBS with
  `-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5`. At run time set
  `MKL_NUM_THREADS=<n> MKL_DYNAMIC=FALSE`. When the two spins run concurrently,
  2 threads per process is fastest (data under "Efficiency" in docs/RESULTS.md).
  wannier90 built with gnu + OpenMPI 5 crashes with several processes; use ifort +
  Intel MPI + MKL. The ready build on the SAI cluster is
  `w90_offline/build/wannier90-3.1.0-impi-thread/` under the LEOS directory; the
  offline workflow is described in `w90_offline/README.md`.

### Common failures

- **NELECT mismatch too large.** The cause is almost always a missing frozen window
  in the disentanglement. Set `LWANNIER90_AUTO_WINDOW = .TRUE.` (together with
  LSCDM); do not write wannier90.win by hand. On abort the patch writes a diagnostic
  EOSCAR (marked `# LEOS_ABORTED`) whose `#NLD` block lists the (spin, k-point,
  band) entries with the largest loss. If the loss is concentrated at the top of the
  occupied bands and spread over all k-points, the frozen window is missing.
- **`kmesh_get_bvector: Not enough bvectors found`.** Small body-centred or oblique
  cells, and the large cells of some molecular crystals, fail the b-vector check of
  wannier90 on particular k-meshes. Try in turn: an even mesh or one that matches
  the cell shape; or add `WANNIER90_WIN = "\nkmesh_tol = 0.00001\n"` to INCAR (this
  relaxes the default 1e-6 to 1e-5; the multi-line string can also carry other
  parameters such as dis_num_iter). Do not switch to Gamma-only to avoid the error:
  localisation behaves differently in Gamma mode and produces many stray WFs.
- **#WP data.** Before the fix of 2026-07-11 (commit 751389a), the #WP blocks of all
  multi-process runs are wrong: with band parallelisation the squared modulus was
  taken before the sum over processes, which dropped the cross terms. Such data must
  be recomputed with the fixed version. In addition, **the #WP of the GPU (nvhpc)
  build and of the CPU (Intel) build differ systematically**: per-atom aggregates
  differ by about a factor of two, for a reason not yet found. #WP from different
  platforms must therefore not be mixed; calculations that need #WP are all done on
  one platform.

### Memory and SCF

- **KPAR and memory (GPU build).** KPAR does not store the wavefunctions separately
  per k-point; instead the memory of each process is proportional to KPAR / number
  of processes. With few processes, KPAR > 1 exhausts the host memory: a Ce system
  on 4 GPUs with KPAR = 2 runs out of memory at start-up. Use KPAR = 1 with 4 GPUs /
  62.5 GB; KPAR = 2 works with 8 GPUs / 125 GB.
- **Give the full MAGMOM for magnetic systems.** This keeps the default rule from
  raising NBANDS (from 332 down to 308 for the Ce system) and makes the SCF converge
  faster.
- **Measured configuration for the Ce18O36Co4CH4 final state.** KPAR = 1,
  NBANDS = 308, full MAGMOM: 5 hours 24 minutes on 4 GPUs, peak memory 49.2 GB (79%
  of the quota). The disentanglement of each spin converges in 100 steps (with a
  frozen window shared by the spins, the second spin ran the full 1000 steps, 1743
  seconds, without converging); the per-spin num_wann is 255 / 246; dNe falls from
  0.193 to 0.0008 and ghost from 6 to 1.
- **SCF settings for the Ce system (four settings compared, Ce initial state).**
  ALGO = ALL is required; Davidson diverges to the order of −10⁷ eV within 176
  steps. Under ALGO = ALL the mixing parameters matter little: 545 steps with
  conservative settings, 526 with default mixing, 511 with mild Kerker, a difference
  of 3–6%. The only effective speed-up is a full initial MAGMOM (final state: from
  539 steps down to 454, 16% fewer). Reusing wavefunctions between images is not yet
  tested (issue #18).

### SLURM jobs

- **Relative paths.** A relative path such as `--output=print-out` is resolved
  against the directory of submission, and so is the working directory of the job.
  If you submit someone else's script from the home directory, VASP looks for INCAR
  there and exits at once, with the error message in ~/print-out. cd to the case
  directory before submitting, or use absolute paths and an explicit cd in the
  script.
- **Put resource options in the script.** The check hook of the SAI cluster matches
  the submission command line: a job given `--cpus-per-task` on the command line is
  cancelled by root within 2 seconds, without even an output file; the same value in
  a `#SBATCH` line of the script works. `--export=ALL,VAR=x` is not affected.
- **Rebuilding the binary.** If queued jobs refer to the same binary path, they may
  be scheduled between the removal of the old directory and the end of the build,
  and fail at once (running jobs are not affected). Check with squeue before
  rebuilding, or deploy to a new path.

### Development notes

- The count argument n of the VASP MPI sum routines M_sum_z8 / M_sum_d8 is
  INTEGER(qi8); call them with `SIZE(x, KIND=qi8)`. With `SIZE(x)` (32-bit) the high
  bits are stack garbage, which causes heap overruns and SIGKILL; it does not show
  in serial runs and is hard to trace.

## Cluster deployment

The patched VASP on the SAI cluster is in
`/org/sjtu-gongxueqing/zhangxianyang/temp/org/apps/vasp.6.4.3_leos_v31_20260711/bin/`
(vasp_std, vasp_gam), built from leos_v3.patch at commit 751389a. The SCDM blocked
assembly added later (commit 1aa7d02) changes only speed, not results, and is not
yet redeployed.

The symlink `/org/sjtu-gongxueqing/zhangxianyang/temp/org/apps/leos_current` and the
old job scripts under the LEOS directory still point to the path from before the
home-directory migration, `/home/sjtu-gongxueqing/zhangxianyang/`; change it to
`/org/sjtu-gongxueqing/zhangxianyang/temp/org/` before use.

## Version history

### After v3.1

| Change | Content |
|---|---|
| SCDM blocked assembly (commit 1aa7d02, issue #17) | When SCDM gathers the orbitals, bands are assembled in groups of about 16, and one collective communication replaces the per-band broadcast and reduction. Each buffer element has a single non-zero contribution, so results are bit-identical |

### v3 → v3.1

| Change | Content |
|---|---|
| Diagnostic EOSCAR before abort | When the electron loss exceeds EOS_NELECT_ABORT, an EOSCAR with the `# LEOS_ABORTED` mark and the full loss record is written before skipping. pyleos explains it when it reads one; single aborted frames of a trajectory do not affect the other frames |
| Locating diffuse loss: `#NLD` | Without a frozen window the occupied bands leak continuously at all k-points; this is the most common kind of loss. The `#NLD` block records where the loss is; MUMNOX gives 1348 records, and 4.35 electrons are located at "spin 1 band 117, spin 2 band 124, spread over all k-points". The `#NL` block for window-type loss is also written to EOSCAR |
| Per-spin num_wann | The minority spin of magnetic systems no longer produces ghosts (the ghost of the second spin of VIRCEF goes to zero and Ce returns to +3). The automatic frozen window is also set per spin; if rewriting the .win copy fails, the shared value is used; a NUM_WANN given explicitly in INCAR keeps its former meaning |
| SCDM no longer double-buffers | Process 0 no longer holds two full-grid matrices at once (for the Ce system that memory increment was about 8 GB); results are bit-identical |
| M/A matrices freed early | Non-master processes free M right after the reduction (−0.9 GB per process); M/A are freed after wannier90 finishes (−1 GB resident per process); A is freed only on the SCDM path, so LOCPROJ and elphon are not affected |
| In situ minimum-cost flow: the OS(flow) column | The last column of the per-atom table of EOSCAR is the minimum-cost flow assignment done inside VASP (β = 10, Allen electronegativity, lanthanides with the same electronegativity table as pyleos); the switch EOS_FLOW is on by default. Identical to pyleos atom by atom: in 7 tests the cost is bit-identical to the Python reference, and CO2 and VIRCEF are verified in real runs |
| #WP band-group reduction fix (commit 751389a) | With multi-process band parallelisation, the complex components of the band groups are summed across groups first, and the squared modulus is taken afterwards. Serial results are bit-identical; after the fix the per-atom aggregates of serial and multi-process runs differ by 2×10⁻⁴ (the level of convergence noise); the regression passes on two clusters |
| Batch runs and offline workflow | See "batch.py: cluster batch runs" in pyleos/README.md, and tools/offline_w90_prep.py |

### v2 → v3

| Change | Content |
|---|---|
| Hungarian assignment branch removed | Design flaw (an anion cannot receive more electrons than the neutral atom); never used |
| Electron-loss diagnosis by (spin, k-point, band) | Reports not only the total; the loss can be located |
| #WF block of EOSCAR | Centre, spread and occupation of every WF; the input of pyleos |
| #WP projection-weight block | Data source of projection-assisted arbitration |
| EOS_SCDM_ONLY | Fast screening mode; on CO2 the wannier90 time is 2 s → 0.023 s, oxidation states unchanged |
| Trajectory mode wired in | The main program had `IF (NSTEP == 1)` hard-coded, so EOS_EVERY_N never took effect; analysis now runs step by step under EOS_DO_THIS_STEP (for users who set only LWANNIER90 it still runs at the first step only), and the dangling reference from the unconditional MLWF_FREE at every step is fixed |
| Warm start | WF centres of the previous step as guiding centres |
| Other | Periodic wrap-around formula, automatic-window tolerance, compiler warnings |

## Verification record

### v3 (2026-07-07, SAI cluster)

1. Source: the diff between stock 6.4.3 and the production tree `vasp.6.4.3_leos_v2`
   (after stripping CCQN). The production tree is the source of all M4/M5 data, so
   it is the reference.
2. Closure: stock + leos_v3.patch is byte-identical to the stripped production tree;
   the difference between the production tree and the stripped version is exactly
   the CCQN file set, with no LEOS content mixed in.
3. Build: a fresh build on the patched stock tree passes (job 594329, about 11
   minutes).
4. Regression: the same CO2 and TiO2 inputs were run with the production binary and
   with the new binary; the integer oxidation states in EOSCAR are identical and the
   numeric fields agree within 1e-6 (the differences are last-digit floating-point
   differences between builds); the #WP and #WF blocks are complete.

Findings during the merge:

- The production tree once also carried the patch of another project (CCQN-TS);
  leos_v3.patch has it stripped.
- The production tree had one uncommitted change, `IF (LEOS) LSCDM = .TRUE.`; it is
  now in leos_v3.patch.
- The trajectory part of main.F in the old leos.patch (v2) in the repository is an
  old development version; it lacks the NSTEP == 1 protection for runs that set only
  LWANNIER90, and does not match production.

### v3.1 (2026-07-10, SAI cluster)

1. Build: stock tree + this patch built successfully three times (jobs 612186,
   612454, 612747, about 5–11 minutes each). The changes went through two rounds of
   code review, recorded in the git history.
2. Regression: the CO2 integer oxidation states are identical and the numeric fields
   differ at the 1e-14 level (the criterion remains "integers identical, numbers
   within 1e-6"; the GPU build has last-digit floating-point differences between
   nodes and cannot be byte-identical). For VIRCEF the ghost of the second spin
   falls from 2 to 0 and Ce returns to +3. MUMNOX verified `#NLD` and the abort
   behaviour under two settings, no window and a wrong window, respectively.
3. Memory: Ce18O36Co4 (63 atoms, ISPIN = 2) on 4 GPUs with KPAR = 1 and
   NBANDS = 308 uses 41 GB of VASP memory in total (66% of the 62.5 GB quota).
4. Regression data: the runs named v31, v31-4gpu and v31-abort2 in cases/INDEX.tsv.

## Archived files

- `leos.patch`: the v2 patch, archive only. Its main.F part is an old version (see
  above); do not use it.
- `v3/`: the v2 → v3 incremental material (full replacement of eos.F and injection
  scripts for main.F / mlwf.F), merged into leos_v3.patch; see
  [v3/README.md](v3/README.md).
