# Fig. 3 – caption draft and provenance

Script: `fig3_gpuw90.py` (run `python3 fig3_gpuw90.py` in this
directory). Output: `fig3_gpuw90.pdf` / `.png`, page 17.1 cm x 11.4 cm.
Tags M1–M6 are defined in panel (a) and reused in panels (b)–(d).

## Caption

**Fig. 3** Accelerating the Wannierisation step. Grey, squares: reference;
blue, circles: with the modification(s) named by the tags M1–M6.
(a) Stages of the calculation of one frame (top row; M(k,b): overlap
matrices between neighbouring k-points, A(k): SCDM projections that start
the Wannierisation, r_n: centres of the Wannier functions, WFs) and, below
each stage, where it runs. Reference: the whole frame runs in one GPU job,
and one of its processes calls the serial wannier90 library. M1: the GPU
job ends after exporting M and A, and wannier90 runs on a CPU node,
parallelised over k-points with MPI; the assignment follows on a CPU.
M2: each MPI process uses a threaded BLAS/LAPACK library (Intel MKL).
M3: the number of WFs, N_WF, is set separately for each spin channel.
M4: gpuw90, a tensor implementation of the two iterative stages that runs
unchanged on CPU or GPU. M5: both iterative stages are run with zero
iterations and the centres of the SCDM start are used (dashed outline).
M6: the centres of the previously analysed frame are passed to the
localisation as guiding centres. Right: wannier90 treats the N_k N_b
overlap matrices one after another, whereas gpuw90 applies each matrix
product or diagonalisation to all of them in one batched operation on an
array of shape [N_k, N_b, N_WF, N_WF] (N_k k-points, N_b neighbours of each
k-point). (b) Wall time per iteration of disentanglement and of localisation
for one spin channel of a 96-atom Li2MnO3 frame (288 WFs, 356 bands,
8 k-points): wannier90 on one core, with 8 MPI processes (M1) and with
8 MPI processes of 4 threads each (M1 + M2), and gpuw90 (M4) on the CPU
cores and on one V100 GPU of a GPU node. Each line runs from the reference
value to the measured value; the columns give the time and the speed-up
over the reference. The times are per iteration, not times to convergence:
gpuw90 localises by fixed-step steepest descent, wannier90 by conjugate
gradients. (c) Wall time per frame for both spin channels of the same
frame: Wannier step (filled bars) and GPU job (open bars) for the reference,
for M1 (8 MPI processes per spin channel, both channels at the same time),
for M1 + M2 (2 threads per process) and for M1 + M3 (240 instead of 288 WFs
in the second spin channel). In the reference the Wannier step is part of
the GPU job; with M1 the GPU job only computes and exports the overlaps, and
its time is one measurement repeated in the three rows. (d) Top: run time of
wannier90 with the full iterative Wannierisation (squares) and in the
SCDM-only mode M5 (circles) for CSe2 (32 WFs, 48 k-points) and GaSe (76 WFs,
49 k-points); numbers above the lines are the ratios. M5 changes the
assigned oxidation states in both cases. Over 38 test cases it reproduced
the element-mean oxidation states of the full calculation in 32 and all
atomic oxidation states in 19, and it reproduced them in none of 12 frames
of the Li2MnO3 cell, where about 288 unlocalised WFs per frame flag the
failure. Bottom: time of the eight analysed steps of two 40-step ab initio
molecular-dynamics runs of an 88-atom Li1.5MnO3 cell, without (squares) and
with (circles) warm start M6; bars and numbers are the means.

Notes for the text (not for the caption):

- (b) M1 alone gives 7.6x (disentanglement) and 7.5x (localisation) per
  iteration on 8 processes; M1 + M2 on 32 cores gives 18.1x and 12.6x;
  gpuw90 on one V100 gives 12.6x and 22.8x, on CPU 2.4x and 3.1x. For
  disentanglement 32-core wannier90 (0.055 s) is faster than gpuw90 on the
  V100 (0.079 s); for localisation the V100 is faster (0.041 s against
  0.074 s). The figure shows both.
- (c) M1 cuts the Wannier step from 506.5 s to 70.0 s and the GPU job from
  816.9 s to 307.1 s (62.4 % less GPU time); the frame as a whole takes
  307.1 + 70.0 = 377.1 s instead of 816.9 s. M2 and M3 each remove a
  further 8 s and 22 s from the Wannier step. Their combination was not
  measured.
- (d) M5: 33x (CSe2) and 50x (GaSe) on the wannier90 run time. CSe2: full
  C 0 / Se 0, SCDM-only C -4 / Se +2 (the reference answer, by accident
  according to RESULTS.md), 4 unlocalised WFs. GaSe: full Ga +2 / Se -2,
  SCDM-only Ga +1 / Se -1, 28 unlocalised WFs. M6: mean 546 s without and
  572 s with warm start (s.d. 29 s and 34 s over eight steps), i.e. no gain;
  a step without analysis takes about 46 s.

## Modifications (how the list was established)

| Tag | Modification | Where documented | Measured in |
|---|---|---|---|
| M1 | Offline wannier90: export-only GPU job (`LWRITE_MMN_AMN`, `LWANNIER90_RUN = .FALSE.`), then wannier90.x with MPI on a CPU node, both spin channels at once | `pyleos/README.md` l. 160–168, `patch/README.md` l. 79–97, `docs/RESULTS.md` l. 635–661 | (b), (c) |
| M2 | Threaded MKL in the offline wannier90 (MPI processes x threads) | `patch/README.md` l. 88–94, `docs/RESULTS.md` l. 648–658 | (b), (c) |
| M3 | Number of WFs per spin channel (patch v3.1) | `patch/README.md` l. 34–35, l. 198 | (c) |
| M4 | gpuw90: disentanglement and localisation as batched `matmul` / `eigh` over k and b; k-dependent window sizes folded into one batched `eigh`; convergence flag read back only every few iterations | `gpuw90/DESIGN.md` l. 67–83, 111–123, 188–212; `gpuw90/README.md` l. 307–349 | (b) |
| M5 | SCDM-only mode: `num_iter = 0`, `dis_num_iter = 0` written to the wannier90 input | `patch/README.md` l. 46, `patch/v3/apply_v3_mlwf.py` l. 32–36, `docs/RESULTS.md` l. 608–633 | (d) |
| M6 | Warm start: centres of the previously analysed step as wannier90 guiding centres | `patch/README.md` l. 47, `patch/v3/apply_v3_mlwf.py` l. 38–54, `docs/PLAN.md` l. 269–275, `docs/RESULTS.md` l. 678–696 | (d) |

The list given in the task was confirmed. Two refinements: (i) the SCDM
start (the A matrix) is computed by VASP inside the GPU job, not by the
wannier90 library, so the chain in (a) is SCF → overlaps → SCDM start →
[disentanglement → localisation] → assignment, and only the bracketed pair
is wannier90; (ii) M3 was introduced to remove surplus ("ghost") WFs in the
minority-spin channel (a correctness fix); the shorter Wannier step is a
side effect, measured in row D of the pipeline table. One further
speed-related change is documented but has no timing in the repository and
is not drawn: the blocked band assembly in the SCDM gather (commit 1aa7d02,
`patch/README.md` l. 175–177 and 190).

## Provenance

All plotted values are printed by the script in a provenance table.

| Panel | Quantity | Source |
|---|---|---|
| b | wannier90, 1 core | `results/w90_threads_scan.tsv`, row `seq / np 1 / nt 1 / spin 1`: 0.997, 0.935 s per iteration (job 654879) |
| b | M1, 8 MPI processes | same file, row `thr / np 8 / nt 1 / spin 1`: 0.132, 0.125 |
| b | M1 + M2, 8 x 4 | same file, row `thr / np 8 / nt 4 / spin 1`: 0.055, 0.074 |
| b | M4, gpuw90 on CPU and on V100 | **hand-copied** block `README_TIMINGS` in the script, from `gpuw90/README.md` l. 321–322 (section "V100 (SAI cluster)", job 644120, 2026-07-11): 0.423, 0.303 and 0.079, 0.041. No data file holds these numbers; the script stops if the two README rows are no longer present |
| c | Wannier step and GPU job | `results/w90_e2e_pipeline.tsv`, columns `w90_wall_s`, `gpu_leg_s`, rows A (reference), B (M1), C (M1 + M2), D (M1 + M3) |
| d | full Wannierisation, CSe2 and GaSe | "Total Execution Time" of `cases/{cse2,sega}/v31-wp/wannier90.wout`: 3.283, 32.068 s |
| d | SCDM-only, CSe2 and GaSe | "Total Execution Time" of `cases/{cse2,sega}/v3-scdm-only/wannier90.wout`: 0.099, 0.646 s. The script checks that the two runs of a pair have the same POSCAR, k grid, number of bands and of WFs, and that the SCDM-only run reports 0 iterations for both stages |
| d | warm start | `results/cathode_md/warmstart_loop_off.txt`, `warmstart_loop_on.txt`: one time per ionic step, 40 steps; analysed steps are every 5th (`EOS_EVERY_N = 5` in `warmstart_INCAR.on`); the two INCAR files differ only in `EOS_WARM_START` (checked by the script) |

Hand-copied values (the block marked `HAND-COPIED VALUES` in the script):

1. gpuw90 timings of panel (b), `gpuw90/README.md` l. 321–322 – plotted.
2. SCDM-only timing of the Li2MnO3 frame, 554 → 8.3 s per frame,
   `docs/RESULTS.md` l. 616 – **not plotted**. `SHOW_DOCUMENTED_SCDM = True`
   adds it to panel (d) with open markers, a dashed line and an asterisk on
   its label; the caption would then need "* value quoted in the project
   documentation, no data file".

The accuracy figures quoted in the caption for M5 (32/38, 19/38, 0/12, about
288 unlocalised WFs) are from `docs/RESULTS.md` l. 615–616 and 622–624; the
script checks that these lines are still there. They cannot be recomputed:
RESULTS.md l. 633 states that the raw data of the 38-case comparison are not
in the repository, and no SCDM-only run of the Li2MnO3 frames is in `cases/`.

## Verified against the documentation

| Statement | From the data files | Status |
|---|---|---|
| RESULTS.md l. 643–646: 506.5 / 70.0 / 62.0 / 48.1 s; GPU 816.9 → 307.1 s, 62 % less | identical; 62.4 % | agrees |
| RESULTS.md l. 648–651: 8 x 4 threads, 0.055–0.065 and 0.074 s per iteration, "about 16x and 13x" | spin 1: 18.1x, 12.6x; spin 2: 15.9x, 12.6x | agrees (the range covers both spin channels; the figure shows spin 1) |
| RESULTS.md l. 670–672, README l. 309–310: V100 "about 14x and 22x" over serial wannier90 | 12.6x and 22.8x against 0.997 / 0.935 s of the data table | the disentanglement factor differs: the README divides the stage time by the iterations of another serial run (41.886 s / 38 = 1.10 s, `gpuw90/testdata/li2mno3/wannier90.1_ref.wout`) |
| RESULTS.md l. 672–673: gpuw90 on CPU "about 4.6x" | 2.4x and 3.1x on the node of the V100 job | 4.6x compares a laptop with the cluster (README l. 341) |
| RESULTS.md l. 687–688: warm start 572 s / 546 s per analysed step | 571.9 / 546.1 s | agrees |
| RESULTS.md l. 687–688: whole runs 6742.6 s / 6532.4 s | sum of the 40 step times: 6735.1 / 6524.1 s | 7.5 s and 8.3 s more in the document than in the files (not reproducible from them) |
| RESULTS.md l. 681–682: warm-start test on the "96-atom" Li2MnO3 cell | `warmstart_INCAR.on`: x = 0.25, MAGMOM for 24 + 16 + 48 = 88 atoms | differs; the caption says 88 atoms |
| RESULTS.md l. 692–693: an analysed step costs about 525 s more than a plain step | 546.1 - 45.8 = 500 s (off), 571.9 - 45.8 = 526 s (on) | agrees for the run with warm start |
| RESULTS.md l. 615–616: SCDM-only 554 → 8.3 s (67x), CO2 2 → 0.023 s, 86x over 38 cases | no file | not supported by data; not plotted |

## Documented, but not supported by a data file

- gpuw90 timings on the V100 and on CPU (README only) – plotted from the
  hand-copied block.
- All SCDM-only timings quoted in RESULTS.md (554 → 8.3 s, 2 → 0.023 s,
  86x) and the agreement rates (32/38, 19/38, 0/12) – timings not plotted,
  rates quoted in the caption.
- "About half a minute for both spin channels with gpuw90" (RESULTS.md
  l. 673–674, README l. 328–330) – an extrapolation from per-iteration times; not
  plotted.
- 72.2 s for 4 threads per process with both spin channels at once
  (RESULTS.md l. 658) – appears only in a comment line of
  `w90_e2e_pipeline.tsv`; not plotted.
- "About 30 s" fixed time to read the overlap files (RESULTS.md l. 655–656) –
  appears only in the note column of row C; not plotted.
- Speed effect of the blocked SCDM assembly (commit 1aa7d02) – no number
  anywhere; not drawn.

## Caveats

1. **Panel (b) compares runs on different machines.** The wannier90 rows are
   from a Xeon Max CPU node (job 654879), the gpuw90 rows from a V100 node
   (job 644120); the number of CPU threads used by gpuw90 on CPU is not
   recorded. wannier90 times are the per-iteration values parsed from the
   .wout files; gpuw90 times are stage time divided by a fixed number of
   iterations (`gpuw90/bench_gpu.py`; 50 / 200 iterations on the GPU and
   10 / 20 on the CPU in the job template `run_bench_v100.sbatch`). All
   rows are spin
   channel 1, whose disentanglement converges in 38 iterations; in spin
   channel 2 it reaches the 200-iteration limit.
2. **Panel (b) is not time to solution.** gpuw90 has no conjugate-gradient
   or line-search step (`gpuw90/DESIGN.md` l. 483–487) and no adapter that
   writes its centres for the assignment (l. 478–479, 499–500). M4 therefore has no
   entry in panel (c): there is no end-to-end timing of a frame with
   gpuw90.
3. **Panel (c), rows.** Row D of the data file is row B plus the smaller
   spin-2 problem, not row C plus it; "M1 + M2 + M3" was not measured.
   Row C uses 2 threads per process because both spin channels run at once
   (4 threads was slower according to the file comment), whereas the 8 x 4
   row of panel (b) was measured with the spin channels one after another.
   The Wannier step of rows B–D runs on another node after the GPU job, so
   the frame is finished after the sum of the two bars.
4. **Panel (d), SCDM-only pairs.** The two runs of each pair are separate
   jobs with different plane-wave cutoffs (500 eV full, 600 eV SCDM-only);
   the size of the wannier90 problem is identical. The node of the full
   runs is not recorded in the repository (the SCDM-only runs are 4-GPU
   jobs with the wannier90 library). Both systems are among the six of 38
   for which RESULTS.md reports that SCDM-only changes the element means;
   the six `cases/*/v3-scdm-only` runs are these six systems.
   The other four have no like-for-like partner: `mocl2/v3-scdm-only` is not
   spin-polarised while the full MoCl2 run is; `wcl{4,5,6}/v3-scdm-only` differ from
   the full runs in the number of WFs and in spin treatment.
5. **Panel (d), warm start.** `warmstart_loop_*.txt` hold one number per ionic step
   and do not say whether it is CPU or elapsed time. The first step
   (735–739 s) is not an analysed step and is not shown.
   RESULTS.md attributes the absence of a gain to the SCDM start being good
   enough already; this is an interpretation, not shown by the data.
6. **Panel (a).** With M1 the assignment is done offline by pyleos (drawn
   as "CPU"). M6 exists only in the in-situ trajectory mode; gpuw90 does
   not support guiding centres (`gpuw90/DESIGN.md` l. 498).
