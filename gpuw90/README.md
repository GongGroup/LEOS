# gpuw90

gpuw90 implements two steps of wannier90 in PyTorch: disentanglement
(Souza–Marzari–Vanderbilt) and Marzari–Vanderbilt localisation. It reads input in
the standard wannier90 formats and returns Wannier centres and spreads. CPU and CUDA
share the same code. Design, algorithms and item-by-item verification results are in
DESIGN.md.

## Overview

gpuw90 runs offline, outside VASP, and takes over the disentanglement and
localisation steps of wannier90.x. It modifies neither wannier90 nor VASP. The call
sequence depends on the type of band group:

| Band group | Criterion | Call sequence |
|---|---|---|
| Isolated bands | `num_bands == num_wann` | `wannierise()` |
| Entangled bands (metal, or conduction bands included) | `num_bands > num_wann` | `disentangle()` → `wannierise()` |

**An entangled band group must be disentangled first. Passing the full-band M and A
straight to `wannierise()` gives wrong centres and raises no error.** Two guards in
the code prevent this; each raises `ValueError`:

- `io_w90.load_case()` requires an explicit `allow_entangled=True` when
  `num_bands > num_wann`.
- `wannierise()` rejects a non-square `M0` or `U0`.

| Item | Content |
|---|---|
| Input | `seedname.mmn`, `seedname.amn`, `seedname.win`; entangled bands also need `seedname.eig`; a `.win` without an explicit kpoints block also needs `seedname.nnkp` |
| Output | Python objects: the centre (Cartesian, Å) and spread (Ų) of each WF, and the Ω decomposition. No file is written |
| Dependencies | `torch`, `numpy` |

In the code comments, Phase 1 means isolated-band localisation and Phase 2 means
disentanglement.

## Installation

gpuw90 needs only `torch` and `numpy`. A CPU build of torch installed locally runs
every check; it computes the same as the CUDA build.

```bash
# Option A: virtual environment (.venv/ is already in .gitignore)
python3 -m venv .venv
source .venv/bin/activate
pip install --index-url https://download.pytorch.org/whl/cpu torch numpy

# Option B: install into the user directory
pip install --user torch numpy --index-url https://download.pytorch.org/whl/cpu
```

`--index-url .../whl/cpu` selects the CPU build; the CUDA build is not needed
locally.

Running on a GPU needs no code change: add `--device cuda` on the command line, or
pass `device="cuda"` in library calls. A V100 restricts the torch version; see
"Running on a V100 of the SAI cluster" under "Measured performance" below.

## Running the checks

The check scripts compare the gpuw90 results with wannier90 reference output
(`testdata/<case>/<seed>_ref.wout`), then print `VERDICT: PASS` or `VERDICT: FAIL`
and exit with code 0 or 1 respectively. Run all commands below in the `gpuw90/`
directory.

### Isolated bands: check_against_w90.py

Both isolated-band cases should PASS.

```bash
python3 check_against_w90.py testdata/gaas --seed gaas          # GaAs, 8 k-points
python3 check_against_w90.py testdata/diamond --seed diamond    # diamond, 64 k-points
```

The script prints, in order:

1. the shell radii, shell weights and B1 residual of the k mesh;
2. the initial state (after Löwdin orthonormalisation, before optimisation) against
   the reference "Initial State";
3. a table of the centre and spread deviations of each WF;
4. the components of Ω (I, D, OD, total) against the reference;
5. the maximum centre deviation and the verdict.

The criterion is a maximum centre deviation ≤ 1e-4 Å. WFs are paired by nearest
centre, not by index, so nearly degenerate WFs may swap order. The measured maximum
centre deviations are GaAs 8.3e-7 Å, Diamond 3.1e-9 Å; the full numbers are in
"Verification results" of DESIGN.md.

The script handles isolated bands only. An entangled-band case is stopped by the
guard at load time.

| Argument | Default | Meaning |
|---|---|---|
| `datadir` (positional) | `testdata/gaas` | Case directory |
| `--seed` | `gaas` | File-name prefix; must match the files in the directory |
| `--device` | `cpu` | `cpu` or `cuda` |
| `--num-iter` | 2000 | Maximum number of iterations; a well-conditioned system converges in about 10 iterations |
| `--step` | 1.0 | Steepest-descent step factor α; 1.0 is the classic MV step |
| `--conv-tol` | 1e-10 | Convergence threshold, same as the wannier90 default |
| `--conv-window` | 3 | Number of consecutive iterations that must meet the threshold |
| `--verbose` | off | Print Ω and the gradient RMS for the first 3 iterations and every 50 iterations |

### Entangled bands: check_dis_against_w90.py

The entangled-band check uses two spin-polarised cases with a frozen window, 4 spin
channels in all; every one should PASS:

- `testdata/vircef`: 180 bands → 142 WFs;
- `testdata/li2mno3`: 356 bands → 288 WFs, the same system size as the LEOS
  production calculations.

First, fetch the large input files. `.win`, `.eig` and the reference `.wout` are in
git. Each `.mmn` and `.amn` is 10–300 MB and is not committed; pull them once from
the SAI cluster:

```bash
cd testdata && sh fetch_entangled.sh && cd ..
```

The script reaches the SAI cluster over ssh. The ssh host alias and the paths on the
cluster are hard-coded in the script; check them against your own ssh configuration
before running.

Second, run the channels one by one. Without `--wannierise`, each channel takes from
a few seconds to about 1 minute:

```bash
python3 check_dis_against_w90.py testdata/vircef  --seed wannier90.1
python3 check_dis_against_w90.py testdata/vircef  --seed wannier90.2
python3 check_dis_against_w90.py testdata/li2mno3 --seed wannier90.1
python3 check_dis_against_w90.py testdata/li2mno3 --seed wannier90.2
```

The outer window, the frozen window and the disentanglement parameters are read
automatically from `.win`; items absent from `.win` take the wannier90 defaults.

The verdict depends on two gates only, and PASS requires both. The measured numbers
of each channel are in "Verification results" of DESIGN.md.

| Gate | Quantity compared | Threshold |
|---|---|---|
| GATE 1 | Difference between Ω_I after disentanglement and the reference `Final Omega_I` | ≤ 1e-4 Ų |
| GATE 2 | Difference between the WF centres of the initial state after disentanglement and the reference "Initial State" | ≤ 1e-4 Å |

The initial state of GATE 2 is the Löwdin-projection gauge in the disentangled
subspace, with no minimisation. The script also prints
`|Omega_I(initial) - Omega_I(dis)|`. Both quantities are gauge-invariant, so the
difference should be close to 0.

`--wannierise` continues with localisation and prints a final-state comparison,
which does not enter the verdict. The localisation part of the reference `.wout`
uses CG with a fixed 100 iterations and no convergence criterion, so none of the
four reference final states is converged. The final-state differences come from the
different run settings of the two sides and cannot tell whether the implementation
is correct. See "Verification results" of DESIGN.md.

| Argument | Default | Meaning |
|---|---|---|
| `datadir` (positional) | `testdata/vircef` | Case directory |
| `--seed` | `wannier90.1` | File-name prefix; `.1` and `.2` are the two spins |
| `--device` | `cpu` | `cpu` or `cuda` |
| `--dis-num-iter` | `dis_num_iter` from `.win`, else 200 | Maximum number of disentanglement iterations |
| `--wannierise` | off | Continue with localisation and print the final-state comparison, for reference only; slow on CPU with 288 WFs |
| `--num-iter` | 3000 | Maximum number of localisation iterations; used only with `--wannierise` |
| `--step` | 1.0 | Localisation step factor α |
| `--omega-tol` | 1e-4 | GATE 1 threshold (Ų) |
| `--tol-centre` | 1e-4 | GATE 2 threshold (Å) |
| `--verbose` | off | Print per-iteration information |

### Guard unit tests and gradient self-check

The repository CI does not run the gpuw90 tests. Run them by hand after changing
gpuw90:

```bash
python3 test_guards.py      # or: python3 -m pytest test_guards.py -q
python3 fd_check.py         # optional arguments [testdata_dir] [seed]; default testdata/gaas gaas
```

- `test_guards.py` has 5 unit tests that cover the entangled-band guards:
  `load_case()` rejects entangled bands by default; entangled bands without
  `eig_path` raise an error; the normal loading path works; `wannierise()` rejects
  non-square input and input with inconsistent nw.
- `fd_check.py` compares, along a random anti-Hermitian direction X, the analytic
  directional derivative `Re Tr[g† X]` with the central difference
  `[Ω(U e^{εX}) − Ω(U e^{−εX})] / 2ε`. ε is 1e-4, 1e-5 and 1e-6; the expected ratio
  is ≈ +1.0000. Run it first after changing the gradient or porting the code. The
  reason is in "Numerical issues and their handling" of DESIGN.md.

## Library use

The examples below assume the current directory is `gpuw90/`; the modules are
imported from that directory.

### Isolated bands (nb == nw)

```python
import torch
import io_w90
import wannierise as wz

case = io_w90.load_case("seed.win", "seed.mmn", "seed.amn", "seed")
# if .win has no explicit kpoints block, add nnkp_path="seed.nnkp"

M0 = torch.from_numpy(case.M)                       # [nk, nnb, nw, nw]
A = torch.from_numpy(case.A).to(torch.complex128)   # [nk, nw, nw]
U0 = wz.lowdin_orthonormalize(A)                    # initial gauge

res = wz.wannierise(
    M0, U0,
    torch.from_numpy(case.bvec),
    torch.from_numpy(case.wb),
    torch.from_numpy(case.kb_idx),
    num_iter=2000, step=1.0, device="cpu",          # on a GPU: device="cuda"
)
print(res.centres)   # [nw, 3], Cartesian, Å
print(res.spreads)   # [nw], Å^2
print(res.omega)     # keys "I", "OD", "D", "total"
```

### Entangled bands (nb > nw)

```python
import torch
import io_w90
import wannierise as wz
import disentangle as dz

case = io_w90.load_case("seed.win", "seed.mmn", "seed.amn", "seed",
                        eig_path="seed.eig", allow_entangled=True)
# without allow_entangled=True this raises ValueError; entangled bands require eig_path

dres = dz.disentangle_case(case, device="cpu")   # windows and iteration parameters from .win, else wannier90 defaults
print(dres.omega_I, dres.converged, dres.n_iter)  # omega_I is gauge-invariant; usable for regression comparison

res = wz.wannierise(
    dres.M0_reduced, dres.U0,                     # already reduced to nw dimensions
    torch.from_numpy(case.bvec),
    torch.from_numpy(case.wb),
    torch.from_numpy(case.kb_idx),
    device="cpu",
)
```

To set the windows, the mixing or the convergence parameters explicitly, call
`dz.disentangle(...)` directly. Its signature and the mapping to the `.win` keywords
are in "Algorithms" of DESIGN.md.

### Main functions and return values

```python
io_w90.load_case(win_path, mmn_path, amn_path, seedname="case",
                 nnkp_path=None, eig_path=None, allow_entangled=False)     # -> W90Case

disentangle.disentangle_case(case, device="cpu", num_iter=None, mix_ratio=None,
                             conv_tol=None, conv_window=None,
                             check_every=1, verbose=False)                  # -> DisentangleResult

wannierise.wannierise(M0, U0, bvec, wb, kb_idx, num_iter=1000, step=1.0,
                      conv_tol=1e-10, conv_window=3, device="cpu",
                      dtype=torch.complex128, verbose=False,
                      check_every=1)                                        # -> WannieriseResult
```

The `num_iter` default of the library function `wannierise()` is 1000; that of
`check_against_w90.py` is 2000.

| Return type | Fields |
|---|---|
| `W90Case` | `seedname`, `M`, `A`, `bvec`, `wb`, `kb_idx`, `kpoints`, `real_lattice`, `recip_lattice`, `kmesh_info`, `win` (the parsed `.win` keywords), `eig`; and `num_bands`, `num_wann`, `num_kpts`, `nntot` |
| `DisentangleResult` | `M0_reduced [nk,nnb,nw,nw]`, `U0 [nk,nw,nw]`, `V [nk,nb,nw]`, `omega_I` (Ų), `history`, `n_iter`, `converged`, `ndimwin [nk]`, `ndimfroz [nk]` |
| `WannieriseResult` | `U [nk,nw,nw]`, `centres [nw,3]` (Å), `spreads [nw]`, `omega`, `history`, `n_iter`, `converged` |

`history` is a list of per-iteration records: each entry is
`(iter, Ω_I(i−1), Ω_I(i), delta)` in `DisentangleResult` and
`(iter, Ω_total, rms_grad)` in `WannieriseResult`.

Both `wannierise()` and `disentangle()` take the argument `check_every=N`. The
convergence counter stays on the device, and the convergence flag is read back to
the host only every N iterations. On a GPU, 10 is recommended. With N > 1 the loop
runs at most N − 1 iterations past actual convergence.

## Test data

`testdata/` holds four check cases, each with the wannier90 reference output
`<seed>_ref.wout`.

| Directory | Type | nb / nw | k-points | Input source | Reference output |
|---|---|---|---|---|---|
| `gaas/` | isolated bands | 4 / 4 | 8 | `testw90_example01` of the wannier90 3.1 test-suite | `benchmark.out` of that example, fully converged |
| `diamond/` | isolated bands | 4 / 4 | 64 | `testw90_example05` of the wannier90 3.1 test-suite | `benchmark.out` of that example, fully converged |
| `vircef/` | entangled bands, two spins | 180 / 142 | 8 | VASP `LWRITE_MMN_AMN` export under `~/LEOS/w90_offline/offline_runs/vircef` on the SAI cluster | serial wannier90.x runs in `run_np1_s{1,2}/` of the same directory |
| `li2mno3/` | entangled bands, two spins | 356 / 288 | 8 | VASP `LWRITE_MMN_AMN` export under `~/LEOS/w90_offline/bench/li2mno3_x000_f000` on the SAI cluster | serial wannier90.x runs in `run_np1_s{1,2}/` of the same directory |

| Directory | Files in git | Files to fetch |
|---|---|---|
| `gaas/` | `gaas.{win,mmn,amn}`, `gaas_ref.wout` | none |
| `diamond/` | `diamond.{win,mmn,amn}`, `diamond_ref.wout` | none |
| `vircef/` | `wannier90.{1,2}.{win,eig}`, `wannier90.{1,2}_ref.wout` | `wannier90.{1,2}.{mmn,amn}` |
| `li2mno3/` | `wannier90.{1,2}.{win,eig}`, `wannier90.{1,2}_ref.wout` | `wannier90.{1,2}.{mmn,amn}` |

- The entangled-band reference runs used symbolic links to the same files, so gpuw90
  and the reference read exactly the same input.
- The files to fetch are excluded by `.gitignore`; download them with
  `testdata/fetch_entangled.sh`.

## Measured performance

On the li2mno3 case, a V100 is about 14 times (disentanglement) and about 22 times
(localisation) faster per iteration than serial wannier90.x on the cluster. All
timings below are times per iteration (s/iteration) at a fixed number of iterations,
so they compare directly across devices.

### V100 (SAI cluster)

Measured on 2026-07-11, SAI cluster job 644120. The case is the li2mno3 production
frame: 288 WFs, 356 bands, 8 k-points.

| Stage | V100 | torch-CPU on the compute node | Serial wannier90.x on the cluster | V100 speed-up over serial |
|---|---|---|---|---|
| disentangle | 0.079 | 0.423 | ~1.10 | ~14× |
| wannierise | 0.041 | 0.303 | ~0.92 | ~22× |

Unit: s/iteration.

- The check also passes on CUDA: |ΔΩ_I| = 4.4e-9 Ų and a maximum initial-state
  centre deviation of 8.3e-7 Å, the same level as the CPU result.
- Taking the 554 s/frame of the "Efficiency" section of docs/RESULTS.md as the
  baseline (wannier90 part), the Wannierisation part for both spins scales to about
  half a minute. This is an estimate scaled by the time per iteration.

### Local CPU baseline

Measured on 2026-07-11. Environment: Apple M-series laptop, torch 2.8.0 CPU build,
FP64 throughout. The serial wannier90.x column is taken from the timings (single
process) in the reference outputs
`testdata/{li2mno3,vircef}/wannier90.1_ref.wout`.

| Stage | Case | gpuw90 CPU (s/iteration) | Serial wannier90.x on the cluster (s/iteration) | Note |
|---|---|---|---|---|
| disentangle | li2mno3, 288 WFs | 0.24 (38 iterations, 9.1 s) | 1.10 (38 iterations, 41.9 s) | same number of iterations, 4.6× |
| wannierise | li2mno3, 288 WFs | 0.196 | ~0.92 (CG, 100 iterations in about 92 s) | same kind of work per iteration |
| disentangle | vircef, 142 WFs | 0.07 | ~0.21 (82 iterations, 17.1 s) | |
| wannierise | vircef, 142 WFs | 0.051 | ~0.18 | |
| `.mmn` read | li2mno3 (258 MB of text) | 3.4 s | 4.7 s | one-off cost |

The local CPU is faster than serial wannier90.x because batched matmul uses
multithreaded BLAS. The gain from a GPU must therefore be judged against both
baselines, torch-CPU and serial wannier90.x; the V100 table lists both columns.

### Reproducing the timings

```bash
python3 bench_gpu.py testdata/li2mno3 --seed wannier90.1                 # local CPU
python3 bench_gpu.py testdata/li2mno3 --seed wannier90.1 --device cuda   # GPU
```

`bench_gpu.py` runs disentanglement and then localisation on one spin channel and
prints the total time and the s/iteration of each stage.

- The number of iterations is fixed: the script sets `conv_tol=0` internally, so it
  never converges early.
- The script does not judge correctness. Pass `check_dis_against_w90.py` before
  timing.
- The Ω_I printed at the end is the value at the fixed number of iterations, not the
  converged value.

| Argument | Default | Meaning |
|---|---|---|
| `datadir` (positional) | required | Case directory |
| `--seed` | `wannier90.1` | File-name prefix |
| `--device` | `cpu` | `cpu` or `cuda` |
| `--dis-iters` | 50 | Fixed number of disentanglement iterations |
| `--wann-iters` | 100 | Fixed number of localisation iterations |
| `--check-every` | 10 | Interval at which the hot loop reads back the convergence flag |

### Running on a V100 of the SAI cluster

`run_bench_v100.sbatch` is the job template for the V100 timing. It runs the
following in order:

1. `check_dis_against_w90.py` on CUDA for li2mno3 spin 1; it exits on failure;
2. timing on CUDA, `--dis-iters 50 --wann-iters 200`;
3. timing on CPU on the same compute node, `--dis-iters 10 --wann-iters 20`.

- The code path defaults to `$HOME/LEOS/gpuw90` (rsync the local `gpuw90/` there)
  and the data path to `$HOME/LEOS/w90_offline/bench/li2mno3_x000_f000`. Change them
  in the `GPUW90` and `DATA` variables at the top of the script.
- The script symlinks the reference `.wout` under `run_np1_s{1,2}/` to
  `<seed>_ref.wout` in the data directory.
- The job output goes to `bench-out` in the submission directory.

The requirements of the cluster environment, also noted in the template comments:

| Item | Requirement |
|---|---|
| Partition and QOS | Every job in the `16V100` partition must request at least 4 GPUs. The template requests 4 and uses only 1 (`CUDA_VISIBLE_DEVICES=0`); the QOS is `flood-gpu` and the time limit 4 hours |
| Resource options | Do not add `--cpus-per-task` or `--mem*`. A site hook automatically cancels GPU jobs that carry these options |
| torch version | The V100 is sm_70. torch 2.9.1+cu128 in the site mace environment does not include sm_70 and fails at run time with "no kernel image". A wheel of the cu124 series is needed (torch ≤ 2.6). The template checks this with a `torch.cuda.get_arch_list()` assertion |
| Python environment | mamba environment `gpuw90` (python 3.11); the one-off install command is in the template comments |

## Files

| File | Purpose |
|---|---|
| `io_w90.py` | Parses `.mmn`, `.amn`, `.eig`, `.win`, `.nnkp`; rebuilds the b vectors and shell weights; parses the reference `.wout`; `load_case()` holds the entangled-band guard |
| `wannierise.py` | Batched MV localisation: Löwdin initial gauge, Ω decomposition, gradient, matrix exponential, steepest-descent loop; rejects non-square input |
| `disentangle.py` | Batched SMV disentanglement: window and frozen masks, initial subspace, fixed-point iteration on the Z matrix; outputs `M0_reduced` and `U0` |
| `check_against_w90.py` | Isolated-band check |
| `check_dis_against_w90.py` | Entangled-band check (two gates: Ω_I and the initial state after disentanglement) |
| `test_guards.py` | Unit tests of the entangled-band guards |
| `fd_check.py` | Finite-difference self-check of the sign and prefactor of the gradient |
| `bench_gpu.py` | Timing script with a fixed number of iterations |
| `run_bench_v100.sbatch` | Job template for the V100 timing on the SAI cluster |
| `DESIGN.md` | Design, algorithms, numerical issues, verification results, GPU memory estimate, connection to LEOS, known limitations |
| `testdata/` | Four check cases, and `fetch_entangled.sh` to fetch the large inputs |
