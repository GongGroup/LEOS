# gpuw90 design

gpuw90 is a post-processing tool independent of pyleos. It reads input in the
standard wannier90 formats, implements Souza–Marzari–Vanderbilt (SMV)
disentanglement and Marzari–Vanderbilt (MV) localisation as batched PyTorch
operations, and outputs the Wannier centres and spreads from which LEOS assigns
oxidation states. Installation, commands, test data and performance numbers are in
README.md.

## Scope

gpuw90 implements only the disentanglement and localisation steps of wannier90. It
does not fork wannier90 and does not change VASP; the `device` argument selects the
device, and CPU and CUDA share the same code.

| | Isolated bands | Entangled bands |
|---|---|---|
| Condition | `nb == nw` | `nb > nw` (metal, or conduction bands included) |
| Entry point | `wannierise()` | `disentangle()` → `wannierise()` |
| Input files | `.mmn` (nb×nb overlaps), `.amn` (nb×nw initial guess), `.win` | plus `.eig` (band energies, which set the outer and frozen windows) |
| Main operations per k per iteration | 1 Hermitian eigendecomposition of size nw×nw, about 20 matrix products of cost nw³ | the disentanglement stage adds 1 Hermitian eigendecomposition of dimension nb (subspace selection from the Z matrix) |
| Iteration method | fixed-step steepest descent (SD) | disentanglement by fixed-point iteration on the Z matrix, with the same mixing as wannier90; then as for isolated bands |

Out of scope:

- The SCDM initial guess. The initial guess is read from the `.amn` that VASP writes.
- Band interpolation. gpuw90 therefore omits the gauge rotation that wannier90
  applies after disentanglement converges; see "Differences from wannier90" under
  "Algorithms".

Disentanglement is a preprocessing step of localisation. It outputs the reduced
`M0_reduced [nk,nnb,nw,nw]` and the initial gauge `U0 [nk,nw,nw]`, which go straight
to `wannierise()`. `wannierise()` is unaware of nb; isolated and entangled bands
share the same localisation code.

## Modules and data flow

The data start as files written by VASP, are read by `io_w90`, follow one of two
paths according to the type of band group, and all end in `wannierise()`.

```
VASP (LWRITE_MMN_AMN)
 │  writes seedname.mmn, .amn, .eig; a .win is supplied as well
 ▼
io_w90.load_case()
 │  parses the files, rebuilds the b vectors and shell weights, returns W90Case
 │
 ├─ isolated bands (nb == nw)
 │    M0 = case.M
 │    U0 = wannierise.lowdin_orthonormalize(A)
 │
 ├─ entangled bands (nb > nw; load_case needs allow_entangled=True and eig_path)
 │    disentangle.disentangle_case()
 │      window and frozen masks → initial subspace → fixed-point iteration on the Z matrix
 │      gives M0_reduced and U0, the M0 and U0 of the next step
 │
 ▼
wannierise.wannierise(M0, U0, bvec, wb, kb_idx)
 │  each iteration: rotate M → gradient g → exp(dW) → update U
 │  returns centres [nw,3], spreads [nw], omega
 ▼
check_against_w90.py (isolated bands) / check_dis_against_w90.py (entangled bands)
    reads the reference <seed>_ref.wout (io_w90.parse_wout_final_state, parse_wout_dis_info)
    pairs by nearest centre, compares minimum-image distances, verdict PASS/FAIL (threshold 1e-4 Å)
```

### Batched layout

k is the batch dimension. All heavy computation is `torch.matmul` and
`torch.linalg.eigh` batched over k (and over the neighbours b); these correspond to
the per-k `zgemm` and `zheev` calls of wannier90. Switching to `device="cuda"`
leaves the mathematics unchanged.

| Tensor | Shape | Type | Meaning |
|---|---|---|---|
| `M0` | `[nk, nnb, nw, nw]` | complex128 | First-principles overlap matrix of each (k, b) |
| `U` | `[nk, nw, nw]` | complex128 | Unitary gauge at each k |
| `bvec` | `[nk, nnb, 3]` | float64 | Cartesian b vectors (Å⁻¹) |
| `wb` | `[nk, nnb]` | float64 | Shell weights |
| `kb_idx` | `[nk, nnb]` | long | Index of the neighbour k, i.e. the k that k+b folds back to in the Brillouin zone |

The input of the disentanglement stage has the full band dimension:
`M0_full [nk,nnb,nb,nb]`, `A_full [nk,nb,nw]`, `eig [nk,nb]`.

## Algorithms

### Localisation (wannierise.py)

Localisation is the same as the MV steepest descent of wannier90 3.1. Each iteration
runs the following steps batched over all k:

1. Rotate the overlaps: `M^{k,b} = U_k† M0^{k,b} U_{k+b}` (`_rotate_M`; the
   neighbour gauges are taken with `U[kb_idx]`).
2. Phases: `Im ln M_nn = atan2(Im, Re)`, principal branch.
3. Centres: `r_n = −(1/N) Σ_{k,b} w_b b · Im ln M_nn`.
4. Gradient, i.e. the true dΩ/dW, anti-Hermitian:
   `g_k = −(1/N) · 4 Σ_b w_b (A[R^{k,b}] − S[T^{k,b}])`,
   where `R_mn = M_mn·conj(M_nn)`, `T_mn = (M_mn/M_nn)·q_n`,
   `q_n = Im ln M_nn + b·r_n`, `A[B] = (B − B†)/2`, `S[B] = (B + B†)/(2i)`.
5. Step: `dW_k = −(α·N / (4 w_tot)) g_k`, `w_tot = Σ_b w_b`. α = 1 is the classic
   MV step.
6. Update: `U_k ← U_k · exp(dW_k)`. The matrix exponential is computed from the
   Hermitian eigendecomposition of `i·dW`: `exp(dW) = V·diag(e^{−iλ})·V†`, as
   wannier90 does with `zheev`.
7. Convergence: `|ΔΩ| < conv_tol` for `conv_window` consecutive iterations, the same
   semantics as wannier90.

The initial gauge for isolated bands is the Löwdin symmetric orthonormalisation
`U0 = A(A†A)^{−1/2}`, which corresponds to the "Initial State" printed by wannier90.

### Convergence check and device synchronisation

The hot loop does not read data back to the host at every iteration; on a GPU that
would cost one synchronisation per iteration.

- The convergence counter `n_below` and the Ω history are kept in device tensors and
  fetched once after the loop ends.
- Only one scalar convergence flag is read back, every `check_every` iterations. With
  `check_every > 1` the loop runs at most `check_every − 1` iterations past actual
  convergence.
- The M rotated after the update of U is the input of the next iteration, so each
  iteration rotates only once.
- The disentanglement loop works the same way.

### Disentanglement (disentangle.py)

Disentanglement selects at each k the nw-dimensional subspace that minimises the
gauge-invariant spread Ω_I. The signature of the low-level function:

```python
res = disentangle(M0_full, A_full, eig, bvec, wb, kb_idx,
                  dis_win_min=None, dis_win_max=None,
                  froz_win_min=None, froz_win_max=None,
                  num_iter=200, mix_ratio=0.5,
                  conv_tol=1e-10, conv_window=3,
                  device="cpu", dtype=torch.complex128,
                  check_every=1, verbose=False)
# res.M0_reduced and res.U0 go straight to wannierise()
```

The convenience entry `disentangle_case(case)` reads the windows and the iteration
parameters from `.win`; items absent from `.win` take the wannier90 defaults:

| `.win` keyword | `disentangle()` argument | Default |
|---|---|---|
| `dis_win_min`, `dis_win_max` | `dis_win_min`, `dis_win_max` | minimum and maximum of all eigenvalues |
| `dis_froz_max` | `froz_win_max` | no frozen window |
| `dis_froz_min` | `froz_win_min` | `dis_win_min` when only `dis_froz_max` is given |
| `dis_num_iter` | `num_iter` | 200 |
| `dis_mix_ratio` | `mix_ratio` | 0.5 |
| `dis_conv_tol` | `conv_tol` | 1e-10 |
| `dis_conv_window` | `conv_window` | 3 |

The implementation follows the source of wannier90 3.1 `disentangle.F90` step by
step, not only the formulas of the literature:

1. Windows (corresponds to `dis_windows`): the outer window selects `ndimwin(k)`
   bands at each k. Reading `.eig` checks that the energies ascend, so the selected
   bands are contiguous. The `ndimfroz(k)` bands inside the frozen window must be
   kept whole. `ndimwin < nw` or `ndimfroz > nw` raises an error, as in wannier90.
2. Initial subspace (corresponds to `dis_project` and `dis_proj_froz`):
   `V0 = polar(A_w)`, the unitary factor of a batched SVD, equivalent to Löwdin
   orthonormalisation. With a frozen window, the non-frozen part takes the leading
   `nw − nf` eigenvectors of `Q_f P_s Q_f` (SMV Eq. 27), and the frozen directions
   are unit vectors.
3. Main loop (corresponds to `dis_extract`): a Jacobi-style iteration in which all k
   update simultaneously, the same semantics as wannier90 when it distributes the
   k-points over MPI.
   - `Z_in⁽¹⁾ = Z[V⁽⁰⁾]`; `Z_in⁽ⁱ⁾ = mix·Z_out⁽ⁱ⁻¹⁾ + (1 − mix)·Z_in⁽ⁱ⁻¹⁾`
   - `V⁽ⁱ⁾` = the frozen states plus the leading `nw − nf` eigenvectors of `Z_in⁽ⁱ⁾`
   - `Z_out⁽ⁱ⁾ = Z[V⁽ⁱ⁾]`
   - Here `Z[V](k) = Σ_b w_b M^{k,b} P_{k+b} M^{k,b†}` with `P = V V†`, then
     projected onto the non-frozen complement at the k-point.

   The mixing acts on the Z matrices, not on the projectors. `Z_out` is assembled
   after the subspace update, at the end of the iteration. Both points are the same
   as in wannier90.
4. Convergence criterion: each iteration computes two values of Ω_I.
   - `Ω_I(i−1)` (SMV Eq. 18) uses the subspaces of the previous iteration for the
     neighbours and equals the frozen contribution plus the sum of the selected
     eigenvalues.
   - `Ω_I(i)` (SMV Eq. 12) uses the new subspaces throughout.
   - `delta = Ω_I(i−1)/Ω_I(i) − 1`; the run is converged when `|delta| < conv_tol`
     for `conv_window` consecutive iterations.
5. Initial gauge (corresponds to `internal_find_u`): `U0 = polar(V† A_w)`, the
   Löwdin-projection initial gauge of wannier90.

### k-dependent dimensions in one batched eigh

`ndimwin` and `ndimfroz` vary with k; for example, nf of li2mno3 spin 1 lies between
284 and 286. A loop over k would break the batching, so the implementation keeps the
full nb dimension and replaces slicing by masks and diagonal level shifts:

- The frozen and the out-of-window rows and columns of Z are set to zero, then
  shifts are added on the diagonal: `+2·w_tot` for frozen directions and `−2·w_tot`
  for out-of-window directions. Because ‖M‖ ≤ 1 and w_b > 0, the spectrum of Z lies
  in `[0, w_tot]`, so the ascending eigenvalues are ordered
  `[out of window | free spectrum | frozen]`.
- The last nw columns of a single `torch.linalg.eigh` are the new subspace. The
  frozen block is exactly `2w_tot·I`, and its eigenvectors are a unitary combination
  of the frozen unit vectors. The subspace is exact, and the gauge within the
  columns has no physical effect.
- The sum of the selected eigenvalues needed by `Ω_I(i−1)` is vectorised with a
  per-k slot mask `nb − nw ≤ pos < nb − nf(k)`.
- The shifts are added to both `Z_in` and `Z_out`. The mixing is linear, so the
  shifts survive it and the slot selection holds at every iteration.
- After eigh, cross-block floating-point noise of about 1e−15 is cleared with masks
  by column type (`_scrub_columns`), so the frozen states stay strictly frozen.
- The `QPQ` of the initial subspace (spectrum in `[0, 1]`) is treated the same way,
  with shifts of ±2. The complement block decouples exactly from the frozen block,
  so the eigenvectors of near-zero eigenvalues are automatically orthogonal to the
  frozen states; this is the property that the ortho-fix of wannier90 enforces.

### Differences from wannier90

After disentanglement converges, wannier90 right-multiplies V by the eigenvectors of
the subspace Hamiltonian, for use in band interpolation. gpuw90 does no band
interpolation and omits this step.

The step changes only the gauge. The right-multiplication `V → VR` cancels exactly
between `M0_reduced = V† M0 V'` and `U0 = polar(V† A)`, because
`polar(R†B) = R†·polar(B)`. Ω_I, the centres, the spreads and the localisation
trajectory are unaffected. The entangled-band initial states agree with the
reference to the print precision of `.wout` (see "Verification results"), consistent
with this derivation.

## Numerical issues and their handling

Every item below has a comment at the relevant place in the code.

### The 1/N prefactor and the overall sign of the gradient

The literature often writes the MV gradient as `4 Σ_b w_b (A[R] − S[T])`. But Ω
itself carries `(1/N) Σ_k`, and an implementation also easily gets the overall sign
wrong.

gpuw90 tests the gradient with a finite-difference directional derivative
(`fd_check.py`): it compares the analytic `Re Tr[g† X]` with
`[Ω(U e^{εX}) − Ω(U e^{−εX})] / 2ε`.

| Form of the gradient | ratio on GaAs |
|---|---|
| `4 Σ_b w_b (A[R] − S[T])` | −1/8, i.e. −1/N (N = 8) |
| `−(1/N) · 4 Σ_b w_b (A[R] − S[T])` (current implementation) | +1.0000, the same for all three ε |

The true gradient is therefore `g = −(1/N) · 4 Σ_b w_b (A[R] − S[T])`, and Ω
decreases monotonically along the descent `dW = −eff·g`. After changing the
gradient, run `fd_check.py` first; do not take the prefactor and sign straight from
a formula.

### Splitting Ω_I and Ω_OD

Ω_I must sum over all m, n; Ω_OD takes only the off-diagonal part:

- `Ω_I = (1/N) Σ_{k,b} w_b (nw − Σ_{m,n} |M_mn|²)`
- `Ω_OD = (1/N) Σ_{k,b} w_b Σ_{m≠n} |M_mn|²`
- `Ω_D = (1/N) Σ_{k,b} w_b Σ_n (Im ln M_nn + b·r_n)²`

A common mistake is to put `Σ_n |M_nn|²` into Ω_I. Ω_OD is then first merged into
Ω_I and afterwards added once more on its own, so the total Ω is too large by one
Ω_OD. The symptom is correct centres (deviation 6e-7 Å) and a total spread too high
by 0.502, which equals Ω_OD. With the formulas above, the initial-state Ω of GaAs
agrees with the reference digit for digit (4.46881212).

### Branch cut of ln

`Im ln M_nn = atan2(Im, Re)` takes values in (−π, π], the same as
`aimag(log(...))` in wannier90. CPU and GPU take the same principal value, so the
branch does not jump from one device to the other.

For the WFs of an insulator, M_nn ≈ 1: the phase is small and far from ±π. A branch
flip is possible only near a tie or when |M_nn| → 0; it would move the centre by the
order of `b·2π/…`. Ω_D of the isolated-band check cases is small and no flip occurs.

### Shell weights of the b vectors

The weights are rebuilt by the kmesh rules of wannier90:

1. Group the b vectors into shells by |b|.
2. Solve by least squares the 6 independent symmetric components of
   `Σ_shell w_s (Σ_{b∈s} b_α b_β) = δ_αβ`.
3. On a uniform mesh every k has the same shell structure; the weights are filled
   back into each (k, b) by matching |b|.

Loading computes the B1 residual `max|Σ_b w_b b⊗b − I|`: 9e-18 for GaAs and 2e-16
for Diamond. w_b agrees with the reference digit for digit (GaAs 0.408635, Diamond
0.527876).

### Phases of degenerate eigenvectors (nothing to handle)

Within a degenerate subspace the eigenvectors of `i·dW` are fixed only up to a
unitary rotation or a phase. But `exp(dW) = V·f(λ)·V†` is a function of the operator
and does not depend on this gauge, so degeneracy does not affect the update. The
same holds for the Löwdin `(A†A)^{−1/2}`. The code notes this; no phase-fixing logic
is needed.

### FP64 throughout

All computation uses complex128 / float64, in line with the double precision of
wannier90.

CPU and GPU differ only in the order of reductions, by about 1e-14. When both
converge to the same minimum, the effect on the centres is far below 1e-4 Å; the
deviations measured in the checks are of the order of 1e-6 to 1e-9 Å (see
"Verification results"). gpuw90 does not require CPU and GPU to agree bit for bit,
only to be physically equivalent.

### WF order and periodic images

The order of the WFs from the minimisation may differ from the reference, and nearly
degenerate WFs may swap.

The check scripts therefore pair the WFs greedily by nearest centre, not by index.
Distances use the minimum image: fold into the real-space lattice first, then take
the length. This tolerates degenerate swaps and lattice translations. The threshold
is 1e-4 Å.

### Column-major order of .mmn

Within each (k, b) block the nb×nb data have the band index m as the fastest-varying
index, i.e. Fortran column-major order. Reading uses `reshape(order="F")` to obtain
`M[m,n]`.

If the data are misread as row-major, the symmetry of the centres is broken while
the B1 residual stays normal. Comparing the initial-state centres with the reference
exposes this error directly.

### Disentanglement: order of the two Ω_I evaluations

The order must be:

1. the frozen contribution to `Ω_I(i−1)` uses the neighbour subspaces of the
   previous iteration and is computed as a whole before any k is updated;
2. `Ω_I(i)` uses the new subspaces;
3. `Z_out` is assembled after the subspace update, at the end of the iteration.

If any of these is out of order, the delta sequence departs from wannier90 and the
number of iterations to convergence no longer matches. With this order, the
iteration counts of the four entangled-band channels equal those of wannier90 one by
one.

### Disentanglement: size of the level shift

If the diagonal shifts of the frozen and out-of-window directions are very large
(e.g. 1e6), the error of the eigh eigenvectors (about eps·‖A‖/gap) is amplified and
mixes into the free subspace.

gpuw90 uses `±2·w_tot`, of the same order as the spectrum. Because ‖M‖ ≤ 1, the
spectrum of Z is bounded by `[0, w_tot]`, so this shift is enough to separate the
blocks and loses no precision. The `QPQ` of the initial subspace (spectrum in
`[0, 1]`) likewise takes ±2.

## Verification results

The 2 isolated-band cases and the 4 entangled-band channels all PASS. Every maximum
centre deviation is below 1e-6 Å; the check threshold is 1e-4 Å. Commands and
arguments are in "Running the checks" of README.md.

### Isolated bands

The reference output is the `benchmark.out` of the wannier90 3.1 test-suite, fully
converged: the GaAs centres no longer change from iteration 2 on.

| Case | nb/nw | nk | nntot | Iterations to convergence | Max centre deviation (Å) | Max spread deviation (Ų) | Ω_total deviation |
|---|---|---|---|---|---|---|---|
| GaAs (example01) | 4/4 | 8 | 8 | 10 | 8.31e-07 | 4.07e-09 | 2.6e-10 |
| Diamond (example05) | 4/4 | 64 | 8 | 12 | 3.13e-09 | 1.98e-09 | 1.3e-10 |

- The threshold is 1e-4 Å; the margin of the two cases is a factor of 10⁵–10⁷.
  Diamond has 64 k-points and exercises the batched-k path.
- The initial-state centres (after Löwdin, before optimisation) differ from the
  reference "Initial State" by 6e-7 Å (GaAs) and 3e-9 Å (Diamond). This item
  involves no minimisation and verifies on its own the file reading, the shell
  weights and the centre formula.

Reproduce: `python3 check_against_w90.py testdata/gaas --seed gaas`; for Diamond use
`testdata/diamond --seed diamond`.

### Entangled bands

The reference output comes from serial runs of wannier90.x 3.1 on the SAI cluster
(`run_np1_s*/wannier90.*.wout`). The reference runs and gpuw90 read the same VASP
export (`.mmn/.amn/.eig/.win`, symbolic links to one copy). Both cases have a frozen
window (`dis_froz_max` from `.win`); the outer window is the default full range.
li2mno3 has the same system size as the LEOS production calculations (288 WFs).

| Case (spin) | nb/nw | nk/nnb | nf | Disentanglement iterations (gpuw90 / reference) | Ω_I (Ų, same printed value on both sides) | \|ΔΩ_I\| | Max initial-state centre deviation (Å) | Max initial-state spread deviation (Ų) |
|---|---|---|---|---|---|---|---|---|
| vircef s1 | 180/142 | 8/12 | 141 | 82 / 82, converged | 81.06033903 | 3.4e-09 | 7.9e-07 | 5.0e-09 |
| vircef s2 | 180/142 | 8/12 | 140 | 188 / 188, converged | 80.32720233 | 4.6e-10 | 7.7e-07 | 4.9e-09 |
| li2mno3 s1 | 356/288 | 8/8 | 284–286 | 38 / 38, converged | 154.86511542 | 4.4e-09 | 8.3e-07 | 5.0e-09 |
| li2mno3 s2 | 356/288 | 8/8 | 240 | 200 / 200, both at the `dis_num_iter` limit, not converged | 135.30457827 | 4.1e-09 | 8.2e-07 | 5.0e-09 |

Gates: |ΔΩ_I| ≤ 1e-4 Ų and a maximum initial-state centre deviation ≤ 1e-4 Å. All
four channels PASS.

Why these quantities are checked:

- The localisation part of the reference `.wout` uses CG with a fixed 100 iterations
  and no convergence criterion (`conv_window = −1`). None of the four reference
  final states is converged; ΔΩ in the last iterations is still 1e-2 to
  1e-1 Ų/iteration.
- The difference between the final-state centres therefore reflects the different
  run settings of the two sides. It measures about 0.3 Å, and the converged Ω_total
  of gpuw90 is lower (vircef s1: 94.2618, reference 94.6541). It cannot tell whether
  the implementation is correct.
- The hard check therefore uses only quantities that do not depend on the run
  settings: Ω_I (the gauge-invariant output of disentanglement), and the centres and
  spreads of the "Initial State" after disentanglement. The latter is the subspace
  plus the Löwdin-projection initial gauge, a deterministic computation with no
  minimisation.
- Before choosing the quantities to check, read the run parameters echoed in the
  reference `.wout` and confirm whether the reference itself is converged.

Reading the results:

- The initial-state centre deviation of about 8e-7 Å and spread deviation of about
  5e-9 Ų correspond to the rounding at the print precision of `.wout` (6 decimals
  for centres, 8 decimals for spreads). The 8.3e-7 Å of isolated-band GaAs has the
  same origin. The actual agreement is better than this.
- The disentanglement iteration counts of the four channels equal those of wannier90
  one by one (38, 82, 188, 200). This shows that the mixing, the order of the two
  Ω_I evaluations and the semantics of the convergence window all match wannier90.
- li2mno3 s2 is truncated at 200 iterations, and its unconverged trajectory is still
  reproduced to 4e-9 Ų, so the floating-point drift of the fixed-point iteration is
  negligible.
- The script checks `|Ω_I(initial state) − Ω_I(disentangled)|` on every run. Ω_I is
  gauge-invariant; the measured difference is about 1e-11.

Check on CUDA: li2mno3 also passes on a V100 (|ΔΩ_I| = 4.4e-9 Ų, initial-state
centre deviation 8.3e-7 Å); the measurement conditions are in "Measured performance"
of README.md.

Reproduce: `python3 check_dis_against_w90.py testdata/vircef --seed wannier90.1`;
for the other channels change the directory and `--seed`. The large inputs are not
committed; run `testdata/fetch_entangled.sh` first.

## GPU memory and tiling

The working set is small at present and GPU memory is not a constraint. The code
does no tiling; all tensors stay on the device.

- The GPU feasibility assessment estimates the whole working set of disentanglement
  and localisation at < 2 GB. Source: the section on the scale profile and baseline
  breakdown of `docs/reports/20260710_gpu_feasibility.md`; the file has been moved
  out of the working tree and can be found in the git history.
- Estimate for a Ce system (nk = 32, nw ≈ 253, nntot ≈ 8–12):
  `M0 = 32 × 10 × 253² × 16 B ≈ 0.33 GB`, `U ≈ 0.1 GB`. A single V100-SXM2-16GB card
  is enough.
- The disentanglement stage keeps the window-masked full-band overlaps
  `M0w [nk,nnb,nb,nb]` for the whole iteration.

For a larger system that needs tiling, the split can be made at two levels:

| Level | Split | Basis |
|---|---|---|
| Spin | Process the two channels of ISPIN = 2 one after the other, or in parallel on two GPUs | The two spin channels are fully independent. gpuw90 itself processes one spin channel at a time, per seed |
| k-points | Stream `M0` from the host in k tiles; U stays whole on the device, and only `M0[k_tile]` and its neighbour tiles are swapped in | Rotating M needs the neighbour `U[k+b]`, and U is small. The sum over b in the gradient and the reduction for the centres involve only the local k; the only cross-tile access is the lookup of the neighbour U |

k tiling is needed only for very large systems.

The largest demand on GPU memory is the `nb × Nr ≈ 5.3 GB` matrix of SCDM, but SCDM
is outside the scope of gpuw90: the initial guess is read from the `.amn` that VASP
writes.

## Connection to the LEOS pipeline

gpuw90 can do the disentanglement and localisation steps of the Wannierisation
offline, outside VASP, in place of the same two steps of the wannier90 library mode
inside VASP.

| Part | Current state |
|---|---|
| Input | VASP with `LWRITE_MMN_AMN` set (an existing switch in mlwf.F) writes `seedname.mmn`, `.amn`, `.eig`; a `.win` that provides the cell and explicit kpoints is also needed, written by VASP or by hand. gpuw90 reads these files offline, does not enter the VASP process, and manages its own GPU memory |
| Output | `centres [nw,3]` (Å) and `spreads [nw]` of each spin channel, returned as numpy arrays in `WannieriseResult` |
| Warm start between frames | `wannierise()` accepts any unitary `U0`, so the converged `res.U` of the previous frame can be passed as the initial gauge. The repository has no test of this use |

There is no adapter yet that writes the output as EOSCAR or hands it to pyleos; see
"Known limitations".

## Known limitations

- Localisation has only fixed-step steepest descent, with no CG and no line search.
  Well-conditioned systems (GaAs, Diamond) converge in 10–12 iterations; on nearly
  degenerate systems such as vircef and li2mno3, SD has not reached 1e-10 after 2000
  iterations. The CG of wannier90.x on the cluster has not converged after 100
  iterations either.
- With a near-degeneracy at the edge of the frozen window, the subspace selection is
  a discrete decision that is sensitive to the order of floating-point operations,
  and the downstream centres may move as a result. The only safeguard at present is
  the regression against the wannier90 reference output
  (`check_dis_against_w90.py`).
- The check scripts pair the WFs greedily by nearest centre
  (`check_against_w90.match_wfs`), not by a globally optimal matching. Mispairing is
  possible when nw is very large.
- When `.win` has no explicit kpoints block, `load_case()` must be given
  `nnkp_path`; it does not generate the k-points from `mp_grid`.
- The `guiding_centres` option of wannier90 is not supported.
- There is no output adapter: results are not written as EOSCAR, and there is no
  file format for pyleos to read.
- Every localisation iteration re-rotates M from `M0` as a whole; there is no
  incremental update as in wannier90.
- The entangled-band check has only the two cases vircef and li2mno3, and no real
  frame of a Ce system.
