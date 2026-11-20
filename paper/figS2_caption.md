# Fig. 4 - caption draft

**Fig. 4** Diagnostics that indicate when an assignment can be trusted.

**(a)** Reliability index *R*<sub>elem</sub> of the 72 scored benchmark cases,
grouped by whether the element-mean oxidation states agree with the reference
(deviation ≤ 0.25; *n* = 64) or deviate from it (*n* = 8). Each symbol is one
case; crosses mark cases with at least one stray Wannier function (WF; spread
> 3 Å² or centre farther than 1.6 Å from every atom). Black bars are the
medians (83.0 and 25.55) and grey boxes the interquartile ranges. All
deviating cases are named, as are the two agreeing cases that lie below the
median of the deviating group. Asterisks mark the two deviating cases with a
high index, for which the reference follows another convention than the
electronegativity rule (Ca<sub>3</sub>AuN: Au<sup>3−</sup> versus the electride
description; WAu<sub>12</sub>: 18-electron rule). "(legacy)" denotes an entry
of the historical data set that shares its formula with a newer entry. Every
deviating case has either *R*<sub>elem</sub> ≤ 25.9 or at least one stray WF;
the converse does not hold, since 16 of the 64 agreeing cases contain stray WFs
and two have *R*<sub>elem</sub> < 25.55.

**(b)** Number of the 72 scored cases that agree with the reference as a
function of the electronegativity weight *β* in the cost
*d*² + *β*(1 − *χ*/*χ*<sub>max</sub>). The dashed line marks the default,
*β* = 10 Å². The count is 62 for *β* = 9–12, 63 at *β* = 8, 58 and 55 at
*β* = 5 and 3, and 61 at *β* = 15 and 20. This scan was carried out before the
tungsten chloride structures of the benchmark were corrected and therefore
gives 62 at *β* = 10, where the benchmark of panel (a) gives 64.

**(c)** Convergence of the Wannierisation for CoO, FeO, MnO, NiO and
ReO<sub>3</sub> (↑, ↓: spin components 1 and 2). Plotted is the excess total
spread Ω − Ω<sub>min</sub>, where Ω<sub>min</sub> is the lowest total spread
obtained for that spin channel in any run (CoO 8.020/8.020, FeO 5.284/7.352,
MnO 5.813/7.685, NiO 7.620/7.620, ReO<sub>3</sub> 10.435/10.434 Å² for ↑/↓).
Diamonds: archived wannier90 run, stopped after the default 100 iterations
without a convergence test. Plus signs: gpuw90 minimisation started from the
same gauge. Triangles: gpuw90 minimisations started from 16 random gauges,
ordered by Ω. Open plus signs and triangles did not meet the convergence
criterion (|ΔΩ| < 10⁻¹⁰ Å² over three iterations) within 2 × 10⁵ iterations.
Numerals give the number of stray WFs of a solution (no numeral: none); the
totals over the ten channels are 13 (archived), 2 (same starting gauge) and 0
(lowest random start), and none of the 160 random-start end points contains a
stray WF. In the shaded channels (†; FeO ↓ and MnO ↓) no run converged, so
Ω<sub>min</sub> is the lowest value reached rather than a minimum. The
vertical axis is linear below 10⁻³ Å² (dotted line) and logarithmic above.

---

## Provenance

Produced by `paper/figS2_reliability.py` (run `python3 figS2_reliability.py`
in that directory; it prints every plotted number with its source file).

| Panel | File | What is read |
|---|---|---|
| (a) | `results/bench74_golden.tsv` | columns `batch`, `case`, `status`, `R_elem`, `n_stray`; the two `run` rows (no reference) are excluded; "deviates" = status `~` or `FAIL` |
| (b) | `results/scan_beta.tsv` (from `tools/beta_scan.py`) | columns `beta`, `n_OK`; 72 = `n_OK + n_mid + n_FAIL`; default β from the column name `changed_vs_beta10`, cross-checked against the default of `pyleos.transport.assign_case` |
| (c) | `<holdout>/a1/multistart_summary.json`, `multistart_summary.tsv`, `compare_standard.json` | total spread Ω of the archived run, of the same-start gpuw90 run and of the 16 random starts; convergence flags; stray counts. The three files are checked against each other at run time. Protocol details (100 iterations, convergence criterion, iteration cap) and the judgement on FeO ↓ / MnO ↓ are from `reports/restore.json`, `reports/multistart.json` and `review/verdict.json` |

Numbers checked against `docs/RESULTS.md`: 64 agreeing / 8 deviating cases
(section 2) and medians 83.0 / 25.55 (section 8); 62 at β = 10 in the scan
(section 7). The stray totals 13 / 2 / 0 are not in `RESULTS.md`; they agree
with `reports/multistart.json` and `review/verdict.json`. No mismatch.

## Caveats

- **Panel (c) is provisional.** Its data are held-out test inputs that are not
  yet under `cases/`; the script reads them from one constant (`HOLDOUT_A1`,
  a session scratch directory) and skips the panel with a printed notice if
  that directory is absent. Only Wannier-quality numbers were read from it;
  no oxidation state was computed, read or reported, and nothing was copied
  into the repository. Before publication the runs need a permanent place in
  `cases/` and the constant must be repointed.
- **Panel (c), convergence.** FeO ↓ and MnO ↓ did not converge from any start
  (16 distinct end points each, 0.38 and 0.46 Å² apart at most); the reviewer
  classifies them as stall points near a singular point of the spread
  functional, usable only for the statement that no stray WF appears. In
  FeO ↑ and MnO ↑ the same-start run and 12 and 13 of the 16 random starts
  stop at the iteration cap just short of the criterion, within
  1.3 × 10⁻⁴ Å² of a minimum that converged starts do reach (FeO ↑ has two
  such minima 0.0064 Å² apart; the same-start run and 5 random starts end in
  the higher one). One random start of NiO ↑ ends 4.06 Å² above the minimum,
  unconverged. The same-start gpuw90 runs of CoO ↑ and NiO ↑ converge to the
  higher minimum of the archived run and keep its stray WF.
- **Panel (b)** uses the stored scan, which predates the corrected tungsten
  chloride structures (62 instead of 64 at β = 10). A check run of
  `tools/beta_scan.py` on the current tree (output kept outside the
  repository, not plotted) gives 57, 60, 65, 64, 64, 64, 64, 63, 63 for
  β = 3, 5, 8, 9, 10, 11, 12, 15, 20: the same curve shifted by +2, with the
  same cases changing status. Regenerating the stored file would therefore
  make panels (a) and (b) refer to the same benchmark without changing the
  shape of the curve; the caption sentence on 62 versus 64 could then be
  dropped and the counts in the caption updated.
- **Panel (a), labels.** `As4S4` and `CSe2` each occur twice in the benchmark
  (batches `v1_mol` and `wave1`/`wave2`); the `v1_mol` entries are labelled
  "(legacy)". The asterisks are assigned in the script from the discussion in
  `docs/RESULTS.md` (sections 3 and 8), not from a data column.
- **No ghost panel.** `results/ghost_ab_summary.tsv` does not
  support the planned panel (d) cleanly: (i) it holds no element means of the
  runs before the change, only `old_status` and `old_dev`; (ii) the deviation
  recomputed from its `flow_elem_avg_v31` column differs from `old_dev` for two
  systems (wcl4 2.00 → 1.00, wcl5 3.50 → 2.50), so the statement of
  `docs/RESULTS.md` section 8 that the element means are unchanged holds for
  the verdicts of the 12 scored systems but not for the means of these two
  (recomputing the earlier cells of commit a778265 with pyleos gives W +2.00
  → +3.00 and +1.50 → +2.50); (iii) its four tungsten
  chloride rows are the truncated cells that were later replaced in the
  benchmark; (iv) the ghost count changes for only 3 of 13 systems
  (Fe<sub>3</sub>O<sub>4</sub> 42 → 20, wcl3 5 → 1, wcl5 7 → 5). The script
  prints this check at every run.
