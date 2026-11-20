# v2 → v3 incremental material (archive)

This directory holds the incremental material of the v3 development relative to v2
(patch/leos.patch). All of it is merged into patch/leos_v3.patch. To install, use
leos_v3.patch directly; see [../README.md](../README.md).

v3 was developed and verified on the production tree `vasp.6.4.3_leos_v2` on the SAI
cluster (build job 593892, verification jobs 593893 and 593894, 2026-07-07).

## Files

| File | Content |
|---|---|
| eos.F | Full replacement of the v2 eos.F |
| apply_v3_mlwf.py | Injects the LEOS variables into mlwf.F, and the handling of SCDM-only and guiding centres in the section that writes .win |
| apply_v3_main.py | Injects the trajectory-mode wiring into main.F: wannierisation and analysis step by step under EOS_DO_THIS_STEP (still only at the first step when only LWANNIER90 is set), and removes the unconditional MLWF_FREE at every step |

Features of eos.F that are new relative to v2:

- `EOS_PROJ_WEIGHTS`: the projection norm of every WF on each atom (#WP block), used
  by projection-assisted arbitration. Weighting by the QION augmentation difference
  was tried first; its contrast in soft-pseudopotential channels is only 1.35:1,
  against about 5:1 with the projection norm, so the former was dropped.
- `EOS_SCDM_ONLY`: skips the iterative localisation of wannier90. On CO2 the
  oxidation states are unchanged and the wannier90 time is 2 s → 0.023 s.
- `EOS_WARM_START` (on by default, effective when EOS_EVERY_N > 0): passes the WF
  centres of the previous step to wannier90 as guiding centres.

## Observations during development

- CO2 at 92 K, 3 MD steps: the greedy result drifts at step 3 (C goes from +4 to +3
  or +2).
- On small systems the guiding centres cost more than they gain (0.58 s → 1.3 s).
- Partial occupations in MD increase num_wann (32 → 40); pyleos excludes the surplus
  WFs with the cap on the NELECT total.
