#!/usr/bin/env python3
"""Apply the LEOS v3 mlwf.F injections (idempotent).

Usage: python3 apply_v3_mlwf.py <path-to-src>/mlwf.F
Adds (1) the LEOS hook variables and (2) the win-writer injection for
SCDM-only mode and trajectory warm-start guiding centres.
"""
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()

if "LEOS_SCDM_ONLY" in s:
    print("mlwf.F already patched - nothing to do")
    sys.exit(0)

old = "   !> Write the LOCPROJ functions as a wavecar file"
new = """   !> LEOS hooks (set by the eos module; mlwf cannot USE eos - circular).
   !> SCDM-only skips the wannier90 minimisation; guide centres implement
   !> the trajectory warm start via wannier90 guiding_centres.
   LOGICAL, PUBLIC, SAVE :: LEOS_SCDM_ONLY = .FALSE.
   LOGICAL, PUBLIC, SAVE :: LEOS_GUIDE = .FALSE.
   REAL(q), ALLOCATABLE, PUBLIC, SAVE :: EOS_GUIDE_CENTRES(:,:,:)

   !> Write the LOCPROJ functions as a wavecar file"""
assert s.count(old) == 1, "module anchor not found"
s = s.replace(old, new)

old = """    ! write unit cell
    WRITE(99,'(A)') 'begin unit_cell_cart'"""
new = """    ! LEOS: SCDM-only fast mode - skip disentanglement and localisation
    IF (LEOS_SCDM_ONLY) THEN
       WRITE(99,'(A)') 'num_iter = 0'
       WRITE(99,'(A)') 'dis_num_iter = 0'
    ENDIF

    ! LEOS trajectory warm start: guide this step with the previous step
    ! centres. The projections block only steers the guiding centres; the
    ! initial U matrix still comes from the VASP-side amn (SCDM).
    IF (LEOS_GUIDE .AND. ALLOCATED(EOS_GUIDE_CENTRES)) THEN
       NI = 1
       IF (INDEX(MLWF%SEED_NAME, '.2') > 0) NI = 2
       IF (SIZE(EOS_GUIDE_CENTRES, 2) >= MLWF%NUM_WANN) THEN
          WRITE(99,'(A)') 'guiding_centres = .true.'
          WRITE(99,'(A)') 'begin projections'
          DO NK = 1, MLWF%NUM_WANN
             WRITE(99,'(A,F12.6,A,F12.6,A,F12.6,A)') 'c=', &
                EOS_GUIDE_CENTRES(1,NK,NI), ',', EOS_GUIDE_CENTRES(2,NK,NI), &
                ',', EOS_GUIDE_CENTRES(3,NK,NI), ':s'
          ENDDO
          WRITE(99,'(A)') 'end projections'
       ENDIF
    ENDIF

    ! write unit cell
    WRITE(99,'(A)') 'begin unit_cell_cart'"""
assert s.count(old) == 1, "win-writer anchor not found"
s = s.replace(old, new)

p.write_text(s)
print("mlwf.F patched")
