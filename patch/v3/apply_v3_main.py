#!/usr/bin/env python3
"""Apply the LEOS v3 main.F trajectory-mode rewiring (idempotent).

Usage: python3 apply_v3_main.py <path-to-src>/main.F
Handles both the pre-v3 server integration (first-step-only calls, tab
indentation) and reports if the hook is already the v3 form.
"""
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()

if "EOS_DO_THIS_STEP" in s:
    print("main.F already patched - nothing to do")
    sys.exit(0)

old = """      IF ((.NOT.L2E4W) .AND. MLWFRUN()) THEN
\t         IF (IO%LORBIT==14) THEN
\t            CALL SPHPRO_OPTPROJ(P,T_INFO,W,WDES,LDIMP,LMDIMP,-1,KPOINTS%EMIN,KPOINTS%EMAX)
\t         ENDIF

\t         CALL MLWF_FREE(MLWF_GLOBAL)
\t         IF (NSTEP == 1) CALL EOS_PREPARE_MLWF(W, INFO, IO)
\t         IF (NSTEP == 1) CALL MLWF_MAIN(WDES,W,P,KPOINTS,CQIJ,T_INFO,LATT_CUR,INFO,IO,MLWF_GLOBAL)
\t         IF (NSTEP == 1) CALL EOS_ANALYZE(MLWF_GLOBAL, W, T_INFO, LATT_CUR, INFO, P, IO, DYN)"""

new = """      IF ((.NOT.L2E4W) .AND. (MLWFRUN() .OR. LEOS)) THEN
\t         IF (IO%LORBIT==14) THEN
\t            CALL SPHPRO_OPTPROJ(P,T_INFO,W,WDES,LDIMP,LMDIMP,-1,KPOINTS%EMIN,KPOINTS%EMAX)
\t         ENDIF

\t         IF (NSTEP == 1) CALL EOS_PREPARE_MLWF(W, INFO, IO)
\t         ! LEOS trajectory mode: (re)wannierise and analyse on every step
\t         ! selected by EOS_EVERY_N; plain LWANNIER90 users keep the old
\t         ! first-step-only behaviour.
\t         IF ((MLWFRUN() .AND. NSTEP == 1) .OR. EOS_DO_THIS_STEP(NSTEP, DYN%NSW, INFO%LSTOP)) THEN
\t            CALL MLWF_FREE(MLWF_GLOBAL)
\t            CALL MLWF_MAIN(WDES,W,P,KPOINTS,CQIJ,T_INFO,LATT_CUR,INFO,IO,MLWF_GLOBAL)
\t         ENDIF
\t         IF (EOS_DO_THIS_STEP(NSTEP, DYN%NSW, INFO%LSTOP)) THEN
\t            CALL EOS_ANALYZE(MLWF_GLOBAL, W, T_INFO, LATT_CUR, INFO, P, IO, DYN)
\t         ENDIF"""

assert s.count(old) == 1, "hook anchor not found (unexpected main.F variant)"
s = s.replace(old, new)
p.write_text(s)
print("main.F patched")
