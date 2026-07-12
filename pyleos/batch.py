#!/usr/bin/env python3
"""LEOS batch driver for sai-pod. Self-contained (stdlib only) so it can be
scp'd to the cluster and run there directly.

Workflow:
  python3 batch.py stage   --dest ~/leos_batch1 --cases cases.json
  python3 batch.py submit  --dest ~/leos_batch1 [--max-queue 8]
  python3 batch.py status  --dest ~/leos_batch1
  python3 batch.py harvest --dest ~/leos_batch1 --out ~/leos_batch1_results

cases.json: [{"name": "co2", "src": "/path/with/POSCAR+POTCAR+KPOINTS",
              "incar": {"ENCUT": 500}}, ...]
Each case directory gets the LEOS INCAR template (per-case overrides merged),
a 4-GPU sbatch, and a state entry in <dest>/batch_state.json.

Failure triage (status/harvest):
  FAIL_DFT   - SCF crashed (EDDDAV/ZHEGV/"I REFUSE") or OUTCAR incomplete
  FAIL_W90   - DFT fine but no wannier90 Final State
  FAIL_LEOS  - Wannierisation fine but no usable EOSCAR: file missing, or it
               carries the LEOS_ABORTED marker (the NELECT-abort path writes
               a diagnostics-only EOSCAR, which is not a result)
  OK         - EOSCAR present without LEOS_ABORTED
Retry ladders (retry command): DFT -> ALGO=Normal then ALGO=All;
W90 -> EOS_SCDM_ONLY=.FALSE. + larger NBANDS margin.

Tiered screening (stage --mode tiered + escalate):
  stage assigns SCDM-only fast mode to every case EXCEPT those containing
  4d/5d/f-block elements (d/f manifolds are sensitive to the initial guess -
  the WCl4/5/6 + MoCl2 lesson), which start on full MLWF directly.
  After runs finish, `escalate` re-queues any SCDM-only case whose EOSCAR
  carries quality flags (stray WFs, |dNe| drift, ghost excess) as full MLWF.

Offline-w90 pipeline (stage --offline-w90 + w90; issue 16):
  The GPU job only exports (LWRITE_MMN_AMN=T, LWANNIER90_RUN=F); the w90
  solve moves off the GPU onto the DSPRHBM CPU partition.  After `submit`,
  repeated `w90` calls pump every offline case through:
    SUBMITTED -> EXPORTED       export finished (OUTCAR complete + .mmn)
    EXPORTED  -> W90_SUBMITTED  per-spin .win generated (offline_w90_prep;
                                a trimmed spin also gets a column-sliced
                                .amn), s1/s2 run dirs staged, one CPU job
                                (both spins concurrent) submitted
    W90_SUBMITTED -> OK         both .wout carry a Final State and are
                                copied back into the case dir, which is then
                                pyleos load_case-ready (wout channel)
  Export-only runs write no usable EOSCAR (no #WF block): offline triage
  keys on the export files and the .wout, never on EOSCAR.
"""
from __future__ import annotations

import argparse
import getpass
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

VASP_BIN = "/home/sjtu-gongxueqing/zhangxianyang/apps/vasp.6.4.3_leos_v2/bin/vasp_std"
# The export switches need the v3.1 EOS_READER (respects an explicit
# LWANNIER90_RUN=F): offline-w90 staging therefore uses the current
# production binary, not the v2 one above.
OFFLINE_VASP_BIN = ("/home/sjtu-gongxueqing/zhangxianyang/apps/leos_current/"
                    "bin/vasp_std")
W90_BIN = ("/home/sjtu-gongxueqing/zhangxianyang/LEOS/w90_offline/build/"
           "wannier90-3.1.0-impi/wannier90.x")
OFFLINE_INCAR_SWITCHES = {"LWRITE_MMN_AMN": ".TRUE.",
                          "LWANNIER90_RUN": ".FALSE."}

INCAR_TEMPLATE = {
    "SYSTEM": "LEOS batch case",
    "PREC": "Accurate",
    "ENCUT": 500,
    "EDIFF": 1e-6,
    "ISMEAR": 0,
    "SIGMA": 0.05,
    "ISYM": 0,
    "NSW": 0,
    "IBRION": -1,
    "LWAVE": ".FALSE.",
    "LCHARG": ".FALSE.",
    "LEOS": ".TRUE.",
    "LSCDM": ".TRUE.",
    "LWANNIER90_AUTO_WINDOW": ".TRUE.",
    "EOS_WF_DETAIL": ".TRUE.",
    "EOS_PROJ_WEIGHT": ".TRUE.",
}

SBATCH_TEMPLATE = """#!/bin/bash
#SBATCH --job-name=leos-{name}
#SBATCH --partition=16V100
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --gpus-per-node=4
#SBATCH --qos={qos}
#SBATCH --time={walltime}
#SBATCH --output=print-out
#SBATCH --error=print-out

set -euo pipefail
set +u
if [[ -f /opt/sai_config/mps_mapping.d/${{SLURM_JOB_PARTITION}}.bash ]]; then
  source /opt/sai_config/mps_mapping.d/${{SLURM_JOB_PARTITION}}.bash
fi
source /etc/profile.d/lmod.sh
module purge
module load fftw/3.3.10
module load nvhpc/24.3-pgi-tuned
module load openmpi/5.0.8-nvhpc24.3-pgi-auto
set -u

mpirun -np ${{SLURM_NTASKS}} -map-by "${{MAP_OPT:-ppr:1:gpu}}" -mca coll_hcoll_enable 0 {vasp}
echo LEOS_BATCH_DONE
"""

# Offline wannier90 solve: one mpirun per spin, both spins concurrent, each
# pinned to a disjoint slice of the actual cgroup CPU set (the job never owns
# cores 0..N-1 on DSPRHBM).  Intel stack only: wannier90 built with
# gnu+OpenMPI5 segfaults on any multi-rank run (w90_offline/README.md).
W90_SBATCH_TEMPLATE = """#!/bin/bash
#SBATCH --job-name=w90-{name}
#SBATCH --partition={partition}
#SBATCH --qos={qos}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --time={walltime}
#SBATCH --output=w90-out
#SBATCH --error=w90-out

set -uo pipefail
set +u
source /etc/profile.d/lmod.sh
module purge
module load intel/compiler-rt/2024.2.1 intel/tbb/2021.13 intel/ifort/2024.2.1 \\
            intel/mpi/2021.13 intel/mkl/2024.2
set -u
export I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=off
AFF=($(python3 -c "import os;print(' '.join(map(str,sorted(os.sched_getaffinity(0)))))"))
slice () {{ local out=("${{AFF[@]:$1:$2}}"); local IFS=,; echo "${{out[*]}}"; }}
NP={np}; NT={threads}; SPAN=$((NP*NT))
{run_block}
wait
echo W90_OFFLINE_DONE
"""

W90_RUN_SPIN = ('( cd s{s} && env OMP_NUM_THREADS=$NT MKL_NUM_THREADS=$NT '
                'MKL_DYNAMIC=FALSE taskset -c "$(slice {off} $SPAN)" '
                'mpirun -np $NP {w90} wannier90.{s} ) > s{s}/run.log 2>&1 &')
W90_RUN_SINGLE = ('( env OMP_NUM_THREADS=$NT MKL_NUM_THREADS=$NT '
                  'MKL_DYNAMIC=FALSE taskset -c "$(slice 0 $SPAN)" '
                  'mpirun -np $NP {w90} wannier90 ) > w90_run.log 2>&1 &')


# 4d/5d/f-block elements whose localised manifolds are sensitive to the SCDM
# initial guess (WCl4/5/6 + MoCl2 benchmark lesson): these skip the fast
# screening tier and start on full MLWF directly.
HEAVY_TM = {
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
    "Tm", "Yb", "Lu", "Th", "U",
}


def poscar_elements(path: Path) -> set:
    lines = path.read_text().splitlines()
    toks = lines[5].split()
    return set() if all(t.isdigit() for t in toks) else set(toks)


def eoscar_quality(path: Path) -> list:
    """Self-contained quality read of an EOSCAR (no pyleos import on the
    cluster): returns the list of flags that warrant escalating a fast
    SCDM-only screen to full MLWF."""
    flags = []
    n_stray = 0
    for line in path.read_text(errors="replace").splitlines():
        s = line.split()
        if line.startswith("#") and "|dNe|" in line:
            try:
                if abs(float(line.split("=")[1])) > 1e-3:
                    flags.append("dNe")
            except ValueError:
                pass
        elif line.startswith("#WF") and len(s) >= 8 and s[1].isdigit():
            # columns: #WF spin idx occ cx cy cz spread ...
            try:
                if float(s[7]) > 3.0:
                    n_stray += 1
            except ValueError:
                pass
    if n_stray:
        flags.append(f"stray:{n_stray}")
    return flags


def write_incar(path: Path, overrides: dict) -> None:
    tags = dict(INCAR_TEMPLATE)
    tags.update(overrides or {})
    with path.open("w") as fh:
        for k, v in tags.items():
            fh.write(f"{k} = {v}\n")


def load_state(dest: Path) -> dict:
    f = dest / "batch_state.json"
    return json.loads(f.read_text()) if f.is_file() else {}


def save_state(dest: Path, state: dict) -> None:
    (dest / "batch_state.json").write_text(json.dumps(state, indent=1))


def stage(dest: Path, cases_file: Path, qos: str, walltime: str,
          mode: str = "full", offline_w90: bool = False) -> None:
    cases = json.loads(cases_file.read_text())
    dest.mkdir(parents=True, exist_ok=True)
    state = load_state(dest)
    for c in cases:
        name = c["name"]
        d = dest / name
        d.mkdir(exist_ok=True)
        src = Path(c["src"]).expanduser()
        for f in ("POSCAR", "POTCAR", "KPOINTS"):
            if not (src / f).is_file():
                print(f"[{name}] missing {f} in {src} - skipped")
                break
            shutil.copy2(src / f, d / f)
        else:
            tier = mode
            if mode == "tiered":
                heavy = poscar_elements(d / "POSCAR") & HEAVY_TM
                tier = "full" if heavy else "scdm"
                if heavy:
                    print(f"[{name}] heavy TM {sorted(heavy)} -> full MLWF")
            overrides = dict(c.get("incar", {}))
            if tier == "scdm":
                overrides.setdefault("EOS_SCDM_ONLY", ".TRUE.")
            vasp = VASP_BIN
            if offline_w90:
                overrides.update(OFFLINE_INCAR_SWITCHES)
                vasp = OFFLINE_VASP_BIN
            write_incar(d / "INCAR", overrides)
            (d / "run.sbatch").write_text(SBATCH_TEMPLATE.format(
                name=name[:12], qos=qos, walltime=walltime, vasp=vasp))
            rec = state.setdefault(name, {"status": "STAGED", "jobid": None,
                                          "retries": 0, "tier": tier})
            if offline_w90:
                rec["w90"] = "offline"
            print(f"[{name}] staged ({tier}"
                  f"{', offline-w90' if offline_w90 else ''})")
    save_state(dest, state)


def escalate(dest: Path) -> None:
    """Re-queue SCDM-only screens whose EOSCAR carries quality flags as full
    MLWF runs (the second rung of the tiered strategy)."""
    state = load_state(dest)
    n = 0
    for name, rec in sorted(state.items()):
        if rec.get("tier") != "scdm" or rec.get("status") == "ESCALATED":
            continue
        d = dest / name
        eoscar = d / "EOSCAR"
        if not eoscar.is_file():
            continue                     # failures go through retry, not here
        flags = eoscar_quality(eoscar)
        if not flags:
            print(f"[{name}] clean screen - keep")
            continue
        incar = d / "INCAR"
        tags = {}
        for line in incar.read_text().splitlines():
            if "=" in line:
                k, v = line.split("=", 1)
                tags[k.strip()] = v.strip()
        tags["EOS_SCDM_ONLY"] = ".FALSE."
        with incar.open("w") as fh:
            for k, v in tags.items():
                fh.write(f"{k} = {v}\n")
        for f in ("EOSCAR", "print-out"):
            p = d / f
            if p.is_file():
                p.rename(d / (f + ".scdm"))   # keep the screen for comparison
        rec.update(status="RETRY", tier="full")
        rec["status_note"] = "escalated: " + ",".join(flags)
        n += 1
        print(f"[{name}] ESCALATE ({','.join(flags)}) -> full MLWF re-queued")
    save_state(dest, state)
    print(f"{n} case(s) escalated; run submit to launch them")


def squeue_ids() -> set | None:
    """Ids of this user's queued/running jobs, or None when the query fails.

    None means "queue state unknown": callers must treat it conservatively -
    change no case state, (re)submit nothing - because an empty answer from a
    failed query would make every still-running job look finished."""
    try:
        user = os.environ.get("USER") or getpass.getuser()
    except OSError:
        user = ""
    if not user:
        print("WARNING: could not determine the username - queue state unknown")
        return None
    try:
        r = subprocess.run(["squeue", "-u", user, "-h", "-o", "%A"],
                           capture_output=True, text=True)
    except FileNotFoundError:
        print("WARNING: squeue not found - queue state unknown")
        return None
    if r.returncode != 0:
        print(f"WARNING: squeue failed (rc={r.returncode}): "
              f"{r.stderr.strip()[:120]} - queue state unknown")
        return None
    return set(r.stdout.split())


def submit(dest: Path, max_queue: int) -> None:
    state = load_state(dest)
    in_queue = squeue_ids()
    if in_queue is None:
        print("queue query failed - submitting nothing (occupancy unknown)")
        return
    active = sum(1 for v in state.values() if v.get("jobid") in in_queue)
    for name, rec in state.items():
        if rec["status"] not in ("STAGED", "RETRY"):
            continue
        if active >= max_queue:
            print(f"queue cap {max_queue} reached; run submit again later")
            break
        r = subprocess.run(["sbatch", "run.sbatch"], cwd=dest / name,
                           capture_output=True, text=True)
        if r.returncode == 0:
            jobid = r.stdout.strip().split()[-1]
            rec.update(status="SUBMITTED", jobid=jobid)
            active += 1
            print(f"[{name}] submitted {jobid}")
        else:
            print(f"[{name}] sbatch failed: {r.stderr.strip()[:80]}")
    save_state(dest, state)


def triage(d: Path) -> str:
    eoscar = d / "EOSCAR"
    if eoscar.is_file():
        with eoscar.open(errors="replace") as fh:
            head = [fh.readline() for _ in range(5)]
        if any("LEOS_ABORTED" in ln for ln in head):
            # NELECT abort writes a diagnostics-only EOSCAR (no OS analysis):
            # a failure to retry, not a success
            return "FAIL_LEOS"
        return "OK"
    log = d / "print-out"
    text = log.read_text(errors="replace") if log.is_file() else ""
    if any(k in text for k in ("I REFUSE", "EDDDAV: Call to ZHEGV", "BAD TERMINATION")):
        return "FAIL_DFT"
    if "Not enough bvectors" in text:
        # wannier90 kmesh failure (skewed/body-centred small cells with an
        # unlucky grid) - a k-grid problem, not an SCF one
        return "FAIL_W90"
    wouts = list(d.glob("wannier90*.wout"))
    if not wouts or not any("Final State" in w.read_text(errors="replace") for w in wouts):
        outcar = d / "OUTCAR"
        if outcar.is_file() and "General timing" in outcar.read_text(errors="replace")[-4000:]:
            return "FAIL_W90"
        return "FAIL_DFT"
    return "FAIL_LEOS"


def offline_export_triage(d: Path) -> str:
    """Triage for an export-only run (no EOSCAR, no solve): did the DFT
    finish and did it leave the mmn/amn export behind?"""
    outcar = d / "OUTCAR"
    if not (outcar.is_file()
            and "General timing" in outcar.read_text(errors="replace")[-4000:]):
        return "FAIL_DFT"
    if not ((d / "wannier90.1.mmn").is_file() or (d / "wannier90.mmn").is_file()):
        return "FAIL_EXPORT"
    return "EXPORTED"


def status(dest: Path) -> None:
    state = load_state(dest)
    in_queue = squeue_ids()
    if in_queue is None:
        print("queue query failed - showing stored statuses, re-triaging nothing")
    counts: dict = {}
    for name, rec in sorted(state.items()):
        offline = rec.get("w90") == "offline"
        if in_queue is None:
            st = rec["status"]          # unknown queue: change no state
        elif rec.get("jobid") in in_queue:
            st = "QUEUED/RUNNING"
        elif rec.get("w90_jobid") in in_queue:
            st = "W90_QUEUED/RUNNING"
        elif rec["status"] in ("SUBMITTED", "DONE", "RETRY_DONE"):
            st = offline_export_triage(dest / name) if offline \
                else triage(dest / name)
            rec["status"] = st
        else:
            st = rec["status"]
        counts[st] = counts.get(st, 0) + 1
        print(f"{name:<28} {st}")
    if in_queue is not None:
        save_state(dest, state)
    print("\nsummary:", "  ".join(f"{k}:{v}" for k, v in sorted(counts.items())))


def _load_prep():
    """offline_w90_prep module: next to this file (cluster deploy) or in the
    repo's script/ directory (local checkout)."""
    try:
        import offline_w90_prep
        return offline_w90_prep
    except ImportError:
        pass
    script_dir = Path(__file__).resolve().parent.parent / "script"
    if (script_dir / "offline_w90_prep.py").is_file():
        sys.path.insert(0, str(script_dir))
        import offline_w90_prep
        return offline_w90_prep
    raise SystemExit("offline_w90_prep.py not found - deploy it next to "
                     "batch.py or keep the repo script/ directory")


def _wout_done(p: Path) -> bool:
    return p.is_file() and "Final State" in p.read_text(errors="replace")


def _stage_spin_dirs(d: Path, summary: dict) -> list:
    """s<N>/ run dir per spin: symlinked mmn/eig, per-spin win copy, and the
    trimmed amn for a narrowed spin (symlink to the shared one otherwise).
    Returns the staged spin numbers."""
    spins = []
    for rec in summary["spins"]:
        s = rec["spin"]
        sd = d / f"s{s}"
        sd.mkdir(exist_ok=True)
        for ext in ("mmn", "eig"):
            ln = sd / f"wannier90.{s}.{ext}"
            ln.unlink(missing_ok=True)
            ln.symlink_to(Path("..") / f"wannier90.{s}.{ext}")
        amn = sd / f"wannier90.{s}.amn"
        amn.unlink(missing_ok=True)
        trimmed = rec.get("trimmed_amn")
        if trimmed and not str(trimmed).startswith("SKIPPED"):
            amn.symlink_to(Path("..") / Path(trimmed).name)
        else:
            amn.symlink_to(Path("..") / f"wannier90.{s}.amn")
        shutil.copy2(rec["win"], sd / f"wannier90.{s}.win")
        spins.append(s)
    return spins


def w90pump(dest: Path, np_: int, threads: int, w90_bin: str, qos: str,
            partition: str, walltime: str, max_queue: int) -> None:
    """Pump offline cases: EXPORTED -> solve submitted -> wouts back."""
    state = load_state(dest)
    in_queue = squeue_ids()
    if in_queue is None:
        print("queue query failed - pumping nothing (occupancy unknown)")
        return
    active = sum(1 for v in state.values() if v.get("w90_jobid") in in_queue)
    prep = None
    for name, rec in sorted(state.items()):
        if rec.get("w90") != "offline":
            continue
        d = dest / name
        # export job finished since the last look?
        if rec["status"] == "SUBMITTED" and rec.get("jobid") not in in_queue:
            rec["status"] = offline_export_triage(d)
            print(f"[{name}] export job done -> {rec['status']}")
        # solve finished since the last look?
        if rec["status"] == "W90_SUBMITTED" and \
                rec.get("w90_jobid") not in in_queue:
            spins = rec.get("w90_spins") or [1, 2]
            wouts = {s: d / f"s{s}" / f"wannier90.{s}.wout" for s in spins} \
                if spins != [0] else {0: d / "wannier90.wout"}
            if all(_wout_done(p) for p in wouts.values()):
                for s, p in wouts.items():
                    if s:                       # per-spin: copy back for pyleos
                        shutil.copy2(p, d / p.name)
                rec["status"] = "OK"
                print(f"[{name}] w90 solve done - wout channel complete")
            else:
                rec["status"] = "FAIL_W90"
                missing = [str(p) for p in wouts.values() if not _wout_done(p)]
                rec["status_note"] = "no Final State: " + ", ".join(missing)
                print(f"[{name}] FAIL_W90 ({rec['status_note']})")
        # ready for a solve?
        if rec["status"] != "EXPORTED":
            continue
        if active >= max_queue:
            print(f"w90 queue cap {max_queue} reached; run w90 again later")
            break
        if prep is None:
            prep = _load_prep()
        try:
            summary = prep.generate(d, out_dir=d, overwrite=True,
                                    do_trim_amn=True, crosscheck=True)
        except prep.PrepError as exc:
            rec["status"] = "FAIL_PREP"
            rec["status_note"] = str(exc)
            print(f"[{name}] FAIL_PREP: {exc}")
            continue
        if summary["ispin"] == 2:
            spins = _stage_spin_dirs(d, summary)
            run_block = "\n".join(
                W90_RUN_SPIN.format(s=s, off=(0 if s == 1 else "$SPAN"),
                                    w90=w90_bin) for s in spins)
            cpus = 2 * np_ * threads
        else:                       # ISPIN=1: solve in place, shared seedname
            spins = [0]
            run_block = W90_RUN_SINGLE.format(w90=w90_bin)
            cpus = np_ * threads
        (d / "w90.sbatch").write_text(W90_SBATCH_TEMPLATE.format(
            name=name[:12], partition=partition, qos=qos, cpus=cpus,
            walltime=walltime, np=np_, threads=threads, run_block=run_block))
        r = subprocess.run(["sbatch", "w90.sbatch"], cwd=d,
                           capture_output=True, text=True)
        if r.returncode == 0:
            rec.update(status="W90_SUBMITTED",
                       w90_jobid=r.stdout.strip().split()[-1],
                       w90_spins=spins)
            active += 1
            print(f"[{name}] w90 solve submitted {rec['w90_jobid']} "
                  f"(np={np_} x {threads} threads x {len(spins)} spin(s))")
        else:
            print(f"[{name}] w90 sbatch failed: {r.stderr.strip()[:80]}")
    save_state(dest, state)


RETRY_LADDER = {
    "FAIL_DFT": [{"ALGO": "Normal"}, {"ALGO": "All", "EDIFF": 1e-5}],
    # NBANDS_SCALE is not an INCAR tag: retry() resolves it to a concrete
    # NBANDS (docstring promise: SCDM off + larger NBANDS margin)
    "FAIL_W90": [{"EOS_SCDM_ONLY": ".FALSE.", "NBANDS_SCALE": 1.3}],
    "FAIL_LEOS": [{"EOS_NELECT_ABORT": 0.5}],  # diagnose-mode rerun
}


def resolve_nbands(d: Path, tags: dict) -> int | None:
    """NBANDS actually in effect for a case: the INCAR tag if set, else the
    value VASP printed in OUTCAR, else one found in the EOSCAR header."""
    try:
        return int(float(tags["NBANDS"]))
    except (KeyError, ValueError):
        pass
    outcar = d / "OUTCAR"
    if outcar.is_file():
        with outcar.open(errors="replace") as fh:
            for line in fh:
                m = re.search(r"NBANDS\s*=\s*(\d+)", line)
                if m:
                    return int(m.group(1))
    eoscar = d / "EOSCAR"
    if eoscar.is_file():
        with eoscar.open(errors="replace") as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                m = re.search(r"NBANDS\s*=\s*(\d+)", line)
                if m:
                    return int(m.group(1))
    return None


def retry(dest: Path) -> None:
    state = load_state(dest)
    for name, rec in state.items():
        st = rec["status"]
        if st not in RETRY_LADDER:
            continue
        ladder = RETRY_LADDER[st]
        if rec["retries"] >= len(ladder):
            print(f"[{name}] ladder exhausted ({st})")
            continue
        overrides = dict(ladder[rec["retries"]])   # copy: rungs are shared
        d = dest / name
        incar = d / "INCAR"
        tags = {}
        for line in incar.read_text().splitlines():
            if "=" in line:
                k, v = line.split("=", 1)
                tags[k.strip()] = v.strip()
        scale = overrides.pop("NBANDS_SCALE", None)
        if scale is not None:
            nb = resolve_nbands(d, tags)
            if nb is None:
                # nothing to scale: fall back to a tag that still guarantees
                # the rerun differs from the failed one, and say so
                overrides.setdefault("EOS_SCDM_ONLY", ".FALSE.")
                print(f"[{name}] NBANDS not found in INCAR/OUTCAR/EOSCAR - "
                      f"skipping the x{scale} margin, retrying with "
                      f"EOS_SCDM_ONLY=.FALSE. only")
            else:
                # strictly larger margin even when rounding would not move
                overrides["NBANDS"] = max(int(round(nb * scale)), nb + 1)
        tags.update({k: str(v) for k, v in overrides.items()})
        with incar.open("w") as fh:
            for k, v in tags.items():
                fh.write(f"{k} = {v}\n")
        rec.update(status="RETRY", retries=rec["retries"] + 1)
        print(f"[{name}] staged retry #{rec['retries']} ({st} -> {overrides})")
    save_state(dest, state)


def harvest(dest: Path, out: Path) -> None:
    state = load_state(dest)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, rec in sorted(state.items()):
        d = dest / name
        st = triage(d) if (d / "run.sbatch").is_file() else rec["status"]
        keep = out / name
        keep.mkdir(exist_ok=True)
        for f in ("EOSCAR", "POSCAR", "INCAR", "vasp.out", "print-out"):
            if (d / f).is_file():
                shutil.copy2(d / f, keep / f)
        for w in d.glob("wannier90*.wout"):
            shutil.copy2(w, keep / w.name)
        rows.append((name, st))
    with (out / "summary.tsv").open("w") as fh:
        fh.write("case\tstatus\n")
        for name, st in rows:
            fh.write(f"{name}\t{st}\n")
    print(f"harvested {len(rows)} cases -> {out}/summary.tsv")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("stage")
    p.add_argument("--dest", type=Path, required=True)
    p.add_argument("--cases", type=Path, required=True)
    p.add_argument("--qos", default="rush-gpu")
    p.add_argument("--walltime", default="02:00:00")
    p.add_argument("--mode", default="full", choices=["full", "scdm", "tiered"])
    p.add_argument("--offline-w90", action="store_true",
                   help="GPU job exports only (LWRITE_MMN_AMN=T, "
                        "LWANNIER90_RUN=F); solve later with the w90 command")
    p = sub.add_parser("escalate")
    p.add_argument("--dest", type=Path, required=True)
    p = sub.add_parser("w90", help="pump offline cases: generate per-spin "
                                   "wins, submit CPU solves, pull wouts back")
    p.add_argument("--dest", type=Path, required=True)
    p.add_argument("--np", type=int, default=8,
                   help="MPI ranks per spin (k-point parallelism caps at nk)")
    p.add_argument("--threads", type=int, default=1,
                   help="MKL/OMP threads per rank (needs a threaded-MKL "
                        "wannier90.x to matter)")
    p.add_argument("--w90-bin", default=W90_BIN)
    p.add_argument("--qos", default="rush-cpu")
    p.add_argument("--partition", default="DSPRHBM")
    p.add_argument("--walltime", default="02:00:00")
    p.add_argument("--max-queue", type=int, default=8)
    p = sub.add_parser("submit")
    p.add_argument("--dest", type=Path, required=True)
    p.add_argument("--max-queue", type=int, default=8)
    p = sub.add_parser("status")
    p.add_argument("--dest", type=Path, required=True)
    p = sub.add_parser("retry")
    p.add_argument("--dest", type=Path, required=True)
    p = sub.add_parser("harvest")
    p.add_argument("--dest", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    if a.cmd == "stage":
        stage(a.dest.expanduser(), a.cases, a.qos, a.walltime, a.mode,
              offline_w90=a.offline_w90)
    elif a.cmd == "escalate":
        escalate(a.dest.expanduser())
    elif a.cmd == "w90":
        w90pump(a.dest.expanduser(), a.np, a.threads, a.w90_bin, a.qos,
                a.partition, a.walltime, a.max_queue)
    elif a.cmd == "submit":
        submit(a.dest.expanduser(), a.max_queue)
    elif a.cmd == "status":
        status(a.dest.expanduser())
    elif a.cmd == "retry":
        retry(a.dest.expanduser())
    else:
        harvest(a.dest.expanduser(), a.out.expanduser())


if __name__ == "__main__":
    main()
