#!/usr/bin/env python3
"""Turn a MACE production trajectory into LEOS DFT single-point jobs.

Reads traj.extxyz, picks every Nth frame, writes one VASP directory per
frame (POSCAR + LEOS INCAR + KPOINTS + POTCAR symlink), plus a chain
sbatch that runs frames sequentially inside one allocation (cheap frames,
one queue slot). Per-system Hubbard U on the metal site follows the MP
conventions; oxygen stays U-free here (the polaron study handles that
question separately).

Run on sai-pod:
  python3 frames_to_dft.py m5b_md/li2mno3/x0.25/prod/traj.extxyz \
      --system li2mno3 --every 5 --out m5b_dft/li2mno3/x0.25 \
      --potcar ~/potcars/POTCAR_li2mno3 [--scdm-only]
"""
from __future__ import annotations

import argparse
from pathlib import Path

from ase.io import read, write

U_TABLE = {  # element -> (LDAUL, LDAUU) applied on the TM only
    "li2mno3": {"Mn": 3.9},
    "li2ruo3": {"Ru": 3.0},
    "nmc": {"Ni": 6.2, "Mn": 3.9, "Co": 3.32},
}
MAG = {"Mn": 4.0, "Ru": 1.0, "Ni": 1.0, "Co": 1.0, "Li": 0.0, "O": 0.3}

SBATCH = """#!/bin/bash
#SBATCH --job-name=fr-{tag}
#SBATCH --partition=16V100
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --gpus-per-node=4
#SBATCH --qos={qos}
#SBATCH --time={hours}:00:00
#SBATCH --output=chain-out
#SBATCH --error=chain-out
set -eo pipefail
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
VASP=/home/sjtu-gongxueqing/zhangxianyang/apps/vasp.6.4.3_leos_v2/bin/vasp_std
RUN="mpirun -np ${{SLURM_NTASKS}} -map-by ${{MAP_OPT:-ppr:1:gpu}} -mca coll_hcoll_enable 0"
cd {root}
for f in frame_*; do
  if [ -f "$f/EOSCAR" ]; then echo "skip $f (done)"; continue; fi
  echo "=== $f ==="
  (cd "$f" && $RUN $VASP) || echo "FRAME_FAILED $f"
done
echo CHAIN_DONE_{tag}
"""


def incar_text(system: str, symbols, scdm_only: bool) -> str:
    species = list(dict.fromkeys(symbols))          # POSCAR species order
    counts = [symbols.count(s) for s in species]
    mag = " ".join(f"{c}*{MAG.get(s, 0.0)}" for s, c in zip(species, counts))
    u = U_TABLE[system]
    ldaul = " ".join("2" if s in u else "-1" for s in species)
    ldauu = " ".join(f"{u.get(s, 0.0)}" for s in species)
    lines = [
        "SYSTEM = M5-B frame single point",
        "PREC = Accurate", "ENCUT = 500", "EDIFF = 1e-6",
        "ISMEAR = 0", "SIGMA = 0.05", "ISYM = 0",
        "NSW = 0", "IBRION = -1",
        "ISPIN = 2", f"MAGMOM = {mag}",
        "LDAU = .TRUE.", "LDAUTYPE = 2",
        f"LDAUL = {ldaul}", f"LDAUU = {ldauu}",
        f"LDAUJ = {' '.join('0' for _ in species)}", "LMAXMIX = 4",
        "LWAVE = .FALSE.", "LCHARG = .FALSE.",
        "LEOS = .TRUE.", "LSCDM = .TRUE.",
        "LWANNIER90_AUTO_WINDOW = .TRUE.",
        "EOS_WF_DETAIL = .TRUE.", "EOS_PROJ_WEIGHT = .TRUE.",
    ]
    if scdm_only:
        lines.append("EOS_SCDM_ONLY = .TRUE.")
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("traj", type=Path)
    ap.add_argument("--system", required=True, choices=sorted(U_TABLE))
    ap.add_argument("--every", type=int, default=5)
    ap.add_argument("--max-frames", type=int, default=30)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--potcar", type=Path, required=True)
    ap.add_argument("--scdm-only", action="store_true")
    ap.add_argument("--qos", default="rush-gpu")
    ap.add_argument("--hours", type=int, default=20)
    args = ap.parse_args()

    frames = read(args.traj, index="::")
    picked = frames[:: args.every][: args.max_frames]
    args.out.mkdir(parents=True, exist_ok=True)
    for k, atoms in enumerate(picked):
        d = args.out / f"frame_{k:03d}"
        d.mkdir(exist_ok=True)
        write(d / "POSCAR", atoms, format="vasp", direct=True, sort=True)
        # sort=True writes species alphabetically; INCAR arrays must match
        symbols = sorted(atoms.get_chemical_symbols())
        (d / "INCAR").write_text(incar_text(args.system, symbols,
                                            args.scdm_only))
        (d / "KPOINTS").write_text("auto\n0\nG\n2 2 2\n0 0 0\n")
        pot = d / "POTCAR"
        if not pot.exists():
            pot.symlink_to(args.potcar.expanduser().resolve())
    tag = f"{args.system}-{args.out.name}"
    (args.out / "chain.sbatch").write_text(
        SBATCH.format(tag=tag[:10], root=args.out.resolve(), qos=args.qos,
                      hours=args.hours))
    print(f"staged {len(picked)} frames -> {args.out} "
          f"(chain.sbatch ready, scdm_only={args.scdm_only})")


if __name__ == "__main__":
    main()
