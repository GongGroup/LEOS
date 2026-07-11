#!/bin/bash
#SBATCH --job-name=gpuw90-bench
#SBATCH --partition=16V100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=4
#SBATCH --qos=flood-gpu
#SBATCH --output=bench-out
#SBATCH --error=bench-out
# gpuw90 Phase-2 CUDA benchmark on sai-pod V100 (issue #14, non-blocking item).
#
# Partition rules (see hpc-sai-1344 notes):
#  * 16V100 requires >=4 GPUs/job -- we request 4 and use 1 (CUDA_VISIBLE_DEVICES).
#    flood-gpu = 4 h walltime, plenty (bench needs minutes).
#  * NO --cpus-per-task/--mem* flags: the site hook auto-cancels GPU jobs with them.
#  * cd into the case dir before sbatch, or the relative --output lands elsewhere.
#
# One-time env (login node):
#   export MAMBA_ROOT_PREFIX=$HOME/bin/miniforge3
#   source $MAMBA_ROOT_PREFIX/etc/profile.d/conda.sh && source $MAMBA_ROOT_PREFIX/etc/profile.d/mamba.sh
#   mamba create -n gpuw90 -y python=3.11
#   mamba activate gpuw90
#   pip install numpy torch --index-url https://download.pytorch.org/whl/cu124
#   python -c "import torch; assert 'sm_70' in torch.cuda.get_arch_list()"
# TRAP: V100 is sm_70. The site mace0.3.15 env ships torch 2.9.1+cu128, whose
# wheels DROP sm_70 ("no kernel image" at runtime). Stay on cu124-era wheels
# (torch<=2.6 series); the assert above is the guard.
#
# Data: the .mmn/.amn/.eig/.win already live on this cluster --
#   li2mno3: ~/LEOS/w90_offline/bench/li2mno3_x000_f000   (seeds wannier90.{1,2})
#   vircef : ~/LEOS/w90_offline/offline_runs/vircef
# Code: rsync the local gpuw90/ tree to ~/LEOS/gpuw90 (adjust GPUW90 below).

set -euo pipefail

GPUW90=$HOME/LEOS/gpuw90
DATA=$HOME/LEOS/w90_offline/bench/li2mno3_x000_f000

export MAMBA_ROOT_PREFIX=$HOME/bin/miniforge3
set +u
source $MAMBA_ROOT_PREFIX/etc/profile.d/conda.sh
source $MAMBA_ROOT_PREFIX/etc/profile.d/mamba.sh
mamba activate gpuw90
set -u

export CUDA_VISIBLE_DEVICES=0
cd "$GPUW90"

echo "== GPU sanity =="
python -c "import torch; print(torch.__version__, torch.version.cuda, torch.cuda.get_device_name(0))"

echo "== correctness gates on V100 (li2mno3 spin 1) =="
# the serial reference .wout lives in run_np1_s1/ next to the (symlinked) inputs;
# the check script expects <seed>_ref.wout in the data dir
ln -sf run_np1_s1/wannier90.1.wout "$DATA/wannier90.1_ref.wout"
ln -sf run_np1_s2/wannier90.2.wout "$DATA/wannier90.2_ref.wout"
python check_dis_against_w90.py "$DATA" --seed wannier90.1 --device cuda \
  || { echo 'GOLD GATE FAILED ON CUDA -- stop, compare against CPU run'; exit 1; }

echo "== timing: li2mno3 (production scale, 288 WF) =="
python bench_gpu.py "$DATA" --seed wannier90.1 --device cuda --dis-iters 50 --wann-iters 200
python bench_gpu.py "$DATA" --seed wannier90.1 --device cpu  --dis-iters 10 --wann-iters 20

echo "done $(date)"
