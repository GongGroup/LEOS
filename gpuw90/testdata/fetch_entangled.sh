#!/bin/sh
# Fetch the entangled gold-set raw inputs (.mmn/.amn) from sai-pod.
#
# The small companions (.win/.eig and the wannier90.x reference *_ref.wout) are
# tracked in git; the 10-300 MB .mmn/.amn are deliberately re-ignored in
# .gitignore and must be pulled from the cluster before running
# check_dis_against_w90.py:
#
#   cd gpuw90/testdata && sh fetch_entangled.sh
#
# Source of truth on sai-pod (see docs/reports/20260710_gpu_feasibility.md):
#   li2mno3 : ~/LEOS/w90_offline/bench/li2mno3_x000_f000/   (VASP LWRITE_MMN_AMN export)
#   vircef  : ~/LEOS/w90_offline/offline_runs/vircef/
# Reference .wout = serial run_np1_s{1,2}/wannier90.{1,2}.wout next to the same
# inputs (symlinked), i.e. bit-identical provenance.

set -e
BASE=/home/sjtu-gongxueqing/zhangxianyang/LEOS/w90_offline

mkdir -p li2mno3 vircef

ssh sai-pod "cd $BASE/bench/li2mno3_x000_f000 && tar czf - \
    wannier90.1.mmn wannier90.1.amn wannier90.2.mmn wannier90.2.amn" \
    | tar xzf - -C li2mno3

ssh sai-pod "cd $BASE/offline_runs/vircef && tar czf - \
    wannier90.1.mmn wannier90.1.amn wannier90.2.mmn wannier90.2.amn" \
    | tar xzf - -C vircef

echo "done:"
ls -l li2mno3/*.mmn li2mno3/*.amn vircef/*.mmn vircef/*.amn
