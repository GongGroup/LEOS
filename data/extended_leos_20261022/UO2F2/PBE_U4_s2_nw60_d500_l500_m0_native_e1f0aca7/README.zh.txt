UO2F2 原生 VASP + LEOS v3 静态计算
源目录：/home/zhangyanglu/zack/leos.sample/UO2F2/relax_ao2021_20261015
计算目录：/home/zhangyanglu/zack/leos.sample/UO2F2/relax_ao2021_20261015/leos_v3_500_500_20261016_01
采用已收敛弛豫的 CONTCAR 作为 POSCAR，保持 KPOINTS 和 POTCAR 字节一致。
原输入仅将 IBRION=2 改为 -1、NSW=200 改为 0，并增加 LEOS 和 Wannier 控制。
解纠缠 dis_num_iter=500；局域化 num_iter=500。解纠缠可提前收敛。
Ueff=4 eV；ISPIN=2；初始 MAGMOM=6*0 6*0 3*2；ENCUT=600；k 网格 7×7×2。
保留 ISTART=0 / ICHARG=2，从头进行静态自洽，不复用弛豫 WAVECAR。
LEOS 使用程序原生默认 beta=10；不套用 FeO2 调参结果。
程序：/home/zhangyanglu/zack/app/vasp.6.4.3-leos-v3/bin/vasp_std
资源：sjtu-gongxueqing / improper-gpu，4V100 或 8V100V0，1 GPU、1 MPI，时限 4 小时。
结束后自动汇总 result_summary.json、centres.tsv、oxidation_states.tsv。
