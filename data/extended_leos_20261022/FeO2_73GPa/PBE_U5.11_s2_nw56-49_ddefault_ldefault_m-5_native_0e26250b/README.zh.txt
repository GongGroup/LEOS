FeO2 73 GPa：原生 VASP 6.4.3 LEOS v3 单次计算

用户授权：对 pbe_u5p11_static_20261015_02 做一次直接 LEOS 计算。
原目录：SAI:/home/zhangyanglu/zack/leos.sample/FeO2-73Gpa/pbe_u5p11_static_20261015_02
新目录：原目录/leos_v3_20261015_01

保留原 POSCAR、POTCAR、KPOINTS 及全部既有 INCAR 字段。
原结构为 a=4.3488 Å 的立方晶胞，8 O + 4 Fe；PBE+U，Fe Ueff=5.11 eV；
ENCUT=600 eV，EDIFF=1e-8，Gamma 5×5×5，ISPIN=2，原 MAGMOM=8*0.0 4*1.0。
只增加 LEOS、EOS_WF_DETAIL、EOS_PROJ_WEIGHT 三个开关。
保持原 ISTART=0、ICHARG=2，重新自洽后由集成 Wannier90 和 LEOS 计算；
未将历史中心、整数氧化态或 num_wann 作为手动输入。
不覆盖历史结果，不优化几何，不进行机器学习。

作业用课题组账户 sjtu-gongxueqing、improper-gpu、单 V100。
原静态 VASP 耗时 114.353 秒；本次加入 Wannier 的耗时待实测。
结果是否科学可用，须结合电子数误差、局域化展宽和原生警告判断。
