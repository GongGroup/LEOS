FeO2 单结构静态计算与 LEOS v3，2026-10-15

用户授权：检查已复制的 FeO2 文件，保持关键参数，以 SAI 的 VASP 6.4.3 完成一次静态单点和 LEOS，提交作业并取得数据。
源目录：/home/zhangyanglu/zack/leos.sample/FeO2-sai
计算目录：/home/zhangyanglu/zack/leos.sample/FeO2-sai/static_leos_v3_20261015_01
程序：/home/zhangyanglu/zack/app/vasp.6.4.3-leos-v3/bin/vasp_std
程序 SHA256：26ebfa0d8c79a170a3b9f04af6b2ef1c9af40293b1a705bf6b51d8438208d9b6

POSCAR、POTCAR、KPOINTS 原样复制。原 INCAR 的唯一替换是输出开关 LWAVE=False 改为 True。
新增 NSW=0、IBRION=-1 指明静态计算；LCHARG、LAECHG 保存电荷数据；LEOS、EOS_WF_DETAIL、EOS_PROJ_WEIGHT 开启原生 LEOS 和详细输出。
保留 ENCUT=520、EDIFF=0.00015、ISPIN=2、MAGMOM=1*5.0 2*0.6、PBE+U(Fe:5.3 eV)、ISMEAR=-5、SIGMA=0.05、NELM=100、LREAL=Auto，以及其余原设置。ISYM 未覆盖，KPOINTS 保留 Gamma 8x8x5。
不读取旧 CHGCAR 或空 WAVECAR。新目录独立进行自洽静态计算。
使用已安装 LEOS v3 的自动 Wannier 默认协议，不传入真实中心、不预先指定 num_wann，不修改 EOS 阈值或参考氧化态。
课题组账户 sjtu-gongxueqing；improper-gpu；4V100 或 8V100V0；1 GPU、1 MPI rank。

结果包括 OUTCAR、CHGCAR、AECCAR*、WAVECAR、EOSCAR、Wannier 文件，以及 result_summary.json、centres.tsv、oxidation_states.tsv。
本次不修改机器学习程序、不训练、不自动纳入训练集。
