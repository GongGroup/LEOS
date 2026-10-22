CuS PBE+U5 离线2000步单轮

用户授权：用离线LEOS增加到2000步，观察四配位Cu是否出现+2。
作业：1888756；SAI / zhangyanglu / sjtu-gongxueqing；单V100，improper-gpu，2小时上限。
远端：/org/sjtu-gongxueqing/zhangyanglu/machinelearning/CuS_PBEU5_gpuw90_2000_20261021_01
源：/org/sjtu-gongxueqing/zhangyanglu/machinelearning/CuS_PBEU5_leos_500_20261021_01/iter_500

使用原.mmn、.amn、.eig、.win，不重算VASP/DFT。
Ueff=5 eV、ISPIN=1、52个Wannier函数、80条带、576个k点和原能窗不变。
保留原中性DFT的102电子作为flow目标，每个非自旋极化Wannier函数提供最多2电子。
不将旧WF的分数占据或投影权重套到新中心；本次未重算分数占据。
解纠缠上限2000步，可按原阈值提前收敛；局域化固定2000步，每500步保存位置和分配。
采用现有GPUw90固定步长最速下降算法，step=1，complex128。
从原AMN初始投影重做离线优化，并非读取VASP第500步CHK再续算1500步。
离线优化与原生Wannier优化方法不同，因此同时保存离线500/1000/1500/2000步供同算法比较。
flow固定beta=10，不强制Cu价态，不修改原LEOS分配算法。
原生与离线参考中的两个S位点差异保留为记录，不阻止新计算。

输出：provenance.json、disentangle_summary.json、checkpoint_500/1000/1500/2000.json、centres.tsv、summary.json。
收敛报告标准：末10步总展宽变化绝对值均<=1e-6 Å²，不宣称全局最优。
估时：按此前125k点GPUw90每步约0.11秒及本次576k点规模粗估10—30分钟，需实际进度修正。
