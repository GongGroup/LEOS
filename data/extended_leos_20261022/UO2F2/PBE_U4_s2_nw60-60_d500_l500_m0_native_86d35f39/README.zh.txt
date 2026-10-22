用户明确修正：放宽 kmesh_tol 至 1e-5（原默认 1e-6）。
前次误将方向理解为缩小至 1e-8，作业 1805461 已请求取消，保留日志。
本次在 WANNIER90_WIN 中增加 kmesh_tol = 1e-5。
解纠缠和局域化均保持 500；结构、k 网格、POTCAR、Ueff=4 eV、自旋设置、beta=10 均保持。
新目录从头静态自洽，不复用未完成的 WAVECAR / CHGCAR。
仍使用 VASP6.4.3 LEOS v3，课题组账户 sjtu-gongxueqing、improper-gpu、单张 V100。
