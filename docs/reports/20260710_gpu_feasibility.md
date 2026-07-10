> **归档注记** — 生成日期：2026-07-10；生成者：主会话派出的 GPU 可行性评估 sub-agent；原路径：`/private/tmp/claude-501/-Users-mastreina-Desktop-EOS/04bfbd96-cc02-4d31-8159-643385b8af3d/scratchpad/gpu_feasibility.md`（原件保留未动）。由成果收集会话拷贝入库，总索引见 [20260710_index.md](20260710_index.md)。

# LEOS Wannier 化管线 GPU 化可行性评估(决策级报告)

- 评估对象:把 Wannier 化管线三段(SCDM 初猜、disentanglement、wannierise 幺正优化)移到 GPU
- 环境:VASP 6.4.3 OpenACC(nvhpc 24.3)+ wannier90 3.1 库模式(串行 CPU);sai-pod 16×V100-SXM2-16GB/节点
- 方法:只读源码分析(本地 mlwf.F、sai-pod 上 wannier90-3.1.0 源码与 ws_bench 实测输出)+ 网络文献/库文档调研。未改任何文件、未提交任何作业。
- 标注约定:**[实证]** = 源码行号 / 实测输出 / 官方文档直接支持;**[推断]** = 由实证数据外推或经验估计,使用前需按文中标注核实。

---

## ① 一页结论

### 三条子路线四格评分

| 子路线 | 可行性 | 收益 | 工程量 | 风险 | 一句话判断 |
|---|---|---|---|---|---|
| 1. SCDM(QRCP)GPU 化 | 中(cuSOLVER 无 geqp3,只有 MAGMA 混合实现或自写随机化 QRCP) | **极低**(上限 ≈8.3 s/帧,占 w90 段 1.5%) | 中~高 | **高**(枢轴翻转改变初猜落入的局部极小 + 5~10 GB 级矩阵挤爆 V100 HBM) | **不做**:收益天花板 8 秒,却背上全管线最重的数值与内存风险 |
| 2. disentanglement GPU 化 | 中高(batched zgemm 结构清晰;~330 维 batched 厄米特征分解在 V100/CUDA 12.3 工具链上没有好现成库) | 中(实测 34~140 s/自旋/帧) | 高(fork wannier90 或重写) | 中(冻结窗口边缘近简并子空间选择翻转,正撞 LEOS 近平局告警区) | 能做,但单独做不划算;只应作为路线 3 的搭车项 |
| 3. wannierise GPU 化 | 高(每迭代 = 全 k 点 batched zgemm + 每 k 一个 NW×NW zheev,教科书式可批量) | 中高(wannierise 迭代是最大单项,~150~170 s/自旋/帧,×2 自旋) | 高(同上;若做成独立后处理工具则中) | 中(非凸多极小 + 浮点顺序差异,可用金标集回归兜住) | 三条里唯一值得考虑的;**且只应以"独立 GPU 后处理工具"形态做,不碰 wannier90 上游、不在 VASP 内重写** |

### 总判断

**"GPU 化 Wannier 管线"这个方向本身可行,但作为当前优先事项不可靠——因为同量级收益存在多条零代码/低代码的 CPU 路线,而这些路线还没有用掉。**

关键事实链:
1. 554 s/帧的 w90 段目前跑在**单核串行库模式**上([实证] wannier_lib.F90:61/229 "library mode currently works ONLY in serial";mlwf.F:1874 只在 masternode 调用),期间整节点 16 块 V100 全部闲置。
2. 而 wannier90 3.1 **standalone 可执行文件的 dis/wann 主循环上游本来就是 MPI k 点并行的**([实证] disentangle.F90/wannierise.F90 中 comms_array_split + nkp_loc 局部循环),VASP 侧也已有 LWRITE_MMN_AMN 文件出口([实证] mlwf.F:378)。即:把 w90 段搬出 VASP 离线 MPI 跑,是**零新代码**的 6~14× (nk=8 体系)到 ~20×(nk=32 体系)。
3. 再叠加帧间 warm start(CP2K `&LOCALIZE USE_HISTORY` 的同类做法,AIMD 领域标准实践)可再省大部分迭代数;SCDM-only 分层(已有,8.3 s/帧 = 67×)覆盖闸门放行的帧。
4. GPU 化在 V100 上对本工作负载的现实加速比是 **5~15×(推断,见 §3/§4)**,不是 50×——它拿到的东西和第 2、3 条 CPU 路线重叠,却要付出 fork 上游或重写的工程与数值风险。

**执行顺序建议**:先做零代码离线 MPI + 自旋并发(第 2 条),再做 warm start(小工具级工程),SCDM-only 分层继续当快速路径;只有当这些做完后 w90 段仍是轨迹吞吐瓶颈、且目标体系 nk 较大时,才立项"路线 3+2 合并的独立 GPU 后处理工具";路线 1(SCDM QRCP GPU)明确不做。

---

## ② 规模画像与基线拆解(实证锚点)

### 2.1 两个体系的矩阵尺寸

| 量 | Ce18O36Co4CH4(63 原子,"大"标定) | Li2MnO3(96 原子,554 s 实测来源) |
|---|---|---|
| NBANDS (NB) | 332 | 352 [实证 .wout] |
| num_wann (NW)/自旋 | ~250 | 288 [实证 .wout] |
| k 网格(w90 侧全网格 nk) | 4×4×2 = 32 | 2×2×2 = 8 [实证 .wout] |
| ISPIN | 2 | 2(wannier90.1/.2 两套文件)[实证] |
| 实空间精细网格 Nr | ~1e6 量级(题设) | 未查 |
| 单个 NW×NW 复矩阵 | 1.0 MB | 1.3 MB |
| M_mn 全集 NB²·nntot·nk·16B(nntot 取 8~12) | 0.45~0.68 GB | 0.13~0.19 GB |
| SCDM 矩阵 NB×Nr(复16B) | ~5.3 GB/自旋 + PSI 同量级镜像 | 未查(同构更小) |

**要点**:单个矩阵(250~350 维)对 V100 太小,喂不饱 GPU;唯一的并行来源是 nk×(每迭代 gemm 数) 的 batch 维度。nk=8 的 Li2MnO3 batch 很薄,nk=32 的 Ce 体系稍好。dis/wann 的全部工作集 ≤2 GB,HBM 不是问题;**SCDM 的 NB×Nr 矩阵(≈5.3 GB/自旋,外加同量级的 PSI 副本)才是 HBM 硬约束的撞点**。

### 2.2 554 s/帧的实测拆解([实证] sai-pod ~/LEOS/ws_bench/on/wannier90.{1,2}.wout)

| 项 | 自旋 1 | 自旋 2 |
|---|---|---|
| Time to disentangle | 33.97 s | 139.62 s |
| Total Execution Time(该 .wout 含两次运行) | 184.65 s / 246.93 s | 312.10 s / 251.27 s |
| `<-- DIS` 行数(含表头) | 59 | 203 |
| `<-- CONV` 行数(含表头) | 936 | 936 |

- 每帧两自旋合计:184.65+312.10 ≈ 497 s 或 246.93+251.27 ≈ 498 s,与题设 554 s/帧吻合(差额 ≈ SCDM 8.3 s + M/A 矩阵构建与写盘)[推断:差额归属未逐项计时]。
- 扣除表头后每次运行约 **26~98 次 dis 迭代、约 460 次 wannierise 迭代/自旋** [推断:按每文件两次运行均分]。
- 由此:dis ≈ 1.1~1.4 s/迭代,wannierise ≈ 0.32 s/迭代。按每 wannierise 迭代 ≈ nk×(2·nntot+~4) 个 NW³ zgemm ≈ 30 Gflop 估算,0.32 s 对应 ~95 GF/s——**超过单核 FP64 峰值,提示库调用可能已隐性吃到线程化 BLAS**[推断:需核实 w90 库链接的 BLAS 与 OMP_NUM_THREADS,若属实,GPU 化的对比基线还要再抬高]。

### 2.3 算法结构确认([实证] 源码)

- **VASP 侧 SCDM**(mlwf.F `CALC_WANNIER_TRAFO_SCDM`,3794–4003 行):仅 node 1 干活(3842 行);在锚点(Γ)对 **NB×Nr 全实空间网格矩阵做一次 xGEQP3**(3887/3917 行,经 mathtools.F 的 RRQR_DECOMPOSER);随后每 k 只取枢轴列组装 NB×NW 小矩阵,xGESDD SVD 正交化,T=U·V†(3950–4000 行)。即 **QRCP 每帧每自旋只做一次,是 SCDM-only 8.3 s 的主体**;每 k 的 SVD 是 332×250 级,开销可忽略。
- **wannier90 disentangle**(disentangle.F90,3139 行):`dis_extract` 大迭代循环(~1700 行起);每迭代对每个(本地)k 组装 Z 矩阵(num_wann×ndimwin×ndimwin 级 zgemm,147–218/1965–1968 行等)+ 一次 ZHPEVX 打包厄米特征分解(维度 ≤ ndimwin ≈ NB);k 点经 comms_array_split 分布,czmat_in_loc 按 nkp_loc 存放。
- **wannier90 wannierise**(wannierise.F90,3507 行):`do iter=1,num_iter`(416 行);每迭代每 k:cdq 反厄米梯度 → **zheev(NW×NW)**(1228 行)重构 exp(dW) → 3 个 NW³ zgemm(1246/1253/1289 行)更新 U → 每 (k,b) 邻居 2 个 NW³ zgemm 旋转 M(1299–1301 行);U 以 u_matrix_loc 按 k 分布,迭代内 comms_bcast(645 行)。
- 结论:**dis/wann 每迭代 = "batch(nk) × 小 zgemm(NW³/NB·NW²)" + "batch(nk) × 一次 ~250–350 维厄米特征分解"**,与题设预判一致;标准 MPI k 并行已在 standalone 路径存在,只是库模式被锁死为串行。

---

## ③ 子路线 1:SCDM 的 GPU 化 —— 结论:不做

**算法与规模**:一次性 QRCP 于 NB×Nr(Ce:332×~1e6,复 16B ≈ 5.3 GB;另有 Nr×NB 的 PSI 副本 ≈ 5.3 GB)。QRCP 之后的每 k SVD 是毫秒级,不值得动。

**库与前例(调研结果)**:
- **cuSOLVER 没有 geqp3(带列主元 QR)**,只有不带主元的 geqrf;这是长期公开的缺口([实证] [cuSOLVER 文档](https://docs.nvidia.com/cuda/cusolver/index.html)无该 API;[PyTorch issue #82092](https://github.com/pytorch/pytorch/issues/82092) 与 [JAX PR #25955](https://github.com/jax-ml/jax/pull/25955) 均因此转向 MAGMA)。
- **MAGMA 有 (z)geqp3 / geqp3_gpu**([文档](https://icl.utk.edu/projectsfiles/magma/doxygen/group__magma__geqp3.html)):混合实现,panel 在 CPU、trailing update 在 GPU;因主元判据要求 trailing update 完成后才能开始下一 panel,**无法通信-计算重叠,基本受内存带宽限制**([UC Davis 分析](https://www.cs.ucdavis.edu/~bai/publications/tomasbaihernandez13.pdf))。预期相对单核 CPU 仅低个位数倍加速[推断]。
- **随机化 QRCP(RQRCP)**:Duersch–Gu([arXiv:1509.06820](https://arxiv.org/abs/1509.06820))用随机投影 B=ΩA 做主元选择,达到接近**无主元 QR 的速度**,秩揭示质量与 QRCP 相当;Martinsson 等([arXiv:2008.04447](https://arxiv.org/abs/2008.04447))同方向。2025 年有人在 H100 上做出**约为无主元 geqrf 65% 速度的定制 GPU QRCP**([arXiv:2507.00976](https://arxiv.org/abs/2507.00976))——说明技术上"能做",但都是定制代码,不是拿来即用的库。
- **领域侧的做法是绕开而不是移植**:ISDF 社区明确把 QRCP 认定为最贵步骤后,用 K-means/质心 Voronoi(CVT)聚类替代之(Dong–Hu–Lin,[arXiv:1711.01531](https://arxiv.org/abs/1711.01531)/JCTC 2018;平面波复值 K-means,[arXiv:2208.07731](https://arxiv.org/html/2208.07731);[ML K-means JPCA](https://pubs.acs.org/doi/10.1021/acs.jpca.0c06019))。这是与 SCDM 同构的"选网格点"问题的最重要前例:**成本压到可忽略靠换算法,不靠 GPU QRCP**。
- QE 的 SCDM(pw2wannier90,6.3 起)是 **MPI 并行**([arXiv:1909.00433](https://arxiv.org/abs/1909.00433));截至本次检索**未找到 pw2wannier90/SCDM 的 GPU 移植证据**(QE 的 GPU 化集中在 pw.x/CP/线性响应:[JCTC exascale](https://pubs.acs.org/doi/10.1021/acs.jctc.3c00249)、[README_GPU](https://gitlab.com/QEF/q-e/-/blob/develop/README_GPU.md))。

**HBM 判定**:V100 16 GB 上,DFT 波函数已占大头,5.3 GB 的 A 矩阵(MAGMA geqp3 要求驻卡)+ 工作区**无法与 VASP 数据共卡**;分块列主元 QRCP 需要全局列范数更新,通信结构决定了没有干净的分块流式方案(RQRCP 的 sketch 本身仍是 ~(NW+p)×Nr ≈ 4 GB 级,还是要流式扫 A)[推断]。专卡/离线可勉强放下,但——

**收益判定(一票否决)**:SCDM-only 全段实测 8.3 s/帧,QRCP 是其主体;即使 GPU 化做到 ∞ 加速,每帧上限省 ~8 s,占 w90 段 554 s 的 **1.5%**。任何工程量在这个分母面前都不成立。若未来真要压这 8 s:CPU 侧 RQRCP(Duersch–Gu)或 CVT/K-means 替代都比 GPU QRCP 便宜一个量级[推断],且 K-means 路线还顺带消掉主元离散翻转的数值风险(代价是初猜定义改变,需重新过 LEOS 验证)。

**数值风险(全管线最高)**:QRCP 主元选择是离散决策;列范数近平局时,浮点顺序差异即可翻转枢轴 → SCDM 初猜(gauge)改变 → wannierise 落入不同局部极小 → **Wannier 中心可发生非微扰级移动(病态近平局情形 0.01~0.1 Å 量级)**[推断]。这正是 LEOS stray/近平局告警体系监控的区域:等价于往闸门最敏感的输入端注入非确定性。

---

## ④ 子路线 2:disentanglement GPU 化 —— 结论:可做但只作为路线 3 的搭车项

**并行性**:每迭代 = 每 k 组装 Z(数个 NB·NW² 级 zgemm)+ 每 k 一次 ~(ndimwin−ndimfroz) ≤ NB 维厄米特征分解 + 标量混合与 Ω_I 归约。k 维 batch(Li:8;Ce:32)是唯一并行来源;单矩阵 332~352 维喂不饱 V100,必须 batched。

**库支持是短板**:
- batched zgemm:cuBLAS `cublasZgemmBatched/StridedBatched` 成熟,无问题。
- **batched 厄米特征分解(n≈300、batch 8~32)在现有工具链上没有好答案**:cuSOLVER `syevjBatched`(Jacobi)只在很小矩阵(n≲32)上占优([第三方实测,arXiv:2308.03310](https://arxiv.org/pdf/2308.03310):n=32 起就被自研实现反超;另见 [Sunny.jl 的调研](https://github.com/SunnySuite/Sunny.jl/issues/219)),对 n≈330 完全不是甜点区;新的 `cusolverDnXsyevBatched` 属于较新 CUDA(≥12.6)([cuSOLVER 13.x 文档](https://docs.nvidia.com/cuda/cusolver/index.html)),而 **nvhpc 24.3 随附 CUDA 12.3,用不上**[推断:以 nvaccelinfo/nvcc --version 核实]。可行替代:逐 k 在多 stream 上循环 `zheevd`(每个 ~3–8 ms,batch 32 → 每迭代 0.1~0.25 s,可能反而成为 GPU 瓶颈)或引入 MAGMA(新增依赖)或特征分解留在 CPU 线程池、只把 gemm 下卡(每迭代 PCIe 往返 ~30–60 MB,几 ms,可行)[均为推断]。
- 每迭代还有 Ω_I 归约、mixing、窗口逻辑等标量段,GPU 化后成为 Amdahl 尾巴。

**收益**:dis 实测 34~140 s/自旋/帧(自旋通道难度差异大)。GPU 化现实加速 5~10×[推断:gemm 部分 20~30×,被特征分解与标量段拉回],即每帧省 ~60~200 s——有肉,但与 wannierise 共用同一套工程载体,不应单独立项。

**HBM**:czmat/M/U 全集 <1 GB,无压力(共卡跑在 VASP 进程里则另说,见风险清单)。

**数值风险**:dis 的子空间选择在冻结窗口边缘近简并处是离散的(取前 NW 个特征向量);特征值近平局时 CPU/GPU 会选出不同子空间 → 下游中心移动。与 LEOS 近平局告警同区,需金标集回归兜底。

---

## ⑤ 子路线 3:wannierise GPU 化 —— 三条中唯一值得考虑,且指定工程形态

**并行性(最好)**:每迭代 = 每 k 的 cdq 组装(gemm)→ zheev(NW) 重构 exp(dW) → 3 个 NW³ gemm 更新 U → 每 (k,b) 2 个 NW³ gemm 旋转 M → 全局 spread/中心归约(小)。gemm 占比最高、结构最规整,batch = nk×(2·nntot+~4) 个 250³~288³ zgemm ≈ 30(Li)~100+(Ce) Gflop/迭代。

**V100 现实加速估计**[推断]:gemm 段 ~10–30 ms/迭代(3 TF/s 有效)对 CPU ~300 ms → 30× 上下;加上每迭代 nk 次 zheev(NW)(同路线 2 的 batched-eig 短板,~20–40 ms)与归约/传输,**wannierise 段整体 5~15×,乐观 20×;50× 在 nk≤32、NW≤300、V100 的组合下不现实**。可选优化:用纯 gemm 的 Padé/scaling-squaring 求 exp(dW) 消掉特征分解——但这改变与上游 wannier90 的逐位对齐性,验证负担转移到物理端[推断]。

**工程路径三选一(关键决策)**:

| 路径 | 说明 | 判定 |
|---|---|---|
| A. 改 wannier90 上游/fork(OpenACC/cuBLAS 进 dis/wann) | 与上游 3.1→4.0 演进永久 rebase;**与项目既定边界直接冲突:LEOS 补丁目前只动 VASP 侧,wannier90 是未修改的上游库——fork 之后论文的可复现性叙事和维护承诺都变重**(显式标注:这是边界冲突项) | 不推荐 |
| B. VASP 内重写 dis/wann | 在 mlwf.F 旁复制一份 Wannier 物理,失去 wannier90 作为参考实现的对照价值,验证负担最大 | 明确不做 |
| C. **独立 GPU 后处理工具**:VASP 置 LWRITE_MMN_AMN([实证] mlwf.F:378 已有开关)落 .mmn/.amn/.eig → 离线 batched 求解(PyTorch/CuPy complex128 batched matmul + eigh 即可,V100 FP64 7.8 TF/s)→ 输出 U/中心给 pyleos(~/LEOS/leos_analysis/pyleos 已存在,对接点待确认[推断]) | 不碰上游、不进 VASP 进程(HBM 自理)、天然支持**帧间 warm start**与多帧流水;可用 standalone wannier90.x 作金标裁判 | **唯一推荐形态** |

**工程量估计**[推断]:路径 C 原型 2~4 周 + 金标验证 2 周+ + pyleos 对接,合计 1~2 个月;路径 A 2~3 个月 + 永久维护税;路径 B ≥3 个月且不可取。

---

## ⑥ 收益账(Amdahl,两场景)

### 场景 ①:MD 轨迹逐帧分析(每帧 w90 段 554 s;50 ps/500 帧;基线 w90 段合计 277,000 s = 76.9 h)

| 方案 | 每帧 w90 段 | 每帧节省 | 500 帧节省 | 占可省上限(554 s)比 |
|---|---|---|---|---|
| GPU 5× | 110.8 s | 443.2 s | 61.6 h | 80.0% |
| GPU 20× | 27.7 s | 526.3 s | 73.1 h | 95.0% |
| GPU 50× | 11.1 s | 542.9 s | 75.4 h | 98.0% |

- **边际收益急坠**:5× 已拿走 80%,20×→50× 只再多 2.3 h/500 帧。结合 §5 的现实加速估计(5~15×),GPU 路线的期望值就是"每帧省 440~530 s"。
- **GPU·h 视角**(sai-pod 整节点独占、w90 段 16 卡全闲):基线浪费 554×16/3600 ≈ **2.46 GPU·h/帧 → 1231 GPU·h/500 帧**;任何把 w90 段搬离 GPU 节点关键路径的方案(离线后处理、GPU 化、CPU 分区跑)都能回收其中 ~95%。

**替代方案对照(同场景)**:

| 替代路线 | 每帧 w90 段(推断除注明外) | 相对基线 | 新代码量 | 备注 |
|---|---|---|---|---|
| SCDM-only 分层(已有) | 8.3 s [实证] | 67×,≈"GPU 50×"档 | 0 | 仅闸门放行的帧可用;金属/混价帧回落到迭代路径 |
| 离线 standalone w90.x MPI(k 并行)+ 两自旋并发 | Li(nk=8):~45–90 s;Ce(nk=32):可到 ~30–50 s | 6~14×(Li)/ up to ~20×(Ce) | ≈0(作业编排 + 中心回读脚本) | 上游自带能力([实证] 源码 + [社区论文](https://arxiv.org/pdf/1907.09788));加速上界受 nk 限制,Γ-only 帧无效 |
| + 帧间 warm start(用上一帧 U 作初猜/投影) | 迭代数 ~460 → 数十,叠加后 ~10–30 s | 额外 3~10× | 小工具级(路径 C 内免费获得) | AIMD 标准实践:CP2K `&LOCALIZE USE_HISTORY`([CP2K 讨论](https://groups.google.com/d/topic/cp2k/2PhCl2sCaDc)、[TRAVIS 教程:每 5 步定位一次即可](http://www.travis-analyzer.de/files/travis_ir_raman.pdf)) |
| 降低分析频率(每 5 帧) | 平均 110.8 s/帧 等效 | 5× 等效吞吐 | 0 | 需物理判断:价态涨落自相关时间 vs 0.05~0.1 ps 步距 |
| 多帧 farming(帧间并行) | 不变(延迟),吞吐 ∝ 节点数 | 吞吐任意 | 0 | 只解吞吐不解延迟与 GPU·h 成本;与以上全部可叠加 |

**场景 ① 结论**:GPU 化(现实 5~15×)拿到的每帧节省,与"离线 MPI + 自旋并发 + warm start"的组合(合成 10~50×)**同一量级**;后者零边界冲突、零数值新风险。GPU 路线只在以下交集才有独立价值:**帧必须走完整迭代(SCDM-only 被闸门拦)且 nk 较大(≥16)且 CPU 组合拳已部署仍不够用且在意单帧延迟而非吞吐**。

### 场景 ②:单点批跑基准(SCF 占大头,w90 段占比 f)

| f(w90 段占单点总时) | GPU 20× 后整体加速 | GPU ∞ 后整体加速上限 |
|---|---|---|
| 5% | 1.05× | 1.05× |
| 10% | 1.10× | 1.11× |
| 15% | 1.17× | 1.18× |

**场景 ② 结论**:上限贡献 ≈ f 本身,个位数百分点。**单点批跑场景下 GPU 化 w90 段没有立项价值**;该场景的优化主战场在 SCF/内存(另一 agent 的内存账与 gpu-map 范围)。

---

## ⑦ 前例调研结果(找到什么 / 没找到什么)

**找到(可借鉴/可引用)**:
- cuSOLVER 无 geqp3;GPU 侧 QRCP 现状 = MAGMA 混合实现或自研:[cuSOLVER docs](https://docs.nvidia.com/cuda/cusolver/index.html)、[MAGMA geqp3](https://icl.utk.edu/projectsfiles/magma/doxygen/group__magma__geqp3.html)、[PyTorch #82092](https://github.com/pytorch/pytorch/issues/82092)、[JAX #25955](https://github.com/jax-ml/jax/pull/25955)、[QRCP 受内存带宽限制的分析](https://www.cs.ucdavis.edu/~bai/publications/tomasbaihernandez13.pdf)。
- 随机化 QRCP 路线成熟(CPU 为主,GPU 定制可达无主元 QR 的 ~65%):[Duersch–Gu arXiv:1509.06820](https://arxiv.org/abs/1509.06820)、[arXiv:2008.04447](https://arxiv.org/abs/2008.04447)、[arXiv:2507.00976(H100,2025)](https://arxiv.org/abs/2507.00976)。
- 同构问题的领域解法是"换算法绕开 QRCP":ISDF 的 CVT/K-means 系列 [arXiv:1711.01531](https://arxiv.org/abs/1711.01531)、[arXiv:2208.07731](https://arxiv.org/html/2208.07731)、[JPCA ML-K-means](https://pubs.acs.org/doi/10.1021/acs.jpca.0c06019)。
- **MLWF 极小化跑在 GPU 上的唯一近前例:JDFTx**(整代码 CUDA 设计,自带 wannier 模块与 CG/L-BFGS 极小化;changelog 提到通过减少矩阵乘与通信显著加速 Wannier 极小化):[wannier-minimize 文档](http://jdftx.org/CommandWannierMinimize.html)、[ChangeLog](http://jdftx.org/ChangeLog.html)。证明数学上就是 batched 稠密线代、GPU 可行;但它是另一个代码体系,无法直接复用。
- Wannier 生态的 GPU 化集中在**消费端(插值)而非生产端(wannier 化)**:EPW exascale(MPI-GPU-OpenMP 插值,至 29×)[arXiv:2603.10295 / npj CM 2026](https://arxiv.org/abs/2603.10295)、Perturbo GPU(OpenACC,~40×)[arXiv:2511.03683](https://arxiv.org/html/2511.03683v1)、另有电声插值 GPU 实现 [arXiv:2306.16493](https://arxiv.org/abs/2306.16493)。这佐证:社区认为 wannier 化本身不是值得 GPU 化的热点。
- wannier90 官方路线图:v4 重做**库模式使之可并行(MPI)调用、免文件传参**——不是 GPU:[Wannier 生态综述 arXiv:2312.10769](https://arxiv.org/pdf/2312.10769)、[社区论文 arXiv:1907.09788](https://arxiv.org/pdf/1907.09788)。**等 v4 落地,"库模式串行"这个当前痛点会被上游自己解决**;现在 fork 是在和上游赛跑。
- CP2K/AIMD 的帧间 warm start 与降频实践:[USE_HISTORY 讨论](https://groups.google.com/d/topic/cp2k/2PhCl2sCaDc)、[TRAVIS 教程](http://www.travis-analyzer.de/files/travis_ir_raman.pdf)、[CP2K wannier 练习](https://www.cp2k.org/exercises:2015_ethz_mmm:wannier)。

**没找到(检索为空,截至 2026-07)**:
- 任何 wannier90 的公开 GPU 分支/fork/论文(官方 [features](https://wannier.org/features/)、[GitHub issues](https://github.com/wannier-developers/wannier90/issues) 检索无果)。
- pw2wannier90 的 SCDM/QRCP GPU 移植证据(QE GPU 文档与 exascale 论文均未提及)。
- cuSOLVER/SLATE 等主流库中的 batched 大尺寸(n~300)厄米特征分解甜点方案(syevjBatched 只适于 n≲32,[实测来源](https://arxiv.org/pdf/2308.03310);XsyevBatched 需较新 CUDA)。

---

## ⑧ 风险清单

1. **数值一致性(最高优先)**:Wannier 中心是 LEOS 的物理输入。三个离散决策点对浮点顺序敏感——SCDM 枢轴选择(路线 1)、dis 窗口边缘近简并子空间选取(路线 2)、wannierise 非凸多极小(路线 3)。GPU/CPU 的归约顺序差异(每 gemm 相对 ~1e-14)在**收敛到同一极小**时对中心的扰动 ≲1e-6 Å(无害)[推断];但在近平局区可放大为不同离散选择 → 中心移动 0.01~0.1 Å 量级 → 直接触发/污染 stray 与近平局告警统计。**兜底要求**:任何 GPU 实现必须带 (a) CPU wannier90 金标裁判模式,(b) 金标集中心回归(阈值建议 1e-4 Å),(c) 告警率 A/B 对比;cuBLAS 同版本同形状 run-to-run 确定,但 CPU↔GPU 永不逐位一致,不要以逐位一致为验收标准。
2. **HBM(硬约束)**:路线 2/3 工作集 ≤2 GB,单卡可容;但**若在 VASP 进程内共卡实现,V100 16 GB 已被 DFT 波函数占用大半(题设),OOM 风险实际存在且随体系增大恶化**——这是选择"独立后处理工具(专卡/专节点)"形态的第二理由。路线 1 的 5~10 GB 矩阵在共卡场景直接不可行。
3. **工具链混用**:VASP OpenACC(nvhpc 24.3,CUDA 12.3)+ cuBLAS/cuSOLVER 本身兼容(VASP GPU 版已链);但 (a) 需要的 `cusolverDnXsyevBatched` 在 CUDA 12.3 不存在,(b) 引入 MAGMA = 新增一个须与 CUDA 12.3/cc70 对齐编译的依赖,(c) wannier90 库须与 nvfortran 工具链共编译,fork 后每次上游同步都要重验。独立 Python 工具(PyTorch/CuPy)把这些全部隔离在 VASP 工具链之外。
4. **"改 wannier90 上游"与项目既定边界的冲突(显式标注)**:LEOS 的补丁边界是 VASP 侧,wannier90 以未修改的上游 3.x 库形态被调用——这是论文可复现性与维护承诺的一部分。路径 A(fork/改上游)直接违反该边界;且上游 v4 正在重做并行库模式,fork 的生命周期价值存疑。若未来确需上游能力,正确姿势是向 wannier90 社区提 issue/PR 而不是私有 fork。
5. **可维护性/人力**:三条路线共用的真实工程是"batched 稠密线代 + 收敛逻辑移植 + 金标验证体系"。以 fork/VASP 内重写形态做,维护税永久;以独立工具形态做,风险被限制在一个可丢弃的仓库里。
6. **收益不确定性**:§2.2 提示当前 554 s 基线可能已含线程 BLAS 效应(待核实 OMP 设置);若核实后基线可通过纯环境变量(线程 BLAS)再降 2~4×,GPU 路线的相对收益进一步缩水。另:Γ-only 或小 nk 的轨迹帧同时削弱 k-batch(GPU)与 k-MPI(CPU)两条路,此时唯一有效杠杆是 warm start + 自旋并发 + 降频。
7. **机会成本**:同一份工程人力投在"SCF 段 GPU 内存账/执行图优化"(另两 agent 的范围)或"warm start 工具"上,单位收益更高;w90 段 GPU 化应排在它们之后。

---

## 附:本次未做/待核实清单

- 未实测 GPU 原型(纯纸面评估);§3–§5 的加速比区间均为 [推断]。
- 待核实:w90 库链接的 BLAS 实现与运行时 OMP_NUM_THREADS(影响基线);nvhpc 24.3 实际随附 CUDA 版本;Ce 体系 w90 段实测耗时(本报告用 Li2MnO3 锚定);pyleos 离线摄入 wannier90 standalone 输出(.wout/chk 中心)的现成程度;ws_bench .wout 内两次运行的归属(两帧或重复调用)。
- 未做:未触碰集群文件写操作,未提交任何作业,未与 gpu-map/内存账 agent 的产出重复。
