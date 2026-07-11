# gpuw90 — GPU Wannier-化后处理工具 · 设计文档 (Phase 1 + Phase 2)

独立于 pyleos 的后处理工具。读 wannier90 标准格式(`.mmn`/`.amn`/`.eig`/`.win`,可选 `.nnkp`),
用 PyTorch 批量实现 Souza–Marzari–Vanderbilt 解纠缠(disentangle)与 Marzari–Vanderbilt(MV)
幺正最小化(wannierise),输出 Wannier 中心与 spread 供 LEOS 氧化态管线消费。
**不 fork wannier90、不改 VASP**;device 参数化,CPU/CUDA 同码。

Phase 1 = **wannierise**(孤立带组,`num_bands == num_wann`,已完成)。
Phase 2 = **disentangle**(纠缠带组 `nb > nw`,本次完成;金标见 §5b,实现细节见 §7)。

---

## 1. 范围与阶段边界

| | Phase 1(已完成) | Phase 2(本次完成) |
|---|---|---|
| 处理对象 | 孤立带组 `nb == nw` | 纠缠带组 `nb > nw`(金属/含导带) |
| 核心运算 | 每 k:1 个 `nw×nw` 厄米特征分解 + ~20 个 `nw³` matmul | 追加每迭代每 k 一个 `nb` 维厄米特征分解(Z 矩阵子空间选取) |
| 输入 | `.mmn`(`nb×nb` 重叠)+ `.amn`(`nb×nw` 初猜)+ `.win` | 追加 `.eig`(能带能量,定冻结/外窗) |
| 优化器 | 最速下降(SD)+ 固定步长 | dis 为 Z 矩阵不动点迭代(w90 同款 mixing);CG+线搜索仍留待 Phase 2+ |

设计原则(已兑现):disentangle 作为 wannierise 的**前处理**,产出降维后的 `M0_reduced [nk,nnb,nw,nw]`
与初始 `U0 [nk,nw,nw]`,直接喂给已验证的 `wannierise()`,后者数学不改(本次仅加防错守卫
与热循环去同步,金标数值逐位复现,见 §5)。

---

## 2. 模块与数据流

```
        VASP (LWRITE_MMN_AMN)                 wannier90.x (参考金标, 仅验证期)
                 │                                        │
      seedname.{mmn,amn,eig,win}                    seedname.wout
                 │                                        │
   ┌─────────────▼─────────────┐                          │
   │ io_w90.py                 │                          │
   │  read_mmn / read_amn      │                          │
   │  read_eig(纠缠时必需)      │                          │
   │  read_win / read_nnkp     │                          │
   │  compute_bvectors         │  b 向量 (Å⁻¹)             │
   │  compute_shell_weights    │  w_b (B1 完备关系)         │
   │  → W90Case(nb>nw 需显式    │                          │
   │    allow_entangled,否则拒载)│                          │
   └─────────────┬─────────────┘                          │
                 │  nb>nw?── 是 ──┐                        │
                 │ 否             │  M0_full, A_full, eig   │
                 │   ┌────────────▼────────────┐           │
                 │   │ disentangle.py (torch)  │           │
                 │   │  窗口/冻结掩码 → 初猜投影  │           │
                 │   │  Z 矩阵不动点迭代(batched │           │
                 │   │  eigh + mixing,§7)      │           │
                 │   │  → M0_reduced, U0        │           │
                 │   └────────────┬────────────┘           │
                 │                │                        │
   ┌─────────────▼────────────────▼──────┐                 │
   │ wannierise.py (torch)               │                 │
   │  lowdin_orthonormalize → U0(孤立带)  │                 │
   │  wannierise():                      │                 │
   │   循环: rotate M → 梯度 g            │                 │
   │        → exp(dW) → 更新 U            │                 │
   │  → centres[nw,3], spreads           │                 │
   └─────────────┬───────────────────────┘                 │
                 │                          ┌─────────────▼──────────────┐
                 │                          │ io_w90.parse_wout_final_…  │
                 │                          │ io_w90.parse_wout_dis_info │
                 │                          └─────────────┬──────────────┘
   ┌─────────────▼──────────────────────────────────────▼─────────────┐
   │ check_against_w90.py(孤立带)/ check_dis_against_w90.py(纠缠带)     │
   │   最近邻(允许简并交换)+ 最小镜像距离对照 · PASS/FAIL(1e-4 Å)         │
   └───────────────────────────────────────────────────────────────────┘
```

### 批量结构(k = 天然 batch 维)
- `M0`  : `[nk, nnb, nw, nw]` complex128 —— 每 (k, 近邻 b) 的 ab-initio 重叠
- `U`   : `[nk, nw, nw]` complex128 —— 每 k 的幺正规范
- `bvec`: `[nk, nnb, 3]` float64 —— 笛卡尔 b 向量(Å⁻¹)
- `wb`  : `[nk, nnb]` float64 —— shell 权重
- `kb_idx`: `[nk, nnb]` long —— 近邻 k 索引(k+b 落回 BZ 的那个 k)

所有重活都是 batched matmul(`torch.matmul` 广播)与 batched 厄米特征分解(`torch.linalg.eigh`),
对应 wannier90 里逐 k 的 `zgemm`/`zheev`。切换 `device="cuda"` 即上卡,数学不变。

---

## 3. 算法(MV 最速下降,wannier90 3.1 同款)

单步(全 k 批量):

1. **旋转重叠** `M^{k,b} = U_k† M0^{k,b} U_{k+b}`(`_rotate_M`,近邻用 `U[kb_idx]` 取)
2. **相位** `Im ln M_nn = atan2(Im, Re)`(主支)
3. **中心** `r_n = -(1/N) Σ_{k,b} w_b b · Im ln M_nn`
4. **梯度**(真实 `dΩ/dW`,反厄米)
   `g_k = -(1/N) · 4 Σ_b w_b ( A[R^{k,b}] − S[T^{k,b}] )`
   其中 `R_mn = M_mn·conj(M_nn)`,`T_mn = (M_mn/M_nn)·q_n`,`q_n = Im ln M_nn + b·r_n`,
   `A[B]=(B−B†)/2`,`S[B]=(B+B†)/(2i)`
5. **步长** `dW_k = −(α·N/(4 w_tot)) g_k`(α=1 即经典 MV 步长;`w_tot=Σ_b w_b`)
6. **更新** `U_k ← U_k · exp(dW_k)`,矩阵指数经 `i·dW` 的厄米特征分解:
   `exp(dW) = V·diag(e^{−iλ})·V†`(即 wannier90 的 `zheev` 求指数)
7. **收敛** `|ΔΩ| < conv_tol` 连续 `conv_window` 次(与 wannier90 语义一致)

初始规范:`U0 = A(A†A)^{-1/2}`(Löwdin 对称正交化)—— 对应 wannier90 的 "Initial State"。

---

## 4. 数值坑与处理(逐项)

> 老大特别关照的部分。每个坑都在代码里就地注释。

1. **梯度的 1/N 前因子与整体符号(最隐蔽,已用有限差分钉死)**
   MV 梯度常被写成 `4 Σ_b w_b (A[R]−S[T])`,但 Ω 自带 `1/N` 求和,且实现里易差一个整体符号。
   用有限差分方向导数检验(`fd_check.py`):比较解析 `Re Tr[g† X]` 与
   `[Ω(U e^{εX})−Ω(U e^{−εX})]/2ε`,修正前 ratio = **−1/8 = −1/N**,修正后 = **+1.0000**(三个 ε 稳定)。
   结论:真实梯度 `g = −(1/N)·4 Σ_b w_b (A[R]−S[T])`,标准下降 `dW=−eff·g` 单调降 Ω。
   **教训:任何梯度实现落地前先过 FD 检验,不要相信纸面公式的前因子/符号。**

2. **Ω_I 与 Ω_OD 的拆分(易双重计数)**
   `Ω_I = (1/N)Σ w_b (nw − Σ_{mn}|M_mn|²)`(对全部 m,n 求和),
   `Ω_OD = (1/N)Σ w_b Σ_{m≠n}|M_mn|²`。
   初版误把 `Σ_n|M_nn|²` 放进 Ω_I,相当于把 Ω_OD 折进 Ω_I 后又加了一遍 → 总 Ω 偏高恰好一个 Ω_OD。
   现象:中心全对(6e-7)但总 spread 偏 0.502(=Ω_OD)。修正后初态 Ω 与金标逐位一致(4.46881212)。

3. **ln 支割线(主支)**
   `Im ln M_nn = atan2(Im, Re) ∈ (−π, π]`,与 wannier90 的 `aimag(log(...))` 一致;CPU/GPU 同一主值,
   无跨设备支跳变。绝缘体 WF 的 `M_nn ≈ 1`(相位小),远离 ±π,安全。近平局或 `|M_nn|→0` 时才有支翻转风险
   (会使中心跳 `b·2π/…`);本阶段金标 O_D 极小,未触发。

4. **b 向量权重(shell 最小二乘 + B1 校验)**
   按 wannier90 kmesh 规则:b 向量按 `|b|` 分 shell,解最小二乘
   `Σ_shell w_s (Σ_{b∈s} b_a b_b) = δ_ab` 的 6 个独立对称分量。均匀网格上各 k 共享 shell 结构,权重按 `|b|` 匹配回填。
   落地即校验 B1 残差 `‖Σ_b w_b b⊗b − I‖`:GaAs 9e-18、Diamond 2e-16,且 w_b 与金标逐位吻合(0.408635 / 0.527876)。

5. **简并 eigh 的相位规范(实为伪坑)**
   `i·dW` 的特征向量在简并子空间内只定到一个幺正旋转/相位。但 `exp(dW)=V·f(λ)·V†` 是算符函数,
   对该规范**不变**——所以"简并 eigh 相位"不影响更新。同理 Löwdin 的 `(A†A)^{-1/2}` 也规范无关。
   代码就地注明,避免后人误加相位固定逻辑。

6. **FP64 全程**
   torch complex128 / float64,对齐 wannier90 的 FP64。CPU↔GPU 只在归约顺序上差 ~1e-14,
   收敛到同一极小时对中心的扰动 ≪ 1e-4 Å(见 §5 实测 ~1e-6~1e-9)。不追求逐位一致(物理等价即可)。

7. **WF 顺序与周期镜像(对照期)**
   最小化可能给出与金标不同的 WF 排序或近简并交换。对照时按**最近中心**贪心匹配(非按索引),
   距离用**最小镜像**(经实空间晶格折叠),从而容忍简并交换与晶格平移镜像。验收阈值 1e-4 Å。

8. **`.mmn` 列主序**
   块内 `nb×nb` 数据以带序 m 为最快索引(Fortran 列主序),`reshape(order="F")` 得 `M[m,n]`。
   诊断办法:若误当行主序,中心的对称性会破而 B1 仍过——本工具用初态中心与金标比对可当场抓出。

9. **(Phase 2)金标参考本身未收敛的坑**
   站点 wout 的 wannierise 段是 CG、100 迭代封顶、无收敛判据;若拿参考"末态中心"当验收目标,
   会把协议差(实测 ~0.3 Å)误判成实现错误。教训:**先读参考跑的参数回显再定金标口径**;
   纠缠阶段的硬验收只用协议无关量(Ω_I + 解纠缠初态),见 §5b。

10. **(Phase 2)双 Ω_I 的时序**
   `Ω_I(i−1)` 的冻结贡献必须用**上一迭代**的邻居子空间(在更新任何 k 之前整体计算),
   `Ω_I(i)` 才用全新子空间;`Z_out` 在迭代尾部、更新之后组装。次序错一处,delta 序列就偏离
   w90,收敛迭代数对不上(对齐后四个金标通道迭代数逐一相同)。

11. **(Phase 2)能级移位的幅度选择**
   冻结/窗外方向的对角移位若取"大数"(如 1e6),eigh 的特征向量误差 ~eps·‖A‖/gap 会被放大进
   自由子空间。取与谱同量级的 `±2·w_tot`(Z 谱有界 `[0, w_tot]`,因 `‖M‖≤1`)即可保证块间分离,
   又不牺牲精度。同理初猜 QPQ(谱 `[0,1]`)用 ±2。

12. **热循环的隐性 host↔device 同步(本次修复)**
   旧 wannierise 每迭代 5 次 `float()`(rms_grad + 4 个 Ω 分量)各强制一次设备同步;GPU 上这是
   每迭代的固定延迟税。现收敛计数(`n_below`)、Ω 历史全部留在设备端张量,每迭代只 poll 一个
   标量收敛旗标(`check_every` 可再稀释);顺带复用更新后的旋转 M(上一版每迭代重复旋转一次,
   数学恒等)。重构后 GaAs/Diamond 金标数字逐位不变。dis 循环同款设计。

---

## 5. 正确性验证(金标对照结果)

金标 = wannier90 3.1 test-suite 的 `benchmark.out`(完全收敛:GaAs 第 2 迭代起中心固定不变)。

| 算例 | nb/nw | nk | nntot | 收敛迭代 | **最大中心偏差 (Å)** | 最大 spread 偏差 (Ų) | Ω_total 偏差 |
|---|---|---|---|---|---|---|---|
| GaAs (example01) | 4/4 | 8 | 8 | 10 | **8.31e-07** | 4.07e-09 | 2.6e-10 |
| Diamond (example05) | 4/4 | 64 | 8 | 12 | **3.13e-09** | 1.98e-09 | 1.3e-10 |

两例均 **PASS**(阈值 1e-4 Å,余量 10⁵~10⁷ 倍)。Diamond 的 64 k 点验证了批量 k 路径。
初态(Löwdin 后未优化)中心与金标 "Initial State" 也吻合到 6e-7 / 3e-9 Å,单独验证了 io + 权重 + 中心公式。

复现:`python check_against_w90.py testdata/gaas --seed gaas`(详见 README.md)。
(Phase 2 的热循环去同步重构后复跑,两例数字逐位不变:8.31e-07 / 3.13e-09 Å。)

### 5b. Phase 2 解纠缠金标(带冻结窗的真实两自旋体系)

金标 = sai-pod 站点 wannier90.x 3.1 串行跑(`run_np1_s*/wannier90.*.wout`),与本工具读**同一套**
VASP 导出 `.mmn/.amn/.eig/.win`(符号链接同源)。两套算例都带冻结窗(`dis_froz_max` 来自 .win),
外窗为默认全范围;li2mno3 是 LEOS 生产尺度(288 WF)。

| 算例(自旋) | nb/nw | nk/nnb | nf 范围 | dis 迭代 (本/金标) | **Ω_I (Ų): 两边打印值相同** | \|ΔΩ_I\| | 初态最大中心偏差 (Å) | 初态最大 spread 偏差 (Ų) |
|---|---|---|---|---|---|---|---|---|
| vircef s1 | 180/142 | 8/12 | 141 | 82 / 82(收敛) | 81.06033903 | 3.4e-09 | 7.9e-07 | 5.0e-09 |
| vircef s2 | 180/142 | 8/12 | 140 | 188 / 188(收敛) | 80.32720233 | 4.6e-10 | 7.7e-07 | 4.9e-09 |
| li2mno3 s1 | 356/288 | 8/8 | 284~286 | 38 / 38(收敛) | 154.86511542 | 4.4e-09 | 8.3e-07 | 5.0e-09 |
| li2mno3 s2 | 356/288 | 8/8 | 240 | 200 / 200(**双双到 dis_num_iter 上限未收敛,协议对齐复现**) | 135.30457827 | 4.1e-09 | 8.2e-07 | 5.0e-09 |

四通道全 **PASS**(门:\|ΔΩ_I\| ≤ 1e-4 Ų、初态中心 ≤ 1e-4 Å)。要点:

- **验收量选择**:金标 wout 的 wannierise 是 CG、固定 100 迭代、**无收敛判据**(conv_window=−1),
  四个参考末态全都未收敛(末端 ΔΩ 还在 1e-2~1e-1 Ų/迭代)——末态中心对比只反映协议差
  (实测 ~0.3 Å,且我们收敛端 Ω_total 更低:vircef s1 94.2618 vs 参考 94.6541),不是正确性信号。
  因此硬验收采用**协议无关量**:Ω_I(dis 的 gauge 不变输出)+ 解纠缠"Initial State" 中心/spread
  (子空间 + Löwdin 投影初值,确定性,不经任何最小化)。这正是"以 Ω_I 和中心集合匹配为准"。
- 初态中心偏差 ~8e-7 Å、spread 偏差 ~5e-9 Ų 是 **wout 打印精度地板**(中心 6 位小数、spread
  8 位小数的舍入),与 Phase 1 GaAs 的 8.3e-7 同源——真实吻合度更高。
- **迭代数逐一相同**(38/82/188/200)= mixing、双 Ω_I、收敛窗口语义全部对齐的强证据;
  li2mno3 s2 在 200 迭代封顶的**未收敛轨迹**上仍复现到 4e-9 Ų,说明不动点映射的浮点漂移可忽略。
- `Ω_I(初态) − Ω_I(dis)` 一致性(gauge 不变量,~1e-11)在脚本里作为内置 sanity 每次自查。

复现:`python check_dis_against_w90.py testdata/vircef --seed wannier90.1`(四通道换 datadir/--seed;
大输入不入库,先跑 `testdata/fetch_entangled.sh`)。

### 5c. 生产尺度计时:本地 CPU 基线(2026-07-11 实测)

环境:Apple M 系笔记本,torch 2.8.0 CPU 轮子,FP64 全程。对照 = sai-pod 站点串行 wannier90.x
(`run_np1_s1` wout 的计时,Xeon/AMD 单进程)。

| 段 | 体系 | 本工具 CPU (s/iter) | 站点 wannier90 串行 (s/iter) | 备注 |
|---|---|---|---|---|
| disentangle | li2mno3 288 WF | **0.24**(38 迭代 9.1 s) | 1.10(38 迭代 41.9 s) | 同迭代数,4.6× |
| wannierise  | li2mno3 288 WF | **0.196** | ~0.92(CG,100 迭代 ~92 s) | 每迭代工作量同型 |
| disentangle | vircef 142 WF | 0.07 | ~0.21(82 迭代 17.1 s) | |
| wannierise  | vircef 142 WF | 0.051 | ~0.18 | |
| `.mmn` 读入 | li2mno3(258 MB 文本) | 3.4 s | 4.7 s | 单次成本 |

(本工具快的来源是多线程 BLAS 的 batched matmul;这也是可行性报告"GPU 化对比基线要抬高"的
现身说法——GPU 加速比要对着这条 CPU torch 基线量,而不是站点串行 Fortran。)
**CUDA 实测待集群空档**:`bench_gpu.py` + `run_bench_v100.sbatch` 就绪,先金标门后计时;
V100 环境坑(cu128 轮子无 sm_70、16V100 分区 ≥4 GPU/作业)已写进模板注释与 §9.3。

---

## 6. HBM / 显存分块策略(V100-SXM2-16GB)

本阶段工作集很小,**HBM 不是约束**(与可行性报告 §2.1 一致:dis/wann 全集 <2 GB)。
以 Ce 体系(nk=32, nw≈253, nntot≈8~12)估:`M0 = 32×10×253²×16 B ≈ 0.33 GB`;`U ≈ 0.1 GB`。单卡绰绰有余。

给未来更大体系预留的两级分块:
- **自旋分块(ISPIN=2)**:两自旋通道完全独立 → 外层顺序处理,或分派到两块 GPU 并发。天然切分,零耦合。
- **k 分块(仅超大体系需要)**:`M0` 按 k 切 tile 从 host 流入;但 `M` 旋转要用近邻 `U[k+b]`,故
  **U 全量常驻卡**(小),只流式换入 `M0[k_tile]` 及其近邻块。梯度对 b 求和、中心归约都是逐 k 局部
  (仅近邻 U 查表跨 tile),因此 k 分块干净可行。
- 说明:真正的显存猛兽是 SCDM 的 `nb×Nr ≈ 5.3 GB` 矩阵——**明确不在本工具范围**(本工具从 VASP 落地的
  `.amn` 读初猜,不做 SCDM)。

---

## 7. Phase 2 (disentangle) — 已实现(disentangle.py)

接口按预留设计落地(返回值升级为 dataclass,多带出诊断量):

```python
res = disentangle(M0_full, A_full, eig, bvec, wb, kb_idx,
                  dis_win_min=None, dis_win_max=None,        # 缺省 = 全能带范围(w90 默认)
                  froz_win_min=None, froz_win_max=None,      # froz_min 缺省 = win_min(w90 默认)
                  num_iter=200, mix_ratio=0.5,
                  conv_tol=1e-10, conv_window=3, device="cpu")
# res.M0_reduced [nk,nnb,nw,nw] / res.U0 [nk,nw,nw] → 原样喂 wannierise()
# res.omega_I / res.history / res.converged / res.ndimwin / res.ndimfroz
# 便捷入口:disentangle_case(case)——窗口/迭代参数自动取自 .win(w90 默认值兜底)
```

与 wannier90 3.1 `disentangle.F90` 的逐步对齐(源码级核对,非仅论文公式):

1. **窗口(dis_windows)**:外窗选 `ndimwin(k)` 条带(能量升序保证连续);冻结窗 `ndimfroz(k)` 条,
   必须整条保留。`ndimwin < nw` / `ndimfroz > nw` 与 w90 同语义报错。
2. **初猜(dis_project + dis_proj_froz)**:`V0 = polar(A_w)`(batched SVD 幺正因子 = Löwdin);
   有冻结窗时非冻结部分取 `Q_f P_s Q_f` 的前 `nw−nf` 特征向量(SMV Eq. 27),冻结方向填单位向量。
3. **dis_extract 主循环(Jacobi 式,全 k 同步更新 = w90 的 MPI k 分布语义)**:
   `Z_in⁽¹⁾ = Z[V⁽⁰⁾]`;`Z_in⁽ⁱ⁾ = mix·Z_out⁽ⁱ⁻¹⁾ + (1−mix)·Z_in⁽ⁱ⁻¹⁾`;
   `V⁽ⁱ⁾ =` 冻结态 ∪ `Z_in⁽ⁱ⁾` 前 `nw−nf` 特征向量;`Z_out⁽ⁱ⁾ = Z[V⁽ⁱ⁾]`,
   其中 `Z[V](k) = Σ_b w_b M^{k,b} P_{k+b} M^{k,b†}`,`P = V V†`,再投影到 k 点非冻结补空间。
   **mixing 混的是 Z 矩阵,不是投影子**;`Z_out` 在更新后、迭代尾部组装——两处都照抄 w90。
4. **双 Ω_I 与收敛判据**:`Ω_I(i−1)`(SMV Eq. 18,邻居用上一迭代子空间 = 冻结贡献 + 选中特征值和)
   与 `Ω_I(i)`(Eq. 12,全部用新子空间);`delta = Ω_I(i−1)/Ω_I(i) − 1`;
   连续 `conv_window` 次 `|delta| < conv_tol` 收敛——四个金标通道的**迭代数与 w90 逐一相同**(§5b)。
5. **U0(internal_find_u)**:`U0 = polar(V† A_w)`(SVD 幺正因子),即 w90 的 Löwdin 投影初值。

### batched 化的关键技巧(不规则维度 → 单个 batched eigh)

`ndimwin/ndimfroz` 逐 k 不同(li2mno3 实测 `nf = 284~286` 随 k 变),按 k 循环会碎掉 batch。
处理:保留全 `nb` 维,对 Z 先零掉冻结/窗外行列,再打**对角能级移位**——冻结方向 `+2·w_tot`、
窗外方向 `−2·w_tot`(Z 谱天然落在 `[0, w_tot]`,因 `‖M‖≤1`、`w_b>0`)。于是升序特征值排布为
`[窗外 | 自由谱 | 冻结]`,**一次 `torch.linalg.eigh` 的末 nw 列就是新子空间**(冻结块是精确的
`2w_tot·I`,其特征向量恰为冻结单位向量的幺正组合——子空间层面精确,列内 gauge 无物理效应);
`Ω_I(i−1)` 的"选中特征值和"用逐 k 槽位掩码(`nb−nw ≤ pos < nb−nf(k)`)向量化。移位烘焙进
`Z_in/Z_out`,mixing 线性保持移位,槽位选择全程有效。eigh 后按列类型掩码擦除 ~1e−15 的跨块
浮点噪声(`_scrub_columns`),冻结态严格冻结。初猜的 `QPQ`(谱在 `[0,1]`,移位 ±2)同款处理,
顺带**结构性复现了 w90 的 ortho-fix**(近零特征值时特征向量须与冻结态正交——我们的补空间块
精确解耦,天然满足)。

### 刻意偏离 w90 的一处(纯 gauge,零物理差)

w90 在 dis 收敛后把 V 右乘子空间哈密顿量的特征向量(供能带插值用)。右旋 `V→VR` 在
`M0_reduced = V†M0V'` 与 `U0 = polar(V†A)` 之间**严格相消**(`polar(R†B) = R†·polar(B)`),
所有物理量(Ω_I、中心、spread、wannierise 轨迹)不变——本工具不做能带插值,故跳过。
金标初态八位小数吻合(§5b)实证了这一点。

数值风险提示不变(可行性报告 §4/§8):冻结窗边缘近简并子空间选取是离散决策,对浮点顺序敏感 →
下游中心可移动。兜底 = 本金标回归 + LEOS 近平局告警 A/B 对比。

---

## 8. 与 LEOS 管线的对接点

- **输入侧(生产路径)**:VASP 置 `LWRITE_MMN_AMN`(mlwf.F 已有开关)落 `seedname.{mmn,amn,eig}`,
  外加一个 `.win`(VASP 侧或手写,提供 cell 与显式 kpoints)。本工具**离线**读取,不进 VASP 进程(HBM 自理)。
- **输出侧(两种,推荐先 b)**:
  - (a) **EOSCAR 适配**:把 `centres[nw,3]`+`spreads[nw]`(按自旋)写成 eos.F 认得的中心表格式,直接顶替
    wannier90 库回读路径。需对齐 EOSCAR 列约定(留待与 pyleos 对接时确认)。
  - (b) **中心表直供(原型首选)**:输出 `.npz`/纯文本(每自旋 `centres`、`spreads`、`U`),由 pyleos 侧
    薄适配器摄入。解耦、可回归、不动 eos.F。
- **帧间 warm start(免费获得)**:`wannierise(U0=...)` 以上一帧收敛的 `U` 作初值即可,天然支持
  AIMD 逐帧的 warm start(可行性报告场景①的关键杠杆)。多帧可流水/farming。

---

## 9. 待办清单(Phase 2 后更新)

1. ~~**disentangle 阶段**~~ **已完成**(§7/§5b):`Z_k` 特征分解 + 冻结/外窗 + `read_eig`;
   金标未用 lead/silicon,直接上了两套带冻结窗的真实两自旋体系(vircef + li2mno3 生产尺度)。
2. **CG 优化器** + 线搜索(Fletcher–Reeves / Polak–Ribière),替换固定步长 SD。Phase 2 金标进一步
   佐证必要性:这类近简并体系上 SD 2000 迭代仍未到 1e-10(站点 CG 100 迭代也没到),良态体系无碍。
3. **CUDA 实测**:脚本与作业模板已就绪(`bench_gpu.py` + `run_bench_v100.sbatch`,含金标先行门),
   等 sai-pod GPU 空档提交。已知坑:站点 mace 环境的 torch 2.9.1+cu128 **无 sm_70**,V100 需
   cu124 系轮子(模板内有 assert 守卫);16V100 分区最少 4 GPU/作业(用 1 空 3,flood-gpu 4h)。
   CPU 基线已录入 §5c。
4. **EOSCAR/pyleos 输出适配器** + 一个 Ce 帧端到端集成测试。
5. **更大金标集**:Ce 真实帧回归(li2mno3 已入金标);大 `nw` 时把贪心中心匹配换 Hungarian。
6. **鲁棒性**:`.win` 无显式 kpoints 时的 MP 网格生成(需核实与 `.mmn` 一致的 k 排序)或强制走 `.nnkp`;
   `guiding_centres` 相位选项(棘手初猜时)。
7. **性能**:热循环已消掉重复旋转与逐迭代同步(§4.12);若剖析仍显示值得,再上增量式 M 更新
   (wannier90 风格)与 dis 大迭代的 `M0w` 常驻/分块策略。

---

## 附:文件清单

```
gpuw90/
  io_w90.py                # .mmn/.amn/.eig/.win/.nnkp 解析 + b 向量/权重重建 + .wout 参考解析
                           #   nb>nw 需显式 allow_entangled(防静默错误结果守卫)
  wannierise.py            # torch 批量 MV 最小化(device 参数化;非方阵输入守卫)
  disentangle.py           # torch 批量 SMV 解纠缠(Phase 2,§7)
  check_against_w90.py     # 孤立带金标对照(Phase 1)
  check_dis_against_w90.py # 纠缠带金标对照:Ω_I + 解纠缠初态两道硬门(Phase 2,§5b)
  test_guards.py           # 纠缠守卫单测(python test_guards.py / pytest)
  bench_gpu.py             # 固定迭代数的 dis+wann 计时(CPU 基线 / CUDA 对比)
  run_bench_v100.sbatch    # sai-pod 16V100 基准作业模板(金标先行,含环境坑注释)
  fd_check.py              # 梯度符号/前因子的有限差分自检
  DESIGN.md / README.md
  testdata/
    gaas/     gaas.{win,mmn,amn} + gaas_ref.wout            (example01, 8 k)
    diamond/  diamond.{win,mmn,amn} + diamond_ref.wout      (example05, 64 k)
    vircef/   wannier90.{1,2}.{win,eig} + *_ref.wout        (纠缠金标;mmn/amn 不入库)
    li2mno3/  wannier90.{1,2}.{win,eig} + *_ref.wout        (纠缠金标,生产尺度)
    fetch_entangled.sh     # 从 sai-pod 重取两套纠缠算例的 .mmn/.amn(10~300 MB)
```
