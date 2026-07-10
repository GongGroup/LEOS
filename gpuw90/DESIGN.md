# gpuw90 — GPU Wannier-化后处理工具 · 设计文档 (Phase 1)

独立于 pyleos 的后处理工具。读 wannier90 标准格式(`.mmn`/`.amn`/`.win`,可选 `.nnkp`),
用 PyTorch 批量实现 Marzari–Vanderbilt(MV)幺正最小化(wannierise),输出 Wannier 中心与
spread 供 LEOS 氧化态管线消费。**不 fork wannier90、不改 VASP**;device 参数化,CPU/CUDA 同码。

本阶段(Phase 1)只做 **wannierise**(孤立带组,`num_bands == num_wann`,如绝缘体价带)。
disentangle(解纠缠)在本文档中只预留接口,留待 Phase 2。

---

## 1. 范围与阶段边界

| | Phase 1(本阶段,已完成) | Phase 2(预留) |
|---|---|---|
| 处理对象 | 孤立带组 `nb == nw` | 纠缠带组 `nb > nw`(金属/含导带) |
| 核心运算 | 每 k:1 个 `nw×nw` 厄米特征分解 + ~20 个 `nw³` matmul | 追加每 k 一个 `≤nb` 维厄米特征分解(子空间选取) |
| 输入 | `.mmn`(`nb×nb` 重叠)+ `.amn`(`nb×nw` 初猜)+ `.win` | 追加 `.eig`(能带能量,定冻结/外窗) |
| 优化器 | 最速下降(SD)+ 固定步长 | 追加共轭梯度(CG)+ 线搜索 |

设计原则:Phase 2 的 disentangle 作为 wannierise 的**前处理**,产出降维后的 `M0_reduced [nnb,nw,nw]`
与初始 `U0 [nw,nw]`,直接喂给本阶段已验证的 `wannierise()`,后者代码不改(见 §7)。

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
   │  read_win / read_nnkp     │                          │
   │  compute_bvectors         │  b 向量 (Å⁻¹)             │
   │  compute_shell_weights    │  w_b (B1 完备关系)         │
   │  → W90Case                │                          │
   └─────────────┬─────────────┘                          │
                 │  M0, A, bvec, wb, kb_idx                │
   ┌─────────────▼─────────────┐                          │
   │ wannierise.py (torch)     │                          │
   │  lowdin_orthonormalize    │  U0                       │
   │  wannierise():            │                          │
   │   循环: rotate M → 梯度 g   │                          │
   │        → exp(dW) → 更新 U   │                          │
   │  → centres[nw,3], spreads │                          │
   └─────────────┬─────────────┘                          │
                 │                          ┌─────────────▼──────────────┐
                 │                          │ io_w90.parse_wout_final_… │
                 │                          └─────────────┬──────────────┘
   ┌─────────────▼──────────────────────────────────────▼─────────────┐
   │ check_against_w90.py   最近邻(允许简并交换)+ 最小镜像距离对照        │
   │   逐 WF 偏差表 · 最大中心偏差 · PASS/FAIL(阈值 1e-4 Å)              │
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

## 7. Phase 2 (disentangle) 接口预留

纠缠情形 `nb > nw`,在 wannierise 之前插入子空间选取,产出降维输入喂给本阶段已验证的 `wannierise()`:

```python
def disentangle(M0_full, A_full, eig, bvec, wb, kb_idx,
                dis_win_min, dis_win_max, froz_win_min, froz_win_max,
                num_iter, device) -> (M0_reduced, U0):
    """
    M0_full : [nk, nnb, nb, nb]   全带重叠(直接来自 .mmn)
    A_full  : [nk, nb, nw]        全带初猜(直接来自 .amn)
    eig     : [nk, nb]            能带能量(来自 .eig),定外窗/冻结窗
    返回:
      M0_reduced : [nk, nnb, nw, nw]  V_k† M0_full V_{k+b}  → wannierise() 的 M0
      U0         : [nk, nw, nw]        由 V_k† A_k 经 Löwdin  → wannierise() 的 U0
    """
```

内部(对应 wannier90 `disentangle.F90` 的 `dis_extract`):每迭代每 k 组装
`Z_k = Σ_b w_b M^{k,b} P_{k+b} M^{k,b†}`(`nb×nb` 厄米),取前 `nw`(去掉冻结态占位)特征向量更新 `V_k`;
迭代最小化 Ω_I,收敛后固定 `V_k`。**关键:`wannierise()` 签名不变**——它已按 `M0[nnb,nw,nw]`/`U0[nw,nw]`
设计,disentangle 只是把 `nb` 压到 `nw`。需要新增 `io_w90.read_eig`(格式:`m k E` 三列)。

数值风险(可行性报告 §4/§8 已标注):冻结窗边缘近简并子空间选取是离散决策,对浮点顺序敏感 →
下游中心可移动。兜底:同一套金标回归 + 近平局告警 A/B 对比。

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

## 9. Phase 2+ 待办清单

1. **disentangle 阶段**(§7):`Z_k` 特征分解 + 冻结/外窗 + `read_eig`;金标用 lead/silicon(纠缠算例)。
2. **CG 优化器** + 线搜索(Fletcher–Reeves / Polak–Ribière),替换固定步长 SD,对齐 wannier90 在硬体系上的迭代数
   (SD 在近简并/大体系收敛慢;GaAs/Diamond 这类良态体系 SD 已 10~12 迭代收敛)。
3. **CUDA 实测**:sai-pod V100 上跑通(数学同码),量加速比;补 batched-eigh 在 n≈250 的现实性能核实
   (可行性报告 §4 的短板)。
4. **EOSCAR/pyleos 输出适配器** + 一个 Ce 帧端到端集成测试。
5. **更大金标集**:Ce / Li2MnO3 真实帧回归;大 `nw` 时把贪心中心匹配换 Hungarian。
6. **鲁棒性**:`.win` 无显式 kpoints 时的 MP 网格生成(需核实与 `.mmn` 一致的 k 排序)或强制走 `.nnkp`;
   `guiding_centres` 相位选项(棘手初猜时)。
7. **性能**:若剖析显示值得,改增量式 M 更新(wannier90 风格)替代每迭代全量重算旋转。

---

## 附:文件清单

```
gpuw90/
  io_w90.py             # .mmn/.amn/.win/.nnkp 解析 + b 向量/权重重建 + .wout 参考解析
  wannierise.py         # torch 批量 MV 最小化(device 参数化)
  check_against_w90.py  # 金标对照,逐 WF 偏差表 + PASS/FAIL
  fd_check.py           # 梯度符号/前因子的有限差分自检
  DESIGN.md             # 本文档
  README.md             # 运行说明(本地 CPU torch 复现)
  testdata/
    gaas/     gaas.{win,mmn,amn} + gaas_ref.wout        (example01, 8 k)
    diamond/  diamond.{win,mmn,amn} + diamond_ref.wout  (example05, 64 k)
```
