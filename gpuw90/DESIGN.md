# gpuw90 设计

gpuw90 是独立于 pyleos 的后处理工具。它读 wannier90 标准格式的输入，用 PyTorch 批量实现
Souza–Marzari–Vanderbilt（SMV）解纠缠和 Marzari–Vanderbilt（MV）局域化，输出 Wannier 中心和 spread，
供 LEOS 判定氧化态。安装、运行命令、测试数据和性能数字见 README.md。

## 范围

gpuw90 只实现 wannier90 中的解纠缠和局域化两步。它不 fork wannier90，不改 VASP；
计算设备由 `device` 参数决定，CPU 与 CUDA 共用一份代码。

| | 孤立带 | 纠缠带 |
|---|---|---|
| 条件 | `nb == nw` | `nb > nw`（金属或含导带） |
| 入口 | `wannierise()` | `disentangle()` → `wannierise()` |
| 输入文件 | `.mmn`（nb×nb 重叠）、`.amn`（nb×nw 初猜）、`.win` | 另加 `.eig`（能带能量，用于确定外窗和冻结窗） |
| 每次迭代每个 k 的主要运算 | 1 个 nw×nw 厄米特征分解，约 20 个 nw³ 矩阵乘 | 解纠缠阶段另有 1 个 nb 维厄米特征分解（Z 矩阵的子空间选取） |
| 迭代方法 | 固定步长最速下降（SD） | 解纠缠用 Z 矩阵不动点迭代，mixing 方式与 wannier90 相同；之后同孤立带 |

不在范围内的内容：

- SCDM 初猜。初猜直接读 VASP 写出的 `.amn`。
- 能带插值。因此省略了 wannier90 在解纠缠收敛后的一步规范旋转，见"算法"中的"与 wannier90 的差异"。

解纠缠是局域化的前处理。它输出降维后的 `M0_reduced [nk,nnb,nw,nw]` 和初始规范 `U0 [nk,nw,nw]`，
直接交给 `wannierise()`。`wannierise()` 不感知 nb，孤立带和纠缠带共用同一套局域化代码。

## 模块与数据流

数据从 VASP 写出的文件开始，经 `io_w90` 读入，按带组类型分两路，最后都进入 `wannierise()`。

```
VASP（LWRITE_MMN_AMN）
 │  写出 seedname.mmn、.amn、.eig；另备 .win
 ▼
io_w90.load_case()
 │  解析文件，重建 b 向量和 shell 权重，返回 W90Case
 │
 ├─ 孤立带（nb == nw）
 │    M0 = case.M
 │    U0 = wannierise.lowdin_orthonormalize(A)
 │
 ├─ 纠缠带（nb > nw；load_case 需传 allow_entangled=True 和 eig_path）
 │    disentangle.disentangle_case()
 │      窗口与冻结掩码 → 初始子空间 → Z 矩阵不动点迭代
 │      得到 M0_reduced 和 U0，作为下一步的 M0、U0
 │
 ▼
wannierise.wannierise(M0, U0, bvec, wb, kb_idx)
 │  每次迭代：旋转 M → 梯度 g → exp(dW) → 更新 U
 │  返回 centres [nw,3]、spreads [nw]、omega
 ▼
check_against_w90.py（孤立带）/ check_dis_against_w90.py（纠缠带）
    读参考 <seed>_ref.wout（io_w90.parse_wout_final_state、parse_wout_dis_info）
    按最近中心配对，用最小镜像距离比较，判定 PASS/FAIL（阈值 1e-4 Å）
```

### 批量布局

k 是批量维。所有重计算都是对 k（以及近邻 b）批量执行的 `torch.matmul` 和 `torch.linalg.eigh`，
对应 wannier90 中逐 k 调用的 `zgemm` 和 `zheev`。切换到 `device="cuda"` 时数学不变。

| 张量 | 形状 | 类型 | 含义 |
|---|---|---|---|
| `M0` | `[nk, nnb, nw, nw]` | complex128 | 每个 (k, b) 的第一性原理重叠矩阵 |
| `U` | `[nk, nw, nw]` | complex128 | 每个 k 的幺正规范 |
| `bvec` | `[nk, nnb, 3]` | float64 | 笛卡尔 b 向量（Å⁻¹） |
| `wb` | `[nk, nnb]` | float64 | shell 权重 |
| `kb_idx` | `[nk, nnb]` | long | 近邻 k 的索引，即 k+b 折回布里渊区后对应的 k |

解纠缠阶段的输入是全带维度：`M0_full [nk,nnb,nb,nb]`、`A_full [nk,nb,nw]`、`eig [nk,nb]`。

## 算法

### 局域化（wannierise.py）

局域化与 wannier90 3.1 的 MV 最速下降相同。每次迭代对全部 k 批量执行以下步骤：

1. 旋转重叠：`M^{k,b} = U_k† M0^{k,b} U_{k+b}`（`_rotate_M`，近邻的规范用 `U[kb_idx]` 取）。
2. 相位：`Im ln M_nn = atan2(Im, Re)`，取主支。
3. 中心：`r_n = −(1/N) Σ_{k,b} w_b b · Im ln M_nn`。
4. 梯度，即真实的 dΩ/dW，反厄米：
   `g_k = −(1/N) · 4 Σ_b w_b (A[R^{k,b}] − S[T^{k,b}])`。
   其中 `R_mn = M_mn·conj(M_nn)`，`T_mn = (M_mn/M_nn)·q_n`，`q_n = Im ln M_nn + b·r_n`，
   `A[B] = (B − B†)/2`，`S[B] = (B + B†)/(2i)`。
5. 步长：`dW_k = −(α·N / (4 w_tot)) g_k`，`w_tot = Σ_b w_b`。α = 1 是经典 MV 步长。
6. 更新：`U_k ← U_k · exp(dW_k)`。矩阵指数通过 `i·dW` 的厄米特征分解计算：
   `exp(dW) = V·diag(e^{−iλ})·V†`，与 wannier90 用 `zheev` 求指数的做法相同。
7. 收敛：`|ΔΩ| < conv_tol` 连续 `conv_window` 次，语义与 wannier90 相同。

孤立带的初始规范是 Löwdin 对称正交化 `U0 = A(A†A)^{−1/2}`，对应 wannier90 输出的 "Initial State"。

### 收敛判断与设备同步

热循环不在每次迭代都把数据读回主机，以免 GPU 上每次迭代都要等一次同步。

- 收敛计数 `n_below` 和 Ω 历史都保存在设备端张量中，循环结束后才一次取回。
- 每 `check_every` 次迭代只读回一个标量收敛标志。`check_every > 1` 时，循环最多比实际收敛多跑 `check_every − 1` 次迭代。
- 更新 U 后旋转得到的 M 直接作为下一次迭代的输入，每次迭代只旋转一次。
- 解纠缠循环用同样的做法。

### 解纠缠（disentangle.py）

解纠缠在每个 k 上选出使规范不变 spread Ω_I 最小的 nw 维子空间。底层函数签名：

```python
res = disentangle(M0_full, A_full, eig, bvec, wb, kb_idx,
                  dis_win_min=None, dis_win_max=None,
                  froz_win_min=None, froz_win_max=None,
                  num_iter=200, mix_ratio=0.5,
                  conv_tol=1e-10, conv_window=3,
                  device="cpu", dtype=torch.complex128,
                  check_every=1, verbose=False)
# res.M0_reduced、res.U0 直接交给 wannierise()
```

便捷入口 `disentangle_case(case)` 从 `.win` 读取窗口和迭代参数，`.win` 中没有的项用 wannier90 默认值：

| `.win` 关键字 | `disentangle()` 参数 | 缺省值 |
|---|---|---|
| `dis_win_min`、`dis_win_max` | `dis_win_min`、`dis_win_max` | 全部本征值的最小值、最大值 |
| `dis_froz_max` | `froz_win_max` | 不设冻结窗 |
| `dis_froz_min` | `froz_win_min` | 只给了 `dis_froz_max` 时取 `dis_win_min` |
| `dis_num_iter` | `num_iter` | 200 |
| `dis_mix_ratio` | `mix_ratio` | 0.5 |
| `dis_conv_tol` | `conv_tol` | 1e-10 |
| `dis_conv_window` | `conv_window` | 3 |

实现逐步对照 wannier90 3.1 `disentangle.F90` 的源码，而不只是文献中的公式：

1. 窗口（对应 `dis_windows`）：外窗在每个 k 上选出 `ndimwin(k)` 条带。`.eig` 读入时检查能量升序，
   所以选出的带是连续的。冻结窗内的 `ndimfroz(k)` 条带必须整条保留。
   `ndimwin < nw` 或 `ndimfroz > nw` 时报错，与 wannier90 相同。
2. 初始子空间（对应 `dis_project` 和 `dis_proj_froz`）：`V0 = polar(A_w)`，即批量 SVD 的幺正因子，
   与 Löwdin 正交化等价。有冻结窗时，非冻结部分取 `Q_f P_s Q_f` 的前 `nw − nf` 个本征向量（SMV Eq. 27），
   冻结方向直接用单位向量。
3. 主循环（对应 `dis_extract`）：Jacobi 式迭代，所有 k 同时更新，与 wannier90 按 MPI 分配 k 点时的语义相同。
   - `Z_in⁽¹⁾ = Z[V⁽⁰⁾]`；`Z_in⁽ⁱ⁾ = mix·Z_out⁽ⁱ⁻¹⁾ + (1 − mix)·Z_in⁽ⁱ⁻¹⁾`
   - `V⁽ⁱ⁾` = 冻结态，加上 `Z_in⁽ⁱ⁾` 的前 `nw − nf` 个本征向量
   - `Z_out⁽ⁱ⁾ = Z[V⁽ⁱ⁾]`
   - 其中 `Z[V](k) = Σ_b w_b M^{k,b} P_{k+b} M^{k,b†}`，`P = V V†`，再投影到 k 点的非冻结补空间。

   mixing 混合的是 Z 矩阵，不是投影算符。`Z_out` 在子空间更新之后、迭代末尾组装。这两点都与 wannier90 相同。
4. 收敛判据：每次迭代计算两个 Ω_I。
   - `Ω_I(i−1)`（SMV Eq. 18）的近邻用上一次迭代的子空间，等于冻结贡献加上选中本征值之和。
   - `Ω_I(i)`（SMV Eq. 12）全部用新子空间。
   - `delta = Ω_I(i−1)/Ω_I(i) − 1`，连续 `conv_window` 次 `|delta| < conv_tol` 即收敛。
5. 初始规范（对应 `internal_find_u`）：`U0 = polar(V† A_w)`，即 wannier90 的 Löwdin 投影初值。

### 逐 k 不等的维度并入一次 batched eigh

`ndimwin` 和 `ndimfroz` 随 k 变化，例如 li2mno3 自旋 1 的 nf 在 284~286 之间。逐 k 循环会破坏批量，
所以实现保留全 nb 维，用掩码和对角能级移位代替切片：

- Z 先把冻结行列和窗外行列置零，再在对角线上加移位：冻结方向 `+2·w_tot`，窗外方向 `−2·w_tot`。
  因为 ‖M‖ ≤ 1 且 w_b > 0，Z 的谱落在 `[0, w_tot]` 内，所以升序本征值排成 `[窗外 | 自由谱 | 冻结]`。
- 一次 `torch.linalg.eigh` 的最后 nw 列就是新子空间。冻结块是精确的 `2w_tot·I`，
  它的本征向量是冻结单位向量的幺正组合。子空间是精确的，列内的规范没有物理效应。
- `Ω_I(i−1)` 需要的选中本征值之和，用逐 k 的槽位掩码 `nb − nw ≤ pos < nb − nf(k)` 向量化计算。
- 移位同时加在 `Z_in` 和 `Z_out` 中。mixing 是线性的，移位在混合后保持不变，槽位选择在每次迭代都有效。
- eigh 之后按列类型用掩码清除约 1e−15 的跨块浮点噪声（`_scrub_columns`），冻结态保持严格冻结。
- 初始子空间的 `QPQ`（谱在 `[0, 1]`）用同样方法处理，移位取 ±2。补空间块与冻结块精确解耦，
  所以近零本征值对应的本征向量自动与冻结态正交，wannier90 中 ortho-fix 要保证的也是这一性质。

### 与 wannier90 的差异

wannier90 在解纠缠收敛后，会把 V 右乘子空间哈密顿量的本征向量，供能带插值使用。gpuw90 不做能带插值，省略了这一步。

这一步只改变规范。右乘 `V → VR` 在 `M0_reduced = V† M0 V'` 与 `U0 = polar(V† A)` 之间严格相消，
因为 `polar(R†B) = R†·polar(B)`。Ω_I、中心、spread 和局域化轨迹都不受影响。
纠缠带初态与参考在 `.wout` 打印精度内一致（见"正确性验证结果"），与这一推导相符。

## 数值问题与处理

下面各项在代码中都有就地注释。

### 梯度的 1/N 前因子和整体符号

文献中 MV 梯度常写成 `4 Σ_b w_b (A[R] − S[T])`。但 Ω 本身带 `(1/N) Σ_k`，实现中也容易差一个整体符号。

gpuw90 用有限差分方向导数检验梯度（`fd_check.py`）：比较解析的 `Re Tr[g† X]` 与
`[Ω(U e^{εX}) − Ω(U e^{−εX})] / 2ε`。

| 梯度写法 | GaAs 上的 ratio |
|---|---|
| `4 Σ_b w_b (A[R] − S[T])` | −1/8，即 −1/N（N = 8） |
| `−(1/N) · 4 Σ_b w_b (A[R] − S[T])`（当前实现） | +1.0000，三个 ε 下相同 |

因此真实梯度是 `g = −(1/N) · 4 Σ_b w_b (A[R] − S[T])`，按 `dW = −eff·g` 下降时 Ω 单调减小。
改动梯度后先跑 `fd_check.py`，不要直接采用公式中的前因子和符号。

### Ω_I 与 Ω_OD 的拆分

Ω_I 必须对全部 m、n 求和，Ω_OD 只取非对角部分：

- `Ω_I = (1/N) Σ_{k,b} w_b (nw − Σ_{m,n} |M_mn|²)`
- `Ω_OD = (1/N) Σ_{k,b} w_b Σ_{m≠n} |M_mn|²`
- `Ω_D = (1/N) Σ_{k,b} w_b Σ_n (Im ln M_nn + b·r_n)²`

常见错误是把 `Σ_n |M_nn|²` 放进 Ω_I。这样 Ω_OD 先被并进 Ω_I，又被单独加一次，总 Ω 多出一个 Ω_OD。
症状是中心正确（偏差 6e-7 Å）而总 spread 偏高 0.502，等于 Ω_OD。
按上式实现时，GaAs 初态 Ω 与参考逐位一致（4.46881212）。

### ln 的支割线

`Im ln M_nn = atan2(Im, Re)`，取值在 (−π, π]，与 wannier90 的 `aimag(log(...))` 相同。
CPU 和 GPU 取同一主值，不会因设备不同而跳支。

绝缘体 WF 的 M_nn ≈ 1，相位小，远离 ±π。只有在接近平局或 |M_nn| → 0 时才可能翻支，
使中心跳动 `b·2π/…` 量级。孤立带验收算例的 Ω_D 很小，没有触发。

### b 向量的 shell 权重

权重按 wannier90 的 kmesh 规则重建：

1. b 向量按 |b| 分 shell。
2. 用最小二乘解 `Σ_shell w_s (Σ_{b∈s} b_α b_β) = δ_αβ` 的 6 个独立对称分量。
3. 均匀网格上各 k 的 shell 结构相同，权重按 |b| 匹配回填到每个 (k, b)。

加载时计算 B1 残差 `max|Σ_b w_b b⊗b − I|`：GaAs 为 9e-18，Diamond 为 2e-16。
w_b 与参考逐位一致（GaAs 0.408635，Diamond 0.527876）。

### 简并本征向量的相位（不需要处理）

`i·dW` 的本征向量在简并子空间内只确定到一个幺正旋转或相位。
但 `exp(dW) = V·f(λ)·V†` 是算符函数，与这个规范无关，所以简并不影响更新。
Löwdin 的 `(A†A)^{−1/2}` 同理。代码中注明了这一点，不需要加相位固定逻辑。

### 全程 FP64

计算全部使用 complex128 / float64，与 wannier90 的双精度一致。

CPU 与 GPU 只在归约顺序上不同，差异约 1e-14。收敛到同一极小时，这对中心的影响远小于 1e-4 Å；
验收实测的偏差在 1e-6~1e-9 Å 量级（见"正确性验证结果"）。
gpuw90 不要求 CPU 与 GPU 逐位一致，只要求物理上等价。

### WF 顺序与周期镜像

最小化得到的 WF 顺序可能与参考不同，近简并的 WF 也可能互换。

验收脚本因此按最近中心贪心配对，不按序号；距离用最小镜像，即先在实空间晶格中折叠再取长度。
这样可以容忍简并交换和晶格平移。阈值为 1e-4 Å。

### .mmn 的列主序

每个 (k, b) 块内的 nb×nb 数据以带序 m 为最快变化的索引，即 Fortran 列主序。
读入时用 `reshape(order="F")` 得到 `M[m,n]`。

如果误按行主序读，中心的对称性会被破坏，但 B1 残差仍然正常。
初态中心与参考的对照能直接发现这个错误。

### 解纠缠：两个 Ω_I 的计算次序

次序必须如下：

1. `Ω_I(i−1)` 的冻结贡献用上一次迭代的近邻子空间，在更新任何 k 之前整体算出；
2. `Ω_I(i)` 用新子空间；
3. `Z_out` 在子空间更新之后、迭代末尾组装。

任何一处次序不对，delta 序列就会偏离 wannier90，收敛迭代数对不上。
按此次序实现时，四个纠缠带通道的迭代数与 wannier90 逐一相同。

### 解纠缠：能级移位的幅度

冻结方向和窗外方向的对角移位如果取很大的数（如 1e6），eigh 本征向量的误差（约 eps·‖A‖/gap）
会被放大并混进自由子空间。

gpuw90 取与谱同量级的 `±2·w_tot`。因为 ‖M‖ ≤ 1，Z 的谱有界于 `[0, w_tot]`，
这个移位已足够把各块分开，又不损失精度。初始子空间的 `QPQ`（谱在 `[0, 1]`）同理取 ±2。

## 正确性验证结果

孤立带 2 个算例和纠缠带 4 个通道全部 PASS，最大中心偏差都在 1e-6 Å 以下，验收阈值为 1e-4 Å。
运行命令和参数见 README.md 的"运行验收"。

### 孤立带

参考输出是 wannier90 3.1 test-suite 的 `benchmark.out`，完全收敛：GaAs 从第 2 次迭代起中心不再变化。

| 算例 | nb/nw | nk | nntot | 收敛迭代数 | 最大中心偏差（Å） | 最大 spread 偏差（Ų） | Ω_total 偏差 |
|---|---|---|---|---|---|---|---|
| GaAs（example01） | 4/4 | 8 | 8 | 10 | 8.31e-07 | 4.07e-09 | 2.6e-10 |
| Diamond（example05） | 4/4 | 64 | 8 | 12 | 3.13e-09 | 1.98e-09 | 1.3e-10 |

- 阈值 1e-4 Å，两例的余量为 10⁵~10⁷ 倍。Diamond 有 64 个 k 点，覆盖了批量 k 路径。
- 初态（Löwdin 后、未优化）的中心与参考 "Initial State" 相差 6e-7 Å（GaAs）和 3e-9 Å（Diamond）。
  这一项不经过最小化，单独验证了文件读取、shell 权重和中心公式。

复现：`python3 check_against_w90.py testdata/gaas --seed gaas`；Diamond 换成 `testdata/diamond --seed diamond`。

### 纠缠带

参考输出来自 SAI 集群上 wannier90.x 3.1 的串行运行（`run_np1_s*/wannier90.*.wout`）。
参考运行与 gpuw90 读的是同一套 VASP 导出文件（`.mmn/.amn/.eig/.win` 通过符号链接指向同一份）。
两套算例都有冻结窗（`dis_froz_max` 取自 `.win`），外窗为默认的全范围。
li2mno3 与 LEOS 生产计算的体系规模相同（288 个 WF）。

| 算例（自旋） | nb/nw | nk/nnb | nf | 解纠缠迭代数（gpuw90 / 参考） | Ω_I（Ų，两边打印值相同） | \|ΔΩ_I\| | 初态最大中心偏差（Å） | 初态最大 spread 偏差（Ų） |
|---|---|---|---|---|---|---|---|---|
| vircef s1 | 180/142 | 8/12 | 141 | 82 / 82，收敛 | 81.06033903 | 3.4e-09 | 7.9e-07 | 5.0e-09 |
| vircef s2 | 180/142 | 8/12 | 140 | 188 / 188，收敛 | 80.32720233 | 4.6e-10 | 7.7e-07 | 4.9e-09 |
| li2mno3 s1 | 356/288 | 8/8 | 284~286 | 38 / 38，收敛 | 154.86511542 | 4.4e-09 | 8.3e-07 | 5.0e-09 |
| li2mno3 s2 | 356/288 | 8/8 | 240 | 200 / 200，两边都到 `dis_num_iter` 上限，未收敛 | 135.30457827 | 4.1e-09 | 8.2e-07 | 5.0e-09 |

验收门限：|ΔΩ_I| ≤ 1e-4 Ų，初态最大中心偏差 ≤ 1e-4 Å。四个通道全部 PASS。

验收量的选择依据：

- 参考 `.wout` 的局域化段用 CG，固定 100 次迭代，没有收敛判据（`conv_window = −1`）。
  四个参考末态都没有收敛，最后几次迭代的 ΔΩ 仍有 1e-2~1e-1 Ų/迭代。
- 因此末态中心的差异反映的是两边运行设置不同。实测这一差异约 0.3 Å，
  而且 gpuw90 收敛后的 Ω_total 更低（vircef s1：94.2618，参考为 94.6541）。它不能用来判断实现对错。
- 所以硬验收只用与运行设置无关的量：Ω_I（解纠缠的规范不变输出），以及解纠缠后的 "Initial State" 中心和 spread。
  后者是子空间加 Löwdin 投影的初值，确定性计算，不经过任何最小化。
- 选验收量之前，应先看参考 `.wout` 中回显的运行参数，确认参考本身是否收敛。

结果解读：

- 初态中心偏差约 8e-7 Å、spread 偏差约 5e-9 Ų，对应 `.wout` 打印精度的舍入（中心 6 位小数，spread 8 位小数）。
  孤立带 GaAs 的 8.3e-7 Å 也是同一来源。实际吻合程度比这更高。
- 四个通道的解纠缠迭代数与 wannier90 逐一相同（38、82、188、200）。
  这说明 mixing、两个 Ω_I 的计算次序和收敛窗口语义都与 wannier90 一致。
- li2mno3 s2 在 200 次迭代截断、未收敛的轨迹上仍复现到 4e-9 Ų，说明不动点迭代中的浮点漂移可以忽略。
- 脚本每次都会检查 `|Ω_I(初态) − Ω_I(解纠缠)|`。这是规范不变量，实测约 1e-11。

CUDA 上的验收：V100 上 li2mno3 同样 PASS（|ΔΩ_I| = 4.4e-9 Ų，初态中心偏差 8.3e-7 Å），
测量条件见 README.md 的"性能实测"。

复现：`python3 check_dis_against_w90.py testdata/vircef --seed wannier90.1`，其余通道换目录和 `--seed`。
大输入不入库，先运行 `testdata/fetch_entangled.sh`。

## 显存与分块

目前的工作集很小，显存不构成约束。代码不做分块，全部张量常驻设备。

- GPU 可行性评估估算，解纠缠和局域化的全部工作集 < 2 GB。出处：
  `docs/reports/20260710_gpu_feasibility.md` 的"规模画像与基线拆解"一节；该文件已移出工作区，在 git 历史中可查。
- 以 Ce 体系（nk = 32，nw ≈ 253，nntot ≈ 8~12）估算：`M0 = 32 × 10 × 253² × 16 B ≈ 0.33 GB`，`U ≈ 0.1 GB`。
  V100-SXM2-16GB 单卡足够。
- 解纠缠阶段在整个迭代期间保留全带的窗口掩码重叠 `M0w [nk,nnb,nb,nb]`。

体系更大、需要分块时，可以按两级切分：

| 级别 | 切法 | 依据 |
|---|---|---|
| 自旋 | ISPIN = 2 的两个通道依次处理，或分到两块 GPU 上并行 | 两个自旋通道完全独立。gpuw90 本身按 seed 一次处理一个自旋通道 |
| k 点 | `M0` 按 k 切块从主机流入；U 全量常驻设备，只换入 `M0[k_tile]` 及其近邻块 | 旋转 M 要用近邻的 `U[k+b]`，而 U 很小。梯度对 b 的求和、中心的归约都只涉及本 k，跨块的只有近邻 U 的查表 |

k 分块只有超大体系才需要。

显存压力最大的是 SCDM 中 `nb × Nr ≈ 5.3 GB` 的矩阵，但 SCDM 不在 gpuw90 的范围内：初猜直接读 VASP 写出的 `.amn`。

## 与 LEOS 管线的对接

gpuw90 可以在 VASP 之外离线完成 Wannier 化中解纠缠和局域化两步，代替 VASP 内 wannier90 库模式的这两步。

| 环节 | 现状 |
|---|---|
| 输入 | VASP 设置 `LWRITE_MMN_AMN`（mlwf.F 中已有的开关），写出 `seedname.mmn`、`.amn`、`.eig`；另需一个提供晶胞和显式 kpoints 的 `.win`，可由 VASP 写出或手写。gpuw90 离线读取这些文件，不进入 VASP 进程，显存由 gpuw90 自行管理 |
| 输出 | 每个自旋通道的 `centres [nw,3]`（Å）和 `spreads [nw]`，以 numpy 数组的形式放在 `WannieriseResult` 中返回 |
| 帧间热启动 | `wannierise()` 接受任意幺正 `U0`，可以传入上一帧收敛的 `res.U` 作为初值。仓库中没有针对这一用法的测试 |

目前没有把输出写成 EOSCAR 或交给 pyleos 的适配器，见"已知限制"。

## 已知限制

- 局域化只有固定步长最速下降，没有 CG 和线搜索。良态体系（GaAs、Diamond）10~12 次迭代收敛；
  在 vircef、li2mno3 这类近简并体系上，SD 跑 2000 次迭代仍未达到 1e-10。集群 wannier90.x 的 CG 100 次迭代也没有收敛。
- 冻结窗边缘近简并时，子空间选取是离散决策，对浮点运算顺序敏感，下游中心可能因此移动。
  目前只有 wannier90 参考输出回归（`check_dis_against_w90.py`）一道防护。
- 验收脚本按最近中心贪心配对 WF（`check_against_w90.match_wfs`），不是全局最优配对。nw 很大时可能配错。
- `.win` 没有显式 kpoints 块时，`load_case()` 必须传 `nnkp_path`，不会根据 `mp_grid` 自动生成 k 点。
- 不支持 wannier90 的 `guiding_centres` 选项。
- 没有输出适配器：结果不会写成 EOSCAR，也没有供 pyleos 读取的文件格式。
- 局域化每次迭代都从 `M0` 整体重新旋转 M，没有 wannier90 式的增量更新。
- 纠缠带验收只有 vircef 和 li2mno3 两套算例，没有 Ce 体系的真实帧。
