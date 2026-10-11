# Comas-Vilà et al. (JCTC 2025) 阅读笔记

本文记录 Comas-Vilà 等人把 effective oxidation state（EOS）分析推广到固体的
方法、他们的可靠性指标 R(%)、这套方法与 LEOS 的异同，以及 LEOS 在其
40 固体测试集上的结果。LEOS 的数字以 docs/RESULTS.md 和 data/ 下的文件为准。

## 文献信息

| 项目 | 内容 |
|---|---|
| 作者 | Gerard Comas-Vilà, Leila Pujal, Alberto Otero-de-la-Roza, Davide Tiana, Julia Contreras-García, Pedro Salvador |
| 标题 | Effective Oxidation State Analysis for Solids |
| 出处 | *J. Chem. Theory Comput.* 2025, 21, 7075–7086 |
| DOI | 10.1021/acs.jctc.5c00482 |
| PDF | docs/refs/Comas-Vilà 等 - 2025 - Effective Oxidation State Analysis for Solids.pdf |
| 分子版前身 | Ramos-Cordoba, Postils, Salvador, *J. Chem. Theory Comput.* 2015, 11, 1501–1508，即原文参考文献 20 |

下文页码都是期刊页码（7075–7086）。

这篇文章与 LEOS 的关系：

- 两者回答同一个问题：从 DFT 波函数为固体中每个原子判定整数氧化态。
- 技术路线不同。他们用 QTAIM 原子域和原子重叠矩阵；LEOS 用 Wannier 中心
  的位置和带容量约束的运输问题（最小费用流）。
- "EOS"这个名字来自分子版方法。
- 他们的 40 固体测试集是 LEOS 扩展基准的来源。LEOS 对其中全部 40 个都有
  判定，见下文"40 固体对照结果"。

## 他们的方法

他们先在每个原子域内构造有效原子轨道并求出占据数，再把全体系的这些轨道
按占据数排序，一次性分配整数个电子。

### 思想

Theory 节（p.7076–7077）按单行列式波函数和 QTAIM 划分描述了以下步骤，
固体版沿用同一算法（p.7078）：

1. 划分空间。用 QTAIM（Bader 电子密度拓扑划分）把空间切成互不重叠的
   原子域 Ω_A，每个原子一个吸引盆。
2. 计算原子重叠矩阵（atomic overlap matrix, AOM）。对每个原子域，计算
   占据轨道在域内的重叠积分：
   S_ij^(σ,A) = ∫_{Ω_A} φ_i^σ*(r) φ_j^σ(r) dr（式 1，p.7077）。
3. 对角化 AOM。本征向量是只在 Ω_A 内定义的有效原子轨道（eff-AO）；
   本征值 λ 在 0 到 1 之间，是该轨道的占据数（式 2–4，p.7077）。
   占据数接近 1 的对应芯层或孤对，占据数中等的对应参与成键的轨道。
   以片段为单位时，同样的对象叫有效片段轨道（EFO）。
4. 全局分配。把所有原子或片段的 eff-AO 按占据数从大到小排成一列，
   前 n^σ 个记为占据，其余记为空；n^σ 是自旋通道 σ 的电子数。
   原子的氧化态等于原子序数减去分到的 α、β 电子数。自旋极化体系的
   两个通道分别做。

第 4 步不把占据数逐个四舍五入，而是让全体系的占据数互相比较后截断
（p.7077；Rb2CuCl4 与 XeO2 的讨论再次强调，p.7081、p.7082）。所以缺电子
体系中最后一个记为占据的 eff-AO，占据数可以低于 0.5；富电子体系中也会有
占据数高于 0.5 却记为空的 eff-AO（p.7077）。

LEOS 在这一点上与他们相同：整数氧化态来自一次全局分配，不对连续量取整。

### 程序链

平面波赝势 DFT 不直接给出 QTAIM 需要的全电子密度，所以他们串联了四个程序
（Results 节开头，p.7078；Computational Details，p.7084）。

| 步骤 | 程序 | 输入 | 输出 | 作用 |
|---|---|---|---|---|
| 1a. DFT 自洽 | Quantum ESPRESSO，norm-conserving 赝势，PBEsol | 晶体结构 | Bloch 态、电子密度 | 提供波函数 |
| 1b. 全电子密度重建 | Quantum ESPRESSO，同几何、同网格的另一次 PAW 计算 | 晶体结构 | 全电子密度格点 | QTAIM 拓扑分析需要全电子密度 |
| 2. Wannier 化 | wannier90 | 1a 的 Bloch 态 | MLWF | 把 Bloch 态变成局域轨道，代替分子中的分子轨道 |
| 3. 原子域积分 | critic2 | 全电子密度、MLWF | 每个原子的 AOM | 先做 QTAIM 拓扑分析找吸引盆，再用 Henkelman 格点法在每个盆内积分 MLWF 的乘积 |
| 4. EOS 分配 | APOST-3D | 全部 AOM | eff-AO 占据数、整数氧化态、R(%) | 对角化 AOM，按上一节第 4 步做全局排序 |

计算参数（p.7084）：PBEsol 泛函；波函数截断 80 Ry，密度截断 320 Ry；
k 点网格 4×4×4。

### 计算代价

主要开销是全电子密度重建和盆内格点积分。原文没有报告计时。

- 两遍 DFT。norm-conserving 计算提供 Bloch 态，PAW 计算提供全电子密度，
  两者的几何和网格必须一致（p.7084）。
- 盆内积分。每个原子域都要算一个 N_wann × N_wann 的 AOM，每个矩阵元是
  一次三维格点数值积分。320 Ry 的密度截断对应较密的实空间格点。
- 金属。部分填充的能带做不了 Wannier 变换，AOM 改由 Bloch 态直接计算
  （p.7084）。作者在 Rb3AuO 处指出，这样效率更低、计算量更大（p.7079）。
- 电子化物。游离电子要靠密度中的非核吸引子（non-nuclear attractor, NNA）
  作为额外的盆来记账。作者发现 Ca2N 中 NNA 的有无和个数对 k 点数和拓扑
  分析的格点很敏感（p.7079）。

### 输入与输出

- 输入：晶体结构。原文 40 个结构取自 COD、CCDC 或原始文献（Table 1，
  p.7079）。
- 输出：每个原子或片段的整数氧化态、各自旋通道的 eff-AO 占据数、全胞
  一个 R(%)。40 个体系的边界 eff-AO 占据数和判定列在 SI 的 Table S1–S16
  （p.7078）。
- 片段：QTAIM 原子域互不重叠，所以片段的 AOM 等于成员原子 AOM 之和
  （p.7078）。

## R(%) 可靠性指标

R(%) 用排序边界上两个 eff-AO 的占据数之差，衡量整数分配有多明确。

### 定义

设自旋通道 σ 排序后，最后一个记为占据的 eff-AO 占据数为 λ_LO^σ
（last occupied），第一个记为空的为 λ_FU^σ（first unoccupied）。
式 5（p.7077）：

    R^σ(%) = 100 · min(1, max(0, λ_LO^σ − λ_FU^σ + 1/2))

全胞取两个自旋通道中较小的值：R(%) = min(R^α, R^β)。

固体版对边界轨道的规定与分子版不同（p.7078）：

- 分子版要求两个边界 eff-AO 属于不同片段。
- 固体版允许它们属于同一个原子。理由是超胞论证：原文举了一个二元固体
  AB 的例子，按分子版规则原胞得 R = 90；换成两倍超胞后，边界轨道落在
  两个等价原子上，R 随之改变。放开这条限制后，原胞和超胞给出相同的 R。

### 取值含义

| R(%) | 含义 | 出处 |
|---|---|---|
| 100 | λ_LO − λ_FU ≥ 0.5，作者认为分配无可争议 | p.7077 |
| 50–100 | 值越高，形式离子图像越接近实际电子分布 | p.7077 |
| 50 | 两个边界 eff-AO 占据数相等，即简并 | p.7078 |
| 低于 50 | 分配没有按占据数从高到低进行。正常运行 EOS 不会出现，只出现在人为指定一套氧化态、检验体系与之符合程度的用法中 | p.7077–7078 |

简并由对称性造成时，作者把最后的电子平摊给简并的 eff-AO，得到分数
氧化态。通则在 p.7078；例子有 GaSe（p.7080）、MH4 氢化物
（p.7081–7082）、YBaFe2O5 高温相（p.7081）。

作者对 R(%) 的含义做了两点说明（p.7078）：

- R 不是概率。R = 65% 不表示判定有 65% 的可能正确；它衡量形式离子图像
  与实际电子分布的相似程度。
- 式 5 的函数形式有一定任意性。更重要的是找出边界 eff-AO，分析它们的
  形状和占据数。

### 原文中的用法

原文用 R(%) 做了四类事，每类都有实例。

| 用法 | 例子 | 出处 |
|---|---|---|
| 标注判定是否明确 | CrO3：Cr 判 +6，QTAIM 部分电荷只有 +1.76，R = 71。作者以此说明部分电荷不能用来推氧化态 | p.7078–7079 |
| 在不同片段划分之间选择 | GaSe：按原子划分 R = 50，把 Ga–Ga 作为片段后 R = 77.7，于是取 (Ga–Ga)⁴⁺ 图像。MH4：Mg/Ca/Sr 用 Ha–Ha 片段更明确；Zr/Ce/Th 按原子划分 R 更高 | p.7079–7080；p.7081–7082 |
| 看趋势 | HgF4：50、100、200 GPa 三个结构的 R 随压力单调升高，Hg(+4) 图像随之更可信 | p.7083 |
| 指出方法在何处不适用 | LiN5：五个 N 的占据数几乎相同，R = 50，给不出原子级判定。XeO3 150 GPa 相：R = 51.2，作者认为判定本质上模糊 | p.7084；p.7083 |

结论节把最后一类称为方法能 "self-assess its limits of applicability"
（p.7084）。

## 与 LEOS 的对比

两个方法都用一次全局分配得到整数氧化态，都不对部分电荷取整。差别在记账
单元、所需数据和计算链。

### 原理

| | Comas-Vilà | LEOS |
|---|---|---|
| 记账单元 | eff-AO，占据数是 0 到 1 之间的连续量 | Wannier 函数（WF），每个装整数个电子：闭壳层 2 个，自旋极化每通道 1 个 |
| 归属依据 | 全体 eff-AO 按占据数排序 | 运输问题（最小费用流），代价 = d² + β·(1 − χ/χ_max)。d 是 WF 中心到原子的最近周期像距离（Å），χ 是 Allen 电负性，β = 10 |
| 整数化 | 排序后取前 n^σ 个 | 供给和容量都是整数，最优解自动是整数 |
| 物理约束 | 没有显式容量约束 | 每个原子有容量上下限：上限是价层的物理容量，下限是半芯电子 |
| 所需数据 | 每个原子一个 AOM，N_wann × N_wann，矩阵元是盆内积分 | WF 中心坐标和 spread（EOSCAR）、晶体结构（POSCAR） |
| 空间划分 | QTAIM 原子域 | 不划分空间，只用中心到各原子的距离 |

他们把化学信息放在域内重叠积分里。LEOS 只用 WF 中心的位置和元素电负性，
再用容量约束排除物理上不可能的分配，例如 S 收 14 个电子。运输形式化的
说明见 docs/PLAN.md「把分配写成最小费用流」。

### 计算代价

| | Comas-Vilà | LEOS |
|---|---|---|
| 程序 | QE（两遍）→ wannier90 → critic2 → APOST-3D | 打了 LEOS 补丁的 VASP 一次运行，以库模式调用 wannier90，写出 EOSCAR；pyleos 可离线重算 |
| 主要开销 | 全电子密度重建、QTAIM 拓扑分析、盆内格点积分 | DFT 和 Wannier 化本身；分配问题只有几十到几百个变量，求解是毫秒级 |
| 金属 | 改用 Bloch 态计算 AOM，更慢 | 仍用 WF 中心；多出的 WF 作为 ghost 排除，是近似处理 |
| 计时 | 原文未报告 | pyleos 在笔记本上重算 74 案例约 25 秒，不含 R 指标 |

ghost 指多出来的、无电子可装的 WF。LEOS 用总运输量上限（真实电子数）把
它们排除。这个处理有已知偏差：封顶时先去掉代价最高的边，而代价含电负性
项，所以上游丢失的电子会集中记在电负性最低的原子上。ghost > 0 时，电正性
金属的氧化态可能被高估 1–2（docs/RESULTS.md「ghost 对金属价态的系统偏差」）。他们对金属有
Bloch 态这条路，计算量更大，但不需要这一近似。

他们的成本主要来自现行实现，不来自 EOS 框架本身：

- 分子版 EOS 对 Hilbert 空间划分和实空间划分都适用（p.7076–7077）。
- 结论节说该方案可推广到任何实空间划分，包括 Voronoi 多面体和
  Hirshfeld 型划分，前提是能算出相应的 AOM（p.7084）。
- 对电子化物，作者提出 ELF 与 Voronoi 结合的划分（Weaver et al.）可能比
  QTAIM 更合适（p.7079）。

换成更便宜的划分，每个原子仍要算一个 N_wann × N_wann 的域内积分矩阵；
LEOS 不做空间积分。docs/PLAN.md「默认代价只用几何和电负性」也因为成本，没有把完整的
原子重叠矩阵纳入 LEOS 的代价函数。

### 功能：轨迹

原文 40 个体系都是单点计算，没有分子动力学内容。按他们的程序链，逐帧分析
意味着每帧重做 PAW 重建和盆内积分。

LEOS 的 EOSCAR 带原位判定列 OS(flow)，分子动力学中可以逐帧判定。已有的
轨迹结果：

| 体系 | 规模 | 结果 | 出处 |
|---|---|---|---|
| Li15O16 | 40 步室温 MD，每 2 步分析一次，共 20 帧 | O 的元素平均氧化态 20 帧都是 −0.94 | docs/RESULTS.md「Li2O2 脱锂：氧的价态」 |
| 富锂正极 Li2MnO3、Li2RuO3、NMC | 3 个体系 × 3 个脱锂点 × 50 ps，共 180 帧 | 价态时间序列 | docs/RESULTS.md「富锂正极：价态时间序列」；data/m5b/ |

逐帧 Wannier 化的耗时，以 96 原子 Li2MnO3 脱锂帧计：完整 MLWF 的 wannier90
段 554 秒/帧；SCDM-only 快速模式 8.3 秒/帧，但在该体系上元素平均一致
0/12，失效帧都被 stray 警告拦下；gpuw90 按实测迭代速度折算，两个自旋的
Wannier 段约半分钟。数据和适用范围见 docs/RESULTS.md「效率」。

### 功能：片段

两边都支持把几个原子合成一个片段，用于"每对电子归到一两个原子"的图像
不成立的场合。

| | Comas-Vilà | LEOS |
|---|---|---|
| 实现 | 片段 AOM = 成员原子 AOM 之和（p.7078） | 成员原子合并成一个仓库：容量相加，代价取最近成员。CLI 参数 `--fragment` |
| 用例 | SrBe(OH)4 的 OH、Ba3Si4 的 Si4、GaSe 的 Ga–Ga、MH4 的 Ha–Ha、XeO3 的 O2、CsH3/CsF3/CsF5 的 H3/F3/F5 | hi3 的三碘单元、WAu12 的 Au12、Ce18O36Co4CH4 的 Co4 |
| 何时用片段 | 比较不同划分下的 R(%) | 某对电子的两个候选原子代价接近时，CLI 提示考虑 `--fragment` |
| 对称简并 | 把最后的电子平摊给等价 eff-AO，得到分数氧化态 | 在元素平均或片段总价的层级下结论 |

LEOS 侧的验证见 docs/RESULTS.md「hi3 与 WAu12：片段」。

### 功能：可靠性指标

LEOS 与 R(%) 对应的指标是 R_site 和 R_elem。定义：解完最小费用流后，强制
某原子（R_site）或某元素总账（R_elem）的电子数 ±1，全局重解要多付的代价
记为 m，按 R = 100·m/(m+1) 归一化到 0–100。全胞的 R_site 取所有原子中的
最小值。

| | R(%) | R_site / R_elem |
|---|---|---|
| 度量对象 | 排序边界上两个 eff-AO 的占据数差 | 强制 ±1 后全局重解多付的代价，链式重排全部计价 |
| 归一化 | 100·min(1, max(0, λ_LO − λ_FU + 1/2)) | 100·m/(m+1) |
| 粒度 | 全胞一个数，两个自旋通道取小 | 两级：R_site 计任何 ±1 重排；R_elem 只计改变某元素总账的重排 |
| 位点简并的表现 | R = 50，例如 GaSe 按原子划分（p.7080） | R_site 接近 0，例如 As4S4 1.8、GaSe 0.2 |
| 已知局限 | 不是概率（p.7078） | 输入有系统性错误时，可以 R 高而判定错，需要和 stray、ghost、dNe 一起看 |

stray 指未定域的 WF：spread 超过 3 Å²，或中心离所有原子都超过 1.6 Å。
dNe 是 Wannier 表示中的总电子数与 DFT 值之差。R_site 能单独反映同核键
对半、等价位点换手这类位点简并；R_elem 只看元素总账，与 74 案例基准按
元素平均计分的方式对应。

74 案例基准中两组案例的 R_elem 中位数：

| 分组（status 列） | 案例数 | R_elem 中位数 |
|---|---|---|
| 判对（OK） | 64 | 83.0 |
| 偏差（FAIL 或 ~） | 8 | 25.55 |

数据：data/v3/reliability_74.tsv，按 status 列分组；As4S4、GaSe 的 R_site
取同一文件的 As4S4、sega 行。定义和读法见 docs/RESULTS.md「可靠性指标 R
的表现」。

### 功能：电子化物

两个方法识别游离电子的方式不同。

| | Comas-Vilà | LEOS |
|---|---|---|
| 依据 | 密度中的 NNA，作为额外的盆参与分配 | stray：WF 中心落在某个原子附近，但 spread 异常大 |
| 高压 K | 找到 2 个 NNA，每个分到 2 个电子，K 为 +1，R 低（p.7079） | K 平均 0，游离电子记在 K 名下，带 stray 警告 |
| Ca3AuN | 没找到 NNA，判 Au(−3)，R = 60.4（p.7079） | Ca₂(+2)、Ca(0)、Au(−1)、N(−3)；6 个 WF 的 spread 为 5–10 Å²，正常值低于 2 Å² |
| 已知问题 | NNA 的有无和个数对数值参数敏感，见 Ca2N（p.7079） | 游离电子记在哪个原子名下是任意的，Ca3AuN 的 R_site = 7.1 |

Ca3AuN 的两种图像：按原文转述，IUPAC 技术报告（Karen et al.）指出，沿用
其他钙钛矿的做法会得到 Au(−3)，这违反 12−N 规则，因此报告认为
Au(−1) + 2e⁻ 的电子化物图像更准确（p.7079）。原文的 EOS 分析判 Au(−3)。
LEOS 给出 Au(−1)，2 个游离电子记在一个 Ca 名下，所以这个 Ca 显示 0 价。
两边的可靠性读数都低：原文 R = 60.4，LEOS R_site = 7.1。

数据：data/v3/reliability_74.tsv 的 ca3aun 行；详见 docs/RESULTS.md
「Ca3AuN：电子化物」。

## 40 固体对照结果

LEOS 对原文全部 40 个固体都有判定。多数体系的元素平均氧化态与原文一致；
不同之处集中在电子化物（Ca3AuN、高压 K）、高压 XeO3 和一个 HgF4 结构。

### 测试集构成与参考依据

原文把 40 个体系分成以下几组（p.7078）：

| 编号 | 体系 | 原文的比较对象 |
|---|---|---|
| 1–15 | IUPAC 技术报告（Karen et al.）中的例子：金属氧化物、钙钛矿、反钙钛矿、过渡金属体系 | IUPAC 的 iBOS 或 BVS 算法。作者报告全部一致（p.7081、p.7084） |
| 16–23 | Bi 和 Zurek 的 MH4 金属氢化物 | 原文献的指认 |
| 24–33 | Xe 化合物：高压氧化物 24–29、氟化物 30–32、高压 MgXe 33 | 原文献的指认 |
| 34–36 | 高压 HgF4 | 原文献的指认 |
| 37–40 | 同核键体系：CsH3、CsF3、CsF5、LiN5 | 原文献的指认 |

16–40 的原文献指认基于部分电荷、ELF 形状或几何判据（p.7078）。例如 MH4
依据 QTAIM 电荷和 H–H 键长（p.7081），Xe 氧化物依据 ELF 形状（p.7082），
HgF4 依据 PDOS（p.7083）。

后 25 个体系没有 IUPAC 参考答案。LEOS 与原文在这些体系上的不同，只能记为
两种方法给出的结果不同，不能据此判断对错。

### LEOS 判定的数据位置

- 30 个体系在 data/v3/comas_gap/<case>/，都是完整 MLWF 计算。
- 10 个体系与 74 案例基准重合，判定取 data/v3/reliability_74.tsv：cro3、
  kmgf3、nabf4、srbe_oh_4、cs3auo、rb3auo、ca3aun、sega（GaSe）、wcl4、
  rb2cucl4。其中 9 个为 OK。ca3aun 为 FAIL，因为基准参考答案取 Au(−3)
  （pyleos/tests/reference_os.json）。WCl4 用的是修正后的结构
  data/v3/m4_wave2/wcl4。
- data/v3/comas_gap/ 另有 6 个 scdm_* 目录，是 SCDM-only 快速模式的运行
  （INCAR 中 EOS_SCDM_ONLY = .TRUE.），不作为对照结果。其中 scdm_cse2、
  scdm_mocl2、scdm_wcl5、scdm_wcl6 不在原文测试集内；scdm_sega、scdm_wcl4
  与 GaSe、WCl4 同名，scdm_wcl4 的胞实为 WCl3（docs/testset_registry.tsv）。
- data/v3/comas_gap/summary.tsv 有 36 行（30 + 6），status 全为 OK。这一列
  是批处理运行状态，表示 EOSCAR 正常产出（pyleos/batch.py 的 triage），
  不表示判定与参考答案一致。
- 结构来源：SI 逐表转录、实验 CIF、数据库同相结构共 29 个新结构，另加
  docs/structures_fix/ 中的高压精确相（MgXe 与 Mg/Sr/Ce/Th 四个氢化物）。
  每个结构的来源与体检结论登记在 docs/testset_registry.tsv；其中 cs2pt、
  csh3、lin5 标为 SUSPECT。
- 判定设置：β = 10，Allen 电负性。comas_gap 下 30 个目录在默认代价
  （几何 + 电负性）和投影辅助判定（TAU、PMARGIN 取默认值）下给出相同的
  元素平均值。10 个重合体系取 74 案例基准的默认代价结果。

### 逐类对照

原文一栏的 R 均指 R(%)。LEOS 一栏的写法沿用 docs/RESULTS.md「Comas-Vilà
对照集（40 个固体）」。带 OK 的行来自 data/v3/reliability_74.tsv：
OK 表示 dev = 0.00，即各元素平均氧化态等于基准参考答案，表中数值取自
pyleos/tests/reference_os.json。

| 编号 | 体系 | 原文判定（页码） | LEOS 判定 | LEOS 数据 |
|---|---|---|---|---|
| 1 | CrO3 | Cr +6，R = 71（p.7078–7079） | Cr+6 / O−2，OK | reliability_74.tsv |
| 2–4 | KMgF3、NaBF4、SrBe(OH)4 | 与 IUPAC 一致，R 均为最大值（p.7079） | K+1/Mg+2/F−1、Na+1/B+3/F−1、Sr+2/Be+2/O−2/H+1，均 OK | reliability_74.tsv |
| 5–6 | Cs3AuO、Rb3AuO | 与 IUPAC 一致，R = 100；Rb3AuO 为金属，AOM 改用 Bloch 态（p.7079） | Cs+1/Au−1/O−2、Rb+1/Au−1/O−2，均 OK | reliability_74.tsv |
| 7 | Ca3AuN | Au −3，R = 60.4，未找到 NNA（p.7079） | Ca₂(+2)、Ca(0)、Au(−1)、N(−3)；对基准参考答案记 FAIL | reliability_74.tsv |
| 8 | K（高压） | 2 个 NNA 各分 2 个电子，K +1，R 低（p.7079） | K 平均 0，游离电子记账伴 stray | comas_gap/k_hp |
| 9 | GaSe | 按原子划分 R = 50，平摊后 Ga +2；Ga–Ga 片段 R = 77.7（p.7080） | Ga+2 / Se−2，OK | reliability_74.tsv 的 sega 行 |
| 10 | Ba3Si4 | Si4 作片段：Ba +2、Si4 −6，R = 100（p.7080） | Ba+2 / Si−1.5，Si4⁶⁻ 蝴蝶簇 | comas_gap/ba3si4 |
| 11 | Cs2Pt | Pt −2，判定无争议（p.7080） | Cs+1 / Pt−2 | comas_gap/cs2pt |
| 12 | WCl4 | W +4，R = 92.6（p.7080） | W+4 / Cl−1，OK | reliability_74.tsv；m4_wave2/wcl4 |
| 13 | Rb2CuCl4 | Cu +2，β 通道 R = 75.4（p.7080–7081） | Rb+1/Cu+2/Cl−1，OK | reliability_74.tsv |
| 14–15 | YBaFe2O5 两相 | 低温相 Fe1 +3、Fe2 +2，R = 55.6；高温相两个 Fe 平摊为 +2.5，R = 61.3（p.7081） | 两相 Y+3 / Ba+2 / Fe+2.5 / O−2；元素级一致，位点电荷序另论 | comas_gap/ybafe2o5_a、_b |
| 16–18 | MgH4、CaH4、SrH4 | M +2，Hb −1；Ha–Ha 作片段时为中性 H2，比按原子划分更明确（p.7081） | M+2 + 2H⁻ + H2；H2 单元内部为位点简并 | comas_gap/mgh4、cah4、srh4 |
| 19–20 | ScH4、YH4 | M +3，Hb −1；最后的电子平摊给等价的 Ha，各 −½（p.7081–7082） | M+3 + 3H⁻ + ½H2 | comas_gap/sch4、yh4 |
| 21–23 | ZrH4、CeH4、ThH4 | M +4，四个 H 均为 −1；按原子划分 R 更高，ZrH4 为 75.4（p.7082） | M+4 + 4H⁻ | comas_gap/zrh4、ceh4、thh4 |
| 24–25 | XeO | 100 GPa：Xe +2，R = 74.6；200 GPa：Xe1 +4、Xe2 0，R = 73.3（p.7082） | 两相 +2 | comas_gap/xeo_a、xeo_b |
| 26–27 | XeO2 | 150 GPa：Xe +4，R = 77.6；200 GPa：Xe +4，R = 58.7（p.7082） | 两相 +4 | comas_gap/xeo2_a、xeo2_b |
| 28 | XeO3（130 GPa） | 按原子划分：Xe +4、O1 −2、O2 0，R = 60.9；把 O2 作片段得混合价 Xe +5，但 R 更低（p.7082–7083） | +5，伴 stray | comas_gap/xeo3_a |
| 29 | XeO3（150 GPa） | Xe1 +6、Xe2 +4，R = 51.2，作者认为本质上模糊（p.7083） | +6 | comas_gap/xeo3_b |
| 30–32 | XeF2、XeF3、XeF4 | +2；XeF3 为 XeF2 与 XeF4 分子的混合，Xe1 +2、Xe2 +4；+4。R 在 89–98 之间（p.7083） | XeF2 +2 / XeF3 +3 / XeF4 +4。XeF3 的 +3 是元素平均，逐位点为 +2 和 +4 | comas_gap/xef2、xef3、xef4 |
| 33 | MgXe | Mg +2、Xe −2，R = 59.6（p.7083） | Mg+2 / Xe−2（200 GPa）；晶格为状态方程推定值，误差约 ±1.5% | comas_gap/mgxe |
| 34–36 | HgF4 | 三个结构同为 I4/m，均为 Hg +4；R 随压力升高，34 号为 56.3（p.7083） | hgf4_a、hgf4_c 为 +4；hgf4_b 的 Hg 平均 +3（+2 到 +4） | comas_gap/hgf4_a、_b、_c |
| 37 | CsH3 | H3 作片段为 −1；按原子划分时中心 H 为 +1，R = 68（p.7083） | H−1/3；结构存疑，见表后说明 | comas_gap/csh3 |
| 38 | CsF3 | F3 作片段为 −1；按原子划分 R = 66（p.7083） | F−1/3，即 F3⁻ | comas_gap/csf3 |
| 39 | CsF5 | F5 作片段为 −1；按原子划分为 Cs +1、F −1、F2 0（p.7084） | F−0.2，即 F5⁻ | comas_gap/csf5 |
| 40 | LiN5 | 五个 N 占据数几乎相同，R = 50，无法判定（p.7084） | N−0.2；结构存疑，见表后说明 | comas_gap/lin5 |

### 一致与不同

元素平均氧化态一致的体系：

- 前 15 个中除 Ca3AuN、高压 K 外的全部体系。YBaFe2O5 只在元素平均层面
  一致，位点电荷序另论。
- MH4 氢化物：金属价态和 H 的元素平均都与原文一致。ScH4、YH4 的位点分布
  不同：原文是 2 个 Hb(−1) 加 2 个 Ha(−½)，LEOS 是 3 个 H⁻ 加 1 个中性 H。
- XeO、XeO2、XeF2/XeF3/XeF4、MgXe，以及 HgF4 的两个结构。XeF3 在位点
  层面也一致，两类 Xe 分别为 +2 和 +4（data/v3/comas_gap/xef3）。
- 贪心基线在 SrH4、CeH4、ThH4 上给出金属 +10、+9、+10（EOSCAR 的贪心
  列），运输分配给出 +2、+4、+4。

不同的体系：

| 体系 | 原文 | LEOS | 说明 |
|---|---|---|---|
| Ca3AuN | Au −3，R = 60.4 | Au −1 加 2 个游离电子 | 见上文"功能：电子化物" |
| 高压 K | K +1，电子记在 NNA 上 | K 平均 0 | LEOS 没有 NNA 这类记账单元，游离电子记在 K 名下，并给出 stray 警告 |
| XeO3 130 GPa | Xe +4，R = 60.9 | +5 | 原文把 O2 作片段时得到的备选解也是 Xe +5，但 R 更低 |
| XeO3 150 GPa | Xe1 +6、Xe2 +4，R = 51.2 | +6 | 原文认为该相判定本质上模糊 |
| HgF4（hgf4_b） | +4 | Hg 平均 +3 | 原文三个结构同为 I4/m；本项目三个结构与原文的对应关系未核实，登记表把 hgf4_c 记为 Pm-3m 型骨架 |

另有两点会影响对照：

- CsH3 和 LiN5 的 LEOS 胞内容与原文描述不同。原文是线性 H3 单元
  （p.7083）和五唑阴离子 N5⁻（p.7084）。本项目 CsH3 胞中的 H 组成两个孤立
  的等边三角形 H3；LiN5 胞中是线性叠氮 N3 单元加拉长的 N2，不是五唑环
  （docs/testset_registry.tsv、docs/testset_findings.md）。这两行的数值
  不是同一结构上的比较。
- 用 pyleos 复算 comas_gap 下 30 个完整 MLWF 目录，12 个带 stray 警告，
  按 stray 个数：ba3si4 36、xeo3_b 9、ceh4 8、thh4 8、k_hp 6、mgxe 5、
  ybafe2o5_b 4、csh3 3、mgh4 3、sch4 3、lin5 1、xeo3_a 1。

复现：

```bash
# 单个 comas_gap 目录：逐元素判定、stray/ghost 警告、R_elem/R_site
python3 pyleos/cli.py case data/v3/comas_gap/xeo3_a
# 74 案例基准（含 10 个重合体系）；必须显式给 --out，否则会覆盖基准参考表
python3 script/full_bench.py --out /tmp/bench.tsv
```

## 原文位置速查

| 内容 | 位置 |
|---|---|
| 分子版 EOS（Ramos-Cordoba et al. 2015）；适用于 Hilbert 空间与实空间划分 | Introduction 末与 Theory 开头，p.7076 |
| eff-AO、EFO 与 AOM 的定义，式 1–4 | Theory，p.7077 |
| 全局排序算法；不逐个取整 | Theory，p.7077 |
| R(%) 定义，式 5；R 低于 50 的含义 | Theory，p.7077–7078 |
| 简并时平摊电子；R 不是概率；式 5 形式有任意性 | p.7078 |
| 程序链：PAW 重建、norm-conserving 计算、MLWF、QTAIM 下的 AOM；片段 AOM 相加 | Results 开头，p.7078 |
| 固体版 R 允许边界轨道同属一个原子（超胞论证） | p.7078 |
| 测试集分组与参考依据；SI Table S1–S16 | p.7078 |
| 40 个体系的化学式与结构来源 | Table 1，p.7079 |
| CrO3：+6，QTAIM 电荷 +1.76，R = 71 | p.7078–7079 |
| 钙钛矿 2–4；反钙钛矿 5–6；Rb3AuO 为金属，改用 Bloch 态 | p.7079 |
| Ca3AuN：Au(−3)，R = 60.4，未找到 NNA；Karen et al. 的 Au(−1) + 2e⁻ | p.7079 |
| 高压 K 的 NNA；Ca2N 的 NNA 敏感性；ELF + Voronoi 划分的建议 | p.7079 |
| GaSe：R = 50，Ga–Ga 片段 R = 77.7 | p.7079–7080 |
| Ba3Si4、Cs2Pt、WCl4（R = 92.6） | p.7080 |
| Rb2CuCl4（β 通道 R = 75.4） | p.7080–7081 |
| YBaFe2O5（R = 55.6、61.3）；前 15 个与 iBOS/BVS 一致 | p.7081 |
| MH4 氢化物与片段划分的选择 | p.7081–7082 |
| XeO、XeO2（R = 74.6、73.3、77.6、58.7） | p.7082 |
| XeO3 两相（28 号 R = 60.9；29 号 R = 51.2） | p.7082–7083 |
| XeF2/XeF3/XeF4（R 89–98）；MgXe（R = 59.6）；HgF4；CsH3、CsF3 | p.7083 |
| CsF5；LiN5（R = 50） | p.7084 |
| 结论：可推广到其他实空间划分；R 用于自评适用范围；与 iBOS/BVS 一致 | Conclusions，p.7084 |
| 计算参数与所用程序 | Computational Details，p.7084 |
