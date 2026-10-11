# LEOS 术语表

本仓库文档统一使用下列术语。方法说明见 docs/PLAN.md，结果见 docs/RESULTS.md。

## 基本概念

**氧化态（OS，oxidation state）**
假想把每条键的电子全部判给电负性更大的一端后，原子相对中性时多出或缺少的
电子数，是整数。它是记账约定，不是可测量的物理量。
不要与化合价（valence，成键能力）或电荷（Bader 等方法给出的连续值）混用。

**Wannier 函数（WF）**
由 DFT 占据态变换得到的局域电子团，有中心坐标和 spread。闭壳层计算
（ISPIN = 1）每个 WF 装 2 个电子；自旋极化计算（ISPIN = 2）每个自旋通道的
每个 WF 装 1 个电子。

**Wannier 中心**
WF 的中心坐标（笛卡尔坐标，Å），是判断归属的几何依据。中心紧贴某原子时通常
是孤对电子，落在两个原子之间时是共价键。

**spread**
WF 的空间展宽（Å²），衡量 Wannier 化的质量。定域良好的 WF 一般在 1–2 Å² 以下。

**EOSCAR**
LEOS VASP 补丁输出的文件，pyleos 的输入。包含：逐原子表（贪心结果列
OS(int)，运输结果列 OS(flow)）、#WF 块（每个 WF 的中心、spread、占据数）、
#WP 块（每个 WF 在各原子 PAW 投影上的权重），以及电子丢失诊断块。

## 分配

**贪心（greedy）**
对照方法：每个 WF 独立判给最近（考虑电负性）的原子，WF 之间互不知情。可能
给出物理上不可能的结果，例如一个 S 收到 14 个电子。结果在 EOSCAR 的
OS(int) 列。

**运输（最小费用流）**
LEOS 的分配方法：把所有 WF 的电子一起分配给原子，每个原子有接收上下限，
求总代价最小的方案。电子数和容量都是整数，所以最优解是整数。论文中首次出现时
写作 "minimum-cost-flow / optimal-assignment formulation"，以免被理解为电子
输运。

**OS(flow)**
EOSCAR 逐原子表的最后一列，是 VASP 补丁内部完成的运输分配结果（INCAR 参数
EOS_FLOW，默认开），与 pyleos 默认设置的结果逐原子相同。

**代价**
把一个 WF 的电子给某个原子的代价：d² + β·(1 − χ/χmax)。d 是 WF 中心到原子的
距离，χ 是电负性（默认 Allen 标度），β 默认 10。

**β**
代价中电负性项的权重，默认 10。8–12 之间结果稳定；太小时 X–H 键会被判错且
没有警告。

**容量**
每个原子接收电子数的上下限 [n_min, n_max]，按价层物理极限设定，只排除不可能
的情况，不排除罕见价态。不按"化学上常见的氧化态范围"设定。

**半芯**
必须记满的内层壳层（如 Na 的 2s2p），对应 n_min 中必须满足的部分。

**投影辅助判定**
可选的判定规则：某个 WF 的前两个候选原子代价差小于 TAU（1.0）、且 #WP 投影
权重领先超过 PMARGIN（0.20）时，按投影改判。默认不开启。不适合含大量 C–H 键
的有机体系。代码中叫 `tier2="arbiter"`；旧文档中的"档 2""仲裁"指的就是它，
"档 1"指默认代价。

**片段（fragment）**
用户指定的一组原子，作为一个整体接收电子（容量相加，代价取离 WF 最近的成员），
输出整体的氧化态。用于离域阴离子和金属簇（I3⁻、Au12、Co4 等）。

**位点简并**
等价原子之间交换电子几乎不改变总代价的情况，例如同核键两端、簇内等价位点、
等价 O2 单元。此时元素平均可信，单个原子的值不唯一，R_site 接近 0，应看元素
平均或片段总价。

**可靠性指标 R（R_site、R_elem）**
解出最优分配后，强制某原子多拿或少拿 1 个电子、全局重新求解所需多付的代价 m，
按 R = 100·m/(m+1) 归一化到 0–100。R_site 计入所有扰动，同核键两端等简并会
使它接近 0；R_elem 只计入改变某元素总电子数的扰动。R 是 Comas-Vilà R(%) 的
对应物。R 低表示解不唯一，R 高不保证正确。

## 质量警告

**stray**
未定域的 WF：spread > 3 Å²，或中心离所有原子都超过 1.6 Å。说明 Wannier 化本身
有问题，分配无法修正，应重做 Wannier 化。

**ghost**
WF 数多于占据的电子对时，多出来的、没有电子可装的 WF（常见于金属和分数占据
体系）。总电子数上限（NELECT）使求解器丢掉代价最高的部分。ghost > 0 时电正性
金属的价态可能被高估。

**dNe**
Wannier 表示中的总电子数与 DFT 真实电子数之差。健康计算低于 1e-6；明显不为零
说明有电子没有进入 Wannier 子空间（电子丢失）。

**诊断版 EOSCAR**
电子丢失超过 EOS_NELECT_ABORT 时，补丁中止前写出的 EOSCAR，带
`# LEOS_ABORTED` 标记和丢失定位块（窗口型 #NL、弥散型 #NLD）。wannier90 自身
失败时不产生这个文件。

## 基准与数据

**74 案例基准**
三批测试体系：历史回归集 36、扩展第一波 13、扩展第二波 25，其中 72 个有参考
答案。

**参考答案**
按化学常识和 IUPAC 规则人工标注的氧化态，存于 pyleos/tests/reference_os.json。
混合价体系标元素平均值（Fe3O4 中 Fe 为 +2.667）。

**dev**
每种元素的平均氧化态与参考值之差的最大绝对值。dev ≤ 0.25 为通过。只检查元素
平均，不检查单个原子。

**基准参考表**
data/v3/reliability_74.tsv。CI 把新跑出的基准结果与它逐行比较（status 和
dev 两列）。gpuw90 中与之对应的是 wannier90 参考输出（.wout）。

**Comas-Vilà**
参考方法：Comas-Vilà et al., JCTC 2025，用 QTAIM 原子域和原子重叠矩阵判定
有效氧化态，测试集为 40 个固体。"EOS"这个名字出自该文。

**v1 / v3**
数据版本。v1 是项目早期程序（贪心）产生的历史数据（data/v1/）；v3 是当前
补丁产生的数据（data/v3/）。

**M4 / M5**
数据批次编号。M4 为扩展基准（data/v3/m4_wave1、m4_wave2）；M5 为应用：
M5-A 为 Li2O2，M5-B 为富锂正极三体系（data/m5b/），M5-C 为 CeO2。

## Wannier 流程

**原位 / 离线**
原位：在 VASP 内完成 Wannier 化和分配（贪心列和 OS(flow) 列），可用于分子
动力学逐帧分析。离线：在 VASP 之外用 pyleos、gpuw90 或独立的 wannier90.x 处理。

**局域化（wannierise）**
Marzari–Vanderbilt 方法：在给定子空间内做幺正变换，使总 spread 最小，得到
Wannier 中心和 spread。

**解纠缠（disentanglement）**
能带数多于 WF 数（num_bands > num_wann，金属或含导带的体系）时，在局域化之前
先选出要局域化的子空间（Souza–Marzari–Vanderbilt 方法）。跳过这一步直接局域化
会得到错误结果且没有报错。

**冻结窗口**
解纠缠时必须原样保留的能量范围，通常覆盖占据态，保证价电子不因子空间选取而
丢失。

**自动冻结窗口（AUTO_WINDOW）**
补丁按占据态自动设定冻结窗口的功能（INCAR 参数 LWANNIER90_AUTO_WINDOW），
不需要手写 dis_froz_max。手写输入缺少窗口是 NELECT mismatch 失败的主要原因。

**SCDM**
Wannier 化的初猜方法（selected columns of the density matrix），不需要手工设定
投影。SCDM-only 是跳过迭代局域化、只用初猜中心的快速模式。

**引导中心（guiding centres）**
热启动的实现方式：把上一帧的 Wannier 中心写进下一帧的 wannier90 输入。

**gpuw90**
本仓库中用 PyTorch 实现的解纠缠和局域化工具，CPU 和 CUDA 用同一份代码，不修改
wannier90 和 VASP。

**失败分类（FAIL_DFT / FAIL_W90 / FAIL_LEOS）**
batch.py 对失败算例的分类：FAIL_DFT 为 SCF 崩溃或 OUTCAR 不完整；FAIL_W90 为
DFT 正常但 wannier90 没有完成；FAIL_LEOS 为 Wannier 化正常但没有可用的 EOSCAR
（缺失或带 LEOS_ABORTED 标记）。每类有各自的自动重试设置。
