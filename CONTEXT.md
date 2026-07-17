# LEOS

LEOS（Lattice Effective Oxidation States）：从 DFT 计算的 Wannier 函数自动判定
固体中每个原子的整数氧化态。本文件是项目的术语表（ubiquitous language）——
只定义概念，不记实现细节；实现与数据见 docs/PLAN.md、docs/RESULTS.md。

## 核心概念

**氧化态（OS, oxidation state）**:
给原子记账的整数：假想把每条键的电子全部判给电负性更大的一端后，原子相对
中性状态的电子盈亏。是记账约定，不是可测量的物理量。
_Avoid_: 化合价（valence，指成键能力）、电荷（charge，Bader 等连续量）

**Wannier 函数（WF）**:
把弥散的 DFT 电子波变换成的一团局域化电子云，有明确的中心坐标和 spread。
闭壳层（ISPIN=1）每个 WF 装 2 个电子；自旋极化（ISPIN=2）每通道每 WF 装 1 个。

**Wannier 中心（centre）**:
WF 的中心坐标（笛卡尔，Å）。归属判定的几何依据：落在原子近旁通常是孤对，
落在两原子中间是共价键电子对。

**spread（弥散度）**:
WF 的空间弥散程度（Å²），Wannier 化质量的量度；正常应低于 1 Å² 量级。

**EOSCAR**:
LEOS VASP 补丁写出的输出文件：Wannier 中心、spread、占据数、贪心判定结果，
以及 #WF（per-WF 明细）、#WP（PAW 投影权重）数据块。pyleos 的唯一主输入。

**stray WF**:
质量不合格的 WF：spread > 3 Å²，或中心离所有原子超过 1.6 Å。标记为 stray
表示 Wannier 化本身失败，分配层无法修复，应触发工作流层重跑。
_Avoid_: 坏函数（口语）、unlocalised WF

**幽灵电子（ghost）**:
WF 总数比实际占据电子对多出的部分（体系有分数占据时出现）。用总运输量
上限（NELECT）排除，求解器自动丢弃最贵的多余货物。

**dNe（电子守恒偏差）**:
Wannier 表示中的总电子数与 DFT 真实值之差。健康计算低于 1e-6；明显非零
说明有电子未进入 Wannier 子空间（电子丢失）。与 stray、ghost 并列为三大
质量警告。

## 分配算法

**运输（transport，运输形式化）**:
本项目的核心方案：把"每团电子归哪个原子"表述为带容量约束的最小费用流
（运输问题），所有电子一起做全局决策；整数供给 + 整数容量保证最优解自动
是整数。论文首现处须写 "minimum-cost-flow / optimal-assignment formulation"，
以免凝聚态读者误解为电子输运。
_Avoid_: optimal transport（不加限定语时歧义）、电子输运

**OS(flow)**:
v3.1 起 EOSCAR 逐原子表行尾的原位运输判定列：VASP 补丁内部直接完成最小
费用流分配（档 1 口径，开关 EOS_FLOW 默认开），与 pyleos 结果逐原子全等。
正式结论与 MD 轨迹逐帧判定不再依赖后处理。

**贪心（greedy）**:
对照基线：每团电子独立找最近（或电负性优先）的原子，互相不通气。VASP 补丁
内置的原位路径，会产生物理上不可能的结果（如 S 收 14 个电子）。

**宽物理界（wide physical bounds）**:
原子容量 [n_min, n_max] 的设置原则：按价电子层的物理极限设定，只排除不可能，
不排除罕见（高压相、新奇价态必须保持可行）。
_Avoid_: 化学先验 / 紧容量（那是另一种尚未实现的诊断模式，不是默认口径）

**半芯（semicore）**:
必须保持满占据的内层壳层，对应容量下限 n_min 的强制部分。

**beta（β，电负性权重）**:
运费公式 d² + β·(1 − χ/χmax) 中电负性折扣的权重，默认 10.0（β 过低会静默
错判 X-H 键，全量 β 扫描的安全平台为 [8,12]，默认取平台腹地 10——
data/v3/beta_scan_20260712.tsv）。电负性 χ 默认用 Allen 标度
（IUPAC 2014/2016 配套）。

**档 1 / 档 2（tier 1 / tier 2）**:
档 1 = 纯几何 + 电负性的运费（正式口径，全部零参数结论出自它）；
档 2 = 追加 #WP 投影成分信息，定位为显式诊断工具，默认关闭。

**仲裁（arbiter）**:
档 2 的选定模式：仅当档 1 的前两名归属代价差小于阈值 TAU 时，才让 #WP
成分信息介入裁决（参数 (TAU, PMARGIN)，默认 (1.0, 0.20)）。含大量 C-H 键
的体系不适用（H 1s 在键中点权重偏大，会翻错）。

**片段（fragment）**:
用户指定的一组原子合并成的单一仓库（容量相加、运费按最近成员），输出片段
总氧化态。用于离域阴离子与金属簇（I3⁻、Au12、Co4）等"每团电子归一两个
原子"图像失效的场景。

**可靠性指标（R_site / R_elem）**:
每笔归属的"票差"：解完最小费用流后，强制某原子（或片段）的电子数 ±1、
全局重解需要多付的运费（链式重排全部计价）。R_site 计任何 ±1 重排（同核
键对半等本征简并会拉低它）；R_elem 只计改变某元素总账的重排。Comas-Vilà
R(%) 指标的对应物。

**位点简并（site degeneracy）**:
等价原子之间互换电子归属的代价近零的局面（同核键两端、簇内等价位点、
等价 O2 单元）：元素平均仍唯一可信，但单原子归属不唯一，R_site 给出近零
分。此时该下结论的层级是元素平均或片段总价，不是单个原子。

## 基准与数据

**74 案例基准**:
三批测试体系（历史回归 36 + 扩展第一波 13 + 扩展第二波 25），其中 72 个
有参考答案（计分分母 72）。金标文件 data/v3/reliability_74.tsv，CI 对它做
逐案例 diff 防回归。

**dev（偏差指标）**:
对每个体系：每种元素的氧化态平均值与参考值之差，取所有元素中最大的绝对
偏差。dev ≤ 0.25 记为通过。是元素平均口径，不检验单个原子。

**参考答案（reference OS）**:
按化学常识和 IUPAC 规则人工标注的理想氧化态（pyleos/tests/reference_os.json）。
混合价体系标元素平均值（Fe3O4 的 Fe = +2.667）。

**金标（golden reference）**:
作为回归对照的既定正确输出。两处用法：74 案例的 reliability_74.tsv（分配
层）；wannier90.x 的参考 .wout（gpuw90 数值验收，中心偏差阈值 1e-4 Å）。

**Comas-Vilà（参考方法）**:
JCTC 2025 的重型对照方法（QTAIM 原子域 + 原子重叠矩阵，四程序串联），
"EOS"名字的出处；其 40 固体测试集是扩展基准的蓝本，已全量覆盖
（data/v3/comas_gap/）。

**v1 / v3**:
数据版本：v1 = 项目接手前贪心时代的历史数据；v3 = 当前补丁（leos_v3.patch）
产出的现行数据。

**M4 / M5（波次）**:
数据生产的里程碑波次编号（如 data/v3/m4_wave2）。M5-B 指论文应用的富锂
正极三体系集（Li2MnO3 / Li2RuO3 / NMC）。

## Wannier 管线

**原位 / 离线（in-situ / offline）**:
原位 = VASP 补丁内部完成 Wannier 化与判定——快速贪心（对照基线）加 v3.1
起的原位最小费用流 OS(flow) 列（EOS_FLOW 默认开，服务分子动力学逐帧
判定）；离线 = pyleos / gpuw90 / 独立 wannier90.x 在 DFT 之外做后处理。

**wannierise**:
Marzari–Vanderbilt 幺正最小化：在固定子空间内旋转规范，最小化总 spread，
得到 Wannier 中心与 spread。

**解纠缠（disentanglement）**:
纠缠带组（num_bands > num_wann，金属/含导带体系）在 wannierise 之前的
子空间选取步骤（Souza–Marzari–Vanderbilt）。把全带直接喂给 wannierise
会得到静默错误结果。

**孤立带组 / 纠缠带组**:
num_bands == num_wann 的带组直接 wannierise；num_bands > num_wann 的带组
必须先解纠缠。

**冻结窗（frozen window）**:
解纠缠中强制原样保留的能量窗口（占据态所在），保证价电子物理不被子空间
选取改变。

**AUTO_WINDOW（自动冻结窗）**:
v3 补丁按占据态自动设定解纠缠冻结窗口的机制，免手写 dis_froz_max。旧版
手写输入缺窗口导致的电子丢失（NELECT mismatch 类失败）由它根治。

**SCDM**:
Wannier 化的初猜方法（selected columns of the density matrix）。SCDM-only
指跳过迭代最小化、只用 SCDM 初猜的快速筛查模式。

**引导中心（guiding centres）**:
轨迹模式热启动的实现：把上一帧的 Wannier 中心写进下一帧的输入作初猜引导。

**gpuw90**:
本仓库的独立 GPU Wannier 化工具（PyTorch，CPU/CUDA 同码）：Phase 1 =
wannierise（孤立带），Phase 2 = disentangle（纠缠带）。不 fork wannier90、
不改 VASP。

**abort-EOSCAR**:
丢失电子超过熔断阈值（EOS_NELECT_ABORT）时补丁先写出的诊断版 EOSCAR
（`# LEOS_ABORTED` 标记 + 窗口型 #NL / 弥散型 #NLD 丢失定位块）再中止，
保证失败也可定责。注意 wannier90 自身的硬失败（FAIL_W90）不产生此文件。

**失败分诊（triage：FAIL_DFT / FAIL_W90 / FAIL_LEOS）**:
批跑对失败案例的三级归因：FAIL_DFT = SCF 崩溃或 OUTCAR 不完整；
FAIL_W90 = DFT 正常但 wannier90 未收尾；FAIL_LEOS = Wannier 化正常但无
可用 EOSCAR（缺失或带 LEOS_ABORTED 标记）。每类各有自动重试梯
（retry ladder）。
