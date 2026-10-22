# 测试集结构核查

**状态（2026-07-12 核查，issue #19；2026-10-21 按 VASP 实际读入的几何复核）：**
四个钨氯化物换用正确结构后已重算，MgXe 和四个 MH4 氢化物换用了文献中的精确
结构。复核又查出 li2o2、re2 两个 v1 案例和 Li2O2 脱锂超胞的几何是坏的，见下文
"几何复核"。其余问题案例保留原结构，判定按胞的实际内容解读。逐案例记录见
[testset_registry.tsv](testset_registry.tsv)；每个结构的来源、判定和依据统一登记在
[cases/STRUCTURES.tsv](../cases/STRUCTURES.tsv)。

## 核查方法

对全部 113 个测试单位（74 案例基准、Comas-Vilà 对照集及 SCDM 对照副本），由
12 个并行的自动核查代理逐个读取结构文件，比对化学式，核对参考答案，并检索
文献中的已知晶体相。

| 判定 | 当前数量 | 核查时数量 | 含义 |
|---|---|---|---|
| OK | 89 | 86 | 组成相符，来源可信，或是常规的理想化近似 |
| SUSPECT | 14 | 15 | 来源不明或化学上存疑：手搭团簇、假想多形体、与已知相有出入 |
| WRONG | 10 | 12 | 结构与化学式或已知相明确冲突，或几何本身坏了 |

当前数量是登记表现在的判定：WCl3–6 换用修正结构后由 WRONG 改为 OK，li2o2 由
SUSPECT、re2 由 OK 改为 WRONG（见"几何复核"）。

登记表各列：case、batch、dir、formula_expected（按案例名应有的组成）、
formula_cell（胞内实际组成）、stoich_match、known_phase（文献已知相）、
source_guess（推断的结构来源）、ref_os_entry（参考答案条目）、ref_os_ok
（参考答案是否电荷平衡）、verdict、note。

## 问题类型

### 1. CIF 只转录了不对称单元（6 例）

晶格参数与实验相符，但胞内只有不对称单元的原子，没有按空间群展开。

| 案例 | 胞内实际 | 应有 |
|---|---|---|
| wcl3 | 1 W + 3 Cl（R-3 相） | 24 原子（W6Cl18）。化学配比碰巧正确，体积 441 Å³/原子暴露了问题 |
| wcl4 | 1 W + 3 Cl（C2/m） | 4 W + 16 Cl |
| wcl5 | 2 W + 9 Cl（C2/m） | 48 原子 |
| wcl6 | 1 W + 1 Cl（W 在 3a、Cl 在 18f 各取一个） | 21 原子 |
| scdm_wcl4/5/6 | 与上面相同的结构 | |

### 2. 胞内物种与案例名不符（5 例）

| 案例 | 胞内实际 |
|---|---|
| CO（历史分子集） | 两个乙烯二酮分子 O=C–C=O（C–C 1.53 Å），密度约为固态 CO 的两倍 |
| hi3 | HI 分子与 I2 分子组成的范德华链（H–I 1.70 Å，I–I 2.78 Å），不是 I3⁻ |
| lin5 | 线性 N3 单元加一个拉长的 N2，不是五唑环 |
| rhcl_co__p_c6h5_3_2 | 多一个配位的 SO2，实际为 RhCl(CO)(SO2)(PPh3)2；参考答案按无 SO2 的化合物写 |
| AuClB 配合物 | 胞内 C42H32，与名称不符；一个苯环通过 C–Au σ 键与金相连 |

### 3. 手搭的理想化胞或假想多形体（8 例）

| 案例 | 问题 | 判定 |
|---|---|---|
| sno | 假想的岩盐相；SnO 实际只有层状 litharge 相 | WRONG |
| ti2o3 | 金红石 TiO2 去掉一个 O 的 5 原子缺陷模型，不是刚玉型 Ti2O3 | WRONG |
| cuo | 假想的岩盐多形体；CuO 实际为单斜黑铜矿，Jahn–Teller 配位被抹掉 | SUSPECT |
| ce2o3 | 理想化胞，Ce 配位数与 A 型相不符，c 轴偏长 7% | SUSPECT |
| reo2 | 取整参数的金红石模板，与基态单斜相和实验相都不符；参考答案缺失 | SUSPECT |
| n4s4 | 密度比已知单斜相低 25% | SUSPECT |
| cs2pt | c/a 比实验短 16.6%，体积小 12.7% | SUSPECT |
| csh3 | H 组成孤立的 H3 三角形，与文献中的线性图像不符 | SUSPECT |

### 4. 未弛豫的孤立团簇（5 例，历史含氧酸盐）

na2cro4、na3vo4、namno4、recl4 是手搭的团簇且未弛豫，键长偏差 26–35%（例如
Cr–O 2.08 Å，实际为 1.66–1.68 Å）。孤立的多价阴离子在气相中束缚不住电子，
高价态的电子结构不可靠。li2o2（历史氧化物集）是盒子里的孤立分子，不是晶体。
nareo4 同样是气相离子对，但键长正确，判为 OK。

### 参考答案缺失

na3vo4 和 reo2 在参考答案表中没有条目，因此不在计分案例中。

## 处理状态

| 案例 | 判定 | 处理 | 当前基准结果 |
|---|---|---|---|
| wcl3、wcl4、wcl5、wcl6 | WRONG | 按空间群展开或换用文献结构后重算，cases/wcl*/v3 中已是新结构 | WCl3、WCl4、WCl6 通过；WCl5 未通过（Wannier 化失败） |
| scdm_wcl4/5/6 | WRONG | 未重算（SCDM-only 对照用的副本） | 不计入 74 案例基准 |
| CO | WRONG | 换用 α-CO（mp-11875）重算，cases/co/v31 | 原运行与新运行都通过 |
| sno、ti2o3 | WRONG | 换用 romarchite（COD 9011217）与刚玉型 Ti2O3（COD 9008082）重算，cases/<体系>/v31 | 原运行与新运行都通过 |
| rhcl_co__p_c6h5_3_2 | WRONG | 保留 | 未通过（dev 4.00） |
| AuClB 配合物 | WRONG | 保留 | 无参考答案，不计分 |
| hi3、lin5、csh3 | SUSPECT | 保留：HI3 没有固相；LiN5、CsH3 要用 Comas-Vilà 原文 SI 中的预测结构才能对照，COD 与 MP 都没有（CsH3 现有胞即 mp-1006815） | hi3 未通过；lin5、csh3 见 Comas-Vilà 对照集 |
| cuo、ce2o3、recl4、na2cro4 | SUSPECT | 换用黑铜矿（COD 9014580）、A 型 Ce2O3（mp-2721）、β-ReCl4（COD 4123951）、Cmcm Na2CrO4（COD 2107101）重算，cases/<体系>/v31 | 原运行与新运行都通过 |
| namno4、n2o、n4s4 | SUSPECT | 换用 AgMnO4 型 NaMnO4（mp-765753，与实验 AgMnO4 COD 2105345 同型）、按实验 Pa-3 胞搭的 P2₁3 有序 N2O、实验 S4N4（COD 9008994）重算，cases/<体系>/v31 | 原运行与新运行都通过 |
| cs2pt | SUSPECT | 改用实验晶格（Karpov 2003，a = 5.676、c = 9.471 Å）重算，cases/cs2pt/v31；原胞即 mp-13548 | 两者都是 Cs +1、Pt −2 |
| li2o2、re2（几何复核） | WRONG | li2o2 改用 Föppl 晶体 cases/li2o2/v3；re2 按本意的 2.40 Å 二聚体重算，cases/re2/v31 | 原运行与新运行都通过 |
| 其余 SUSPECT | SUSPECT | 保留，数据来源中注明 | |

另外，Comas-Vilà 对照集中的 MgXe 和 MH4 氢化物（mgh4、srh4、ceh4、thh4）原来
使用 0 GPa 的近似结构，核查时判为 OK，之后换用了文献中的精确结构（Bi & Zurek
2021 SI；MgXe 按 Miao 2015 推定晶格）并重算，数据在 cases/<体系>/v3。

所有替换结构的来源、空间群和自检记录在 cases/<体系>/structure-fixed/SOURCE.md，
汇总表为 results/structures_fixed_manifest.tsv。

**金标表（74 案例原运行）中仍有结构错误的案例：** CO、sno、ti2o3、li2o2、re2，以及
存疑的 cuo、ce2o3、recl4、na2cro4、namno4、n2o、n4s4。它们在正确结构上都已重算并通过，见 RESULTS.md
「结构修正版」和 results/bench74_corrected.tsv；引用基准成绩时可以两表并列。

## 几何复核（2026-10-21）

按 VASP 实际读入的几何（POSCAR 用笛卡尔坐标时坐标也要乘缩放因子）扫描全部
结构的最近原子间距，查出三处坏几何：

| 结构 | 问题 | 依据 | 影响 |
|---|---|---|---|
| li2o2/v1 | 缩放因子 10 作用到笛卡尔坐标上：两个 Li 重合于 (0, 5, 0)，两个 O 相距 5 Å，不是过氧根 | cases/li2o2/v1/OUTCAR 的离子坐标 | 74 案例中判为通过（v1_oxide） |
| re2/v1 | 缩放因子 12 作用到笛卡尔坐标上：两个 Re 相距 4.8 Å，不是 2.4 Å 的二聚体 | cases/re2/v1/OUTCAR | 74 案例中判为通过（v1_oxide） |
| li16o16、li15o16、li14o16 的 v3，以及 li15o16/v3-md | 2×2×1 超胞的第二组 Li 放在过氧根所在列上，落在 O–O 键中点（Li–O 0.78 Å） | POSCAR；单点计算 NSW = 0 | RESULTS.md「Li2O2 脱锂」的三行表格用的是这组几何；v3-md 是在第一步后崩溃的 MD 尝试 |

原核查按未缩放的坐标读 li2o2、re2，得出的键长不成立，两行已在登记表中更正。三处都已在正确几何上重算（见上文处理状态和 RESULTS.md）。
li15o16 的 v3-u、v3-md2 和极化子系列来自从坏几何出发、做了 31 步弛豫后的结构
（SAI: LEOS/inputs/m5a/li15o16_relax），与 Föppl 型 Li2O2 相符，几何正确。

历史批次（v1）以手搭结构为主；新批次（Comas-Vilà 对照集）以文献和数据库结构
为主，质量明显更好。
