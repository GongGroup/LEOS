# cases：全部计算数据

一个化学体系一个目录，体系下每个子目录是一次计算（运行）：

```
cases/<体系>/<运行>/
```

每个运行在 [INDEX.tsv](INDEX.tsv) 中占一行。代码只通过这张表找数据
（pyleos/cases.py），不按目录名猜。

## INDEX.tsv 的列

| 列 | 含义 |
|---|---|
| case | 体系目录名：化学式小写，括号换成下划线，例如 cr_co_6、srbe_oh_4 |
| run | 运行目录名，见下表 |
| role | 用途：bench74（74 案例基准）、comas（Comas-Vilà 对照集）、test-regression / test-homonuclear（留出测试集）、clean-override（clean 集改用的重算）、app-li2o2 / app-ceo2 / app-fe3o4 / app-cathode / app-methane（应用）、legacy（旧程序失败体系的重算）、variant（同一体系的另一种设置）、diag（诊断与重算）、structure（只有结构）、v1-failed（没有产出 EOSCAR） |
| bench_set | 74 案例基准中的分组：v1_mol、v1_oxide、wave1、wave2；与金标表 batch 列一致 |
| ref_key | 该案例在 pyleos/sets/reference_os.json 中的键 |
| origin | 改成本布局之前的路径，旧文档、旧日志里的路径可按它对回 |

## 运行目录名

| 运行名 | 含义 |
|---|---|
| v1 | 项目早期 v1 程序（贪心）的计算；v1-mol 用于 ReO3，它在分子集和氧化物集各有一份 |
| v3 | 当前补丁（v3）的标准计算 |
| v3-recompute | 留出测试集的重算（gpuw90 从归档的重叠矩阵重新 Wannier 化） |
| v3-scdm-only | SCDM-only 快速模式，不计入对照 |
| v3-ispin2、v3-relaxed、v3-u、v3-md、v3-md2、v3-verwey | v3 的变体：自旋极化、弛豫晶格、+U、MD 快照、Verwey 相 |
| polaron-* | Li15O16 极化子测试（U 值、PBE/HSE 单点） |
| v31、v31-wp、v31-ghost、v31-mag*、v31-nowin、v31-abort2、v31-4gpu | v3.1 补丁的重算与诊断：带 #WP 投影块、ghost A/B、磁矩初值、不写 win、中止复现、4 GPU |
| w90-library、w90-offline | wannier90 库模式与离线模式的对照 |
| md-x0.00 等 | 富锂正极 MD：初始结构 POSCAR、生产日志 md.log、验收 verdict.json |
| structure、structure-hand、structure-fixed | 只有结构：原始结构与 Bader 分析、手搭结构、替换后的正确结构（含 SOURCE.md） |

## v1 数据的说明

v1 数据由项目早期的 v1 程序产生，分配方法是贪心，是 74 案例基准中"历史回归集"
（36 个：10 个分子晶体、26 个氧化物和卤化物）的来源。每个氧化物体系原有多次
运行，这里只保留基准实际使用的那一次；其余运行，以及从未进入基准的 zno、
na3vo4、li2o2、reo2 的 v1 运行，只在 git 标签 `pre-reorg` 中
（`git show pre-reorg:data/v1/oxides/...`）。results/v1_oxide_runs.tsv 是 v1
当年的汇总表，它的 workdir 列指向每个体系最新的一次运行，与基准使用的那次多数
不同，只作历史记录。

v1 汇总表中的"状态"只表示程序是否跑完、电子数是否守恒，不表示氧化态是否正确。
例如 As4S4 记为 OK，但其中 S 的氧化态低到 −8。Ca3AuN 的 v1 运行缺 EIGENVAL，
没有产出 EOSCAR。"co" 的胞内实际是两个乙烯二酮分子（O=C–C=O），不是一氧化碳，
见 docs/testset_findings.md。
