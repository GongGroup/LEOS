# LEOS 项目全面审查报告（2026-07-10）

> 审查对象：仓库 HEAD = `874d67d`（"v3.1 complete + transport solver ported in-situ"）+ 未跟踪的 `gpuw90/` + sai-pod 集群侧生产环境。
> 审查方式：只读审查。9 路并行深审（transport 求解器 / eos.F / Wannier 线 / 算法科学性 / 工程化 / 结果复算 / 文档对齐 / 集群核查 / CUDA 首测）+ 对每条 Blocker/Major 发现逐条对抗性核实 + 三次独立复算复现，共 35 个审查代理。
> 完整性声明：审查全程未修改仓库任何被跟踪文件（起止两次冻结校验：HEAD、`git diff` 哈希、未跟踪文件 sha256 清单均一致）；本报告是审查新增的唯一仓库文件。全部证据日志与复算产物归档在会话 scratchpad 的 `audit_evidence/`（见附录 A，注意该目录随会话清理，需要长期保留请自行拷贝）。

---

## 0. 执行摘要

**总体判定：核心科学主张站得住，工程与文档欠一轮"失同步"大扫除。**

- **零 Blocker**。求解器（Python 与 Fortran 两个实现）的正确性经数学论证 + 400 组随机对抗测试 + 27 个真实案例暴力枚举对拍全面通过；74 案例基准表、仲裁参数扫描表、legacy 六体系、Ce 三像逐原子账本、Comas 补齐集、M5-B 时间序列、gpuw90 金标——**全部关键入库数据都在本机逐字节复现**（其中 full_bench 与 arbiter_scan 各三次独立重算，MD5 三方一致）。gpuw90 还首次在真实 V100 上跑通（CUDA 与 CPU 达机器精度一致）。
- **15 条 Major，全部经红队对抗核实 CONFIRMED**；40 条 Minor（1 条降级 Info）、44 条 Info。Major 集中在四类：
  1. **论文数字失准**：主表"运输 61/73 (84%)"实为 **62/72 = 86%**（E-1/F-1）；贪心基线"约 42 (58%)"按仓库自身口径重算为 **49/72 = 68%** 且 58% 无出处（F-18）——论文主标题对比"84% vs 58%"两个数字都要重定；"β∈10-15 平台"被实测证伪，真实平台是 **[8,12]**（C-2）。
  2. **一个真代码缺陷（未编译验证，PLAUSIBLE）**：`EOS_PROJ_WEIGHTS` 在带并行（NB_PAR>1）下先取模方再跨 rank 求和，丢失跨带组交叉项——多 rank 生产运行的 **#WP 块（仲裁/档 2 的数据供给）数值可疑**（A2-1）。档 1 主结论不受影响，但论文引用仲裁 31/38 前需做单 rank vs 多 rank A/B 定谳。
  3. **批跑驱动 batch.py 三个实证缺陷**：FAIL_W90 重试是空操作（D-1）、v3.1 abort 诊断 EOSCAR 被 triage 判成 OK（D-2）、squeue 失效时静默空集可毒化账本并改写运行中作业的 INCAR（D-4）。
  4. **溯源与合规**：regress_v31 四件套中两件出自未存档的中间构建（S-1）；v3.1 二进制从未晋升生产路径，新作业仍跑 v3 物理（D-3）；200 个许可 POTCAR 仍被 GitHub 仓库跟踪，公开前是硬合规债（D-5）；20260710_index.md 与同提交内其他文档自相矛盾（F-4）。
- **版本裁决结论（第 2 节）**：既有 74 案例 / comas / m5b / legacy 数据**不需要**用 v3.1 全量重跑；建议只对"ISPIN=2 且 ghost>0"的约 20 个体系做定点重跑作为论文稳健性附录。

**论文前必办 Top-5**：
1. 重定基准主数字：运输 62/72=86%，贪心侧重跑入库一份评分表（本审重算 49/72=68%），同步 RESULTS §1/§3.1 与 PLAN 里程碑表。
2. A2-1 定谳：编译一次，跑单 rank vs 4 rank 的同体系 #WP 逐行 diff；坐实则修复后重产仲裁扫描的数据基础。
3. 用构建 612888（或按入库补丁新构建）重产 ce_fs_4gpu 与 mumnox_abort2，恢复 v3.1 验收证据的字节级可追溯。
4. β 平台句改写为 [8,12] 并把本审的 8 点 β 扫描表固化入 data/；仲裁"安全平台矩形"按逐点数据改写（矩形内含 unsafe 点 (1.5,0.10)）。
5. 公开仓库前 git filter-repo 清除 200 个 POTCAAR 的全部历史；同时把 m5b 180 帧 EOSCAR（集群单副本、约 MB 量级文本）做异地备份。

---

## 1. 审查范围、方法与限制

**覆盖六个用户点名维度 + 扩展**：代码正确性（transport.py 546 行、eos.F 1539 行逐语义）、Wannier 计算与优化（gpuw90、w90_offline、SCDM/warm-start 证据链）、分配算法定义与科学性（费用函数、β、仲裁参数、R 指标、诊断阈值、幽灵偏置、OS_BOUNDS）、工作负载分配与工程化（batch.py、sbatch 模板、测试/CI、数据治理）、全部计算结果可靠性（本机复算 + 集群对账）、文档对齐（6 份主文档 + 5 份报告共 60 条可检验声明逐条核）。

**方法**：发现 → 对抗核实 → 综合三段式。每条 Major 由独立红队以"推翻它"为目标重跑证据；full_bench 与 arbiter_scan 由三个互不知情的代理各自重算并三方 diff。

**限制（如实声明）**：
- A2-1（#WP 多 rank 归约）为静态代码分析结论，本机无 Fortran 编译器，未运行验证，判定 **PLAUSIBLE** 而非确证。
- Fortran↔Python 交叉对拍（A4.1）通过，但样本仅全仓库唯二带 OS(flow) 列的 EOSCAR（co2、vircef）；MD 轨迹级 flow 输出尚无入库样本（E-7）。
- Comas-Vilà "40 全集"映射需原文 SI（本地无法核，29 目录 + 6 scdm 已逐一核实）。
- gpuw90 生产尺度（nw≈288）GPU 加速比未测（金标 nw=4 不可外推，S7-3）；这是项目自身的开放项，非本审失职。

---

## 2. 版本与数据归属权威裁决

三方证据（本地 git 考古 + 集群 md5/sacct + 数据指纹）联判，结论如下（全文：audit_evidence/verify/version_adjudication.md）：

| 对象 | 版本 | 证据 |
|---|---|---|
| 本地 HEAD `874d67d` | **v3.1**（eos.F md5 `5a0fdcea`，1539 行；patch 内嵌段逐字节一致） | git 历史中 eos.F 仅两个内容版本：`887bb60`=v3、`874d67d`=v3.1 |
| sai-pod 生产树 `~/apps/vasp.6.4.3_leos_v2` | **v3**（`a579285b`）+ CCQN 叠加，07-07 09:36 构建 | v3.1 只存在于 `leos_merge/rebuild_v31`，**从未晋升生产，且无 vasp_gam** |
| 补丁闭环 | **两级全部独立复现成立**：拉 vanilla 源码 + 本地重打补丁，v3 与 v3.1 各 4 文件 md5 与集群逐字节一致 | 正面结论，patch/README 的闭环声明可信 |

**数据归属**：
- 74 案例中 v1 时代 36 例输入是 v1 代码产物（**论文需披露输入横跨 v1/v3 两代的异质性**）；m4 两波、comas_gap、m5b 180 帧、legacy 六目录全部出自 v3 生产二进制。
- v3.1 产物 = regress_v31 四件套 + vircef_off/EOSCAR。其中**仅 co2/vircef 出自与入库补丁一致的构建 612888**；ce_fs_4gpu（构建 612454）与 mumnox_abort2（构建 612747）出自源码未存档的中间构建，四次构建二进制尺寸互不相同实证源码有改动（S-1，Major）。patch/README 恰好漏记了唯一可追溯的 612888。
- w90_offline/li2mno3 两份 EOSCAR 实为 v3 库模式 frame_000 的拷贝（非离线产物，S-4）；vircef_off/EOSCAR 是中间构建输出的**本地孤本**（远端已被覆盖）。

**旧数据是否要用 v3.1 重跑？不需要全量重跑。** v3.1 的 +796 行中，abort 诊断、OS(flow) 列、内存优化对旧数据"只加不改"；唯一改数值口径的 per-spin num_wann 幽灵修复只影响 ISPIN=2 且自旋不对称的体系。暴露集可枚举：74 案例中 ghost>0 共 31 例 + vircef/Ce 三像；且 pyleos 电子数封顶已吸收幽灵总量效应，残余风险仅电子落位偏置（vircef Ce +4→+3 是已实证的判定级变化，方向正确、重跑件已入库）。**建议对 ISPIN=2 且 ghost>0 的约 20 个体系做一轮 v3.1 定点重跑作稳健性附录。**

---

## 3. 复算核验总表

| # | 声明 | 复算结果 | 判定 |
|---|---|---|---|
| 1 | reliability_74.tsv（74 案例含 R 指标） | 三次独立重算（E/A1/红队 rep3）与入库表逐字节一致，MD5 三方同 `acc29cd7…` | ✅ |
| 2 | arbiter_scan_20260710.tsv（6×5 网格） | 三次独立重算逐字节一致，MD5 同 `d853bfa9…`；13/30 SAFE、灾难点 (1.0,0.05)=37 氢错、默认点 (1.0,0.20) 三项全优复现 | ✅ |
| 3 | "运输 61/73 (84%)" | tsv 与重算均为 **OK=62**（wave1 10/12 非 9/12）；分母两个无参考案例口径自相矛盾；权威口径 **62/72=86.1%** | ❌ E-1/F-1 |
| 4 | "贪心 ~42/73 (58%)" | 同口径重算 **49/72=68%**；cli regress 自打印 "better in 6, worse in 1" 与 ~24 内部矛盾；58% 无出处 | ❌ F-18 |
| 5 | legacy 6/6（含 vircef herr=1） | 六目录逐原子 63/63 与 2026-07-09 存档全等；herr=1 可指认为原子 4 | ✅ |
| 6 | Ce 三像 C-H 活化叙事 | 价态 + 几何双通道复现：H60 唯一氢负离子（ce_ts ne=2、OS −1），C-H 1.10→2.35→2.91 Å、Co-H 1.60/1.66 Å；Co4 片段 0→+2→+2 | ✅ |
| 7 | Comas 39/40 | 29 目录 + 6 scdm 逐案例元素平均全部命中 §3.13 表；**但 summary.tsv 内 4 个 FAIL_DFT 标签陈旧**（按它算是 35/40） | ✅（E-3 陈旧文件） |
| 8 | regress_v31 三件 | co2 整数精确（dNe 7.1e-15）、vircef ghost 0/Ce+3、mumnox abort 指路报错——均复现；但 mumnox 归档缺 POSCAR、#NLD 落盘仅 33 条（计数字段 1348） | ✅（E-5 归档瑕疵） |
| 9 | M5-B 整数阶梯 vs NMC 涨落 | 9×20 帧本地在档；TM 序列恰 4.000/4.500/5.000、漂移 0.000；NMC Li/O 钉死、Ni 端点涨落吻合（"随脱锂加深增强"仅端点成立，E-6） | ✅ |
| 10 | gpuw90 金标 GaAs 8.31e-7 / Diamond 3.13e-9 Å | 本机 CPU、集群 CPU、集群 **CUDA（V100，首测，作业 613251）** 三处逐位复现；fd_check ratio 全 +1.0000；CPU↔CUDA 机器精度一致 | ✅ |
| 11 | Fortran OS(flow) vs pyleos | co2 12/12、vircef 74/74 原子零分歧（含 H6=0、C34=0 非平凡值）；常量表 102 元素 0 不一致；patch↔eos.F 逐字节同源 | ✅ |
| 12 | SCDM-only 轨迹失效 stray≈288、67× | 集群原始 EOSCAR 核数成立；唯"完整迭代侧 stray 0-1"实为 **0-4**（B-5） | ✅（数字微修） |
| 13 | OOM 六连、MaxRSS 56.2/106.7 GB、编译 594329≈11min | sacct 逐字吻合（Elapsed 16-21 s 精确命中） | ✅ |
| 14 | "集群队列已清空" | HEAD 提交时刻 scf-G0/G1/G2 在跑、scf-G3 与 w90-export 当日各两次 FAILED，scf_recipe 仓库零记载 | ❌（S-2，声明时点已过期） |
| 15 | 4 GPU 内存 "41 GB=66%" | sacct MaxRSS 49.2 GB=79%（SCDM 段峰值）；41 GB 是稳态口径，两处未互注（F-17） | ⚠️ |

---

## 4. 代码正确性

### 4.1 transport.py 最小费用流求解器 —— 全面通过

- **负费用边 + Dijkstra 势函数**：初始负弧清单完整（n_min 奖励 −1e3；linear γ>0 弧；arbiter 翻案下压），初始残量图为 DAG 无负环，Bellman-Ford bootstrap 必精确收敛；逐轮势更新（不可达节点抬 dmax）满足 Johnson/SSP reduced-cost 非负不变量——在最大三个案例（rbpb 736e/mocl2/srbe_oh_4）解后遍历全部残量弧，最负 reduced cost 1e-15 量级，对 1e-9 容差裕度 5 个数量级。
- **BIG_NEG=1e3 量级分离**：74 案例全局最大单弧费用 279.86（sega），分离 3.57 倍；ULP(1e3)=1.14e-13 ≪ 1e-9。注释里 1e5/1e-12 旧组合失效的故事定量成立。
- **整数解**：supply/capacity/target 全整数路径核实；400 组随机小网络 + 200 组全负费用网络与暴力枚举**零失配**；不可行/超容量/简并平台/max_rounds 边界全部行为正确。
- **reliability_margins**：用 DP 精确枚举整数流多面体，覆盖全部 na≤5 的 27 个真实案例约 120 个 margin 值**全部匹配**（典型偏差 1e-11）；delta 断言全量零触发。
- **仲裁模式**：实现与 RESULTS §3.8 / HANDOFF 任务三逐句一致；38 个 #WP 案例 tier1=29 OK → arbiter=31 OK（恰 hi3 与 Rh 配合物翻正、零反向）精确复现。
- 附加核验：最小镜像 ±1 邻胞对本数据集精确（±2 扫描零差异）。

值得记录的缺陷：**C-1（Major）签名默认 β=3.0**——这是文档记录在案会翻错 X-H 键的取值（β=3 时 74 案例只 55 OK，CH2F2 给出 H=−1/C=+4 且静默无告警）；当前所有生产调用点都靠显式传 10 规避，但任何新代码 `assign_case(case)` 裸调即中招。另有：cli 无 `--tier2` 入口（A1-2）、arbiter 分支静默忽略 gamma（A1-5）、正确性依赖的两处 assert 在 `python -O` 下失效（A1-6）、银行家舍入与 Fortran NINT 的 exact-.5 口径差（A1-3/A2-5，当前数据零触发）。

### 4.2 eos.F —— 移植质量高，两个实质问题

**通过项**：Fortran MCF 与 Python 逐语义等价（势/BF/dmax/容差/BIG_NEG/奖励弧/target/平局键/保险丝全一致）；自研堆逐分支推演无误；EOS_ANALYZE 全部 8 条退出路径 × 已分配数组矩阵无泄漏/重复释放；EOSCAR 仅 rank0 写、abort 分支全员同向无死锁；EOS_INTEGERIZE 守恒代数证明 + 手推三例；abort 协议与 parsers.py 两端实测闭环；电负性/容量常量表 102 元素与 Python 0 不一致；patch 内嵌 eos.F 与独立副本逐字节相同。

**问题**：
- **A2-1（Major，PLAUSIBLE）**：`EOS_PROJ_WEIGHTS`（eos.F:986-1040）在带并行下对局部带部分和先取 |·|² 再 M_sum_d，数学上得 Σ_带组|Σ_组内带|²，丢跨带组交叉项；正确做法是 CW 先在 COMM_INTER 归约再取模方。GPU 版 NCORE=1 强制 → 生产多 GPU 运行几乎必然 NB_PAR>1，Ce 三像 EOSCAR 的 #WP 正是 arbiter_scan 的 legacy 腿输入。**影响面**：仲裁/档 2 的数据基础（(1.0,0.20) 标定）；OS(flow)、贪心列、NELECT 诊断不受影响。**需编译 + 单/多 rank A/B 定谳。**
- **A2-2（Major）**：敌对评审报告的归档注记声称"4 个 MAJOR 已同日修复"，但 **M1 实际未修**（EOS_PREPARE_MLWF 仍只在 NSTEP==1 调用，per-spin 占据单调抬升在轨迹步≥2 是死代码——某自旋占据数增长时轨迹功能不可用，现在的缓解是响亮 abort 而非无声），M3 仅部分修（>255 列 .win 行仍静默截断，A2-3）。RESULTS §7 与 HANDOFF 的口径反而是对的——**是报告注记错，不是代码谎报**，但会误导后续维护者。

### 4.3 解析层与 Fortran↔Python 交叉

parsers.py 的 abort 帧识别、#NLD、flow_os 列、分块解析（评审 MAJOR-2 修复）实测闭环。交叉对拍：全仓库唯二带 OS(flow) 的 EOSCAR 逐原子零分歧；费用打印位内一致。注意 A2-6：离线复算吃 EOSCAR 6 位小数打印值，简并几何上平局理论上可翻转——对拍偶发单原子差异应先归因于此。

---

## 5. Wannier 线：gpuw90 与优化路径

**正确性交付真实可信**：金标、梯度自检、数值静审（MV 梯度前因子/Im ln 主支/反厄米指数/Löwdin/.mmn 列主序/b 向量 B1 完备性）全部与 wannier90 3.1 约定一致；testdata 与集群上游 test-suite 逐字节同源；性能估算复核无数量级错误（每迭代 30.6 Gflop、0.32 s/iter≈95 GF/s、GPU 现实 5-15× 的算术全部成立）。**本审补上了 CUDA 首测**（S7）：V100 上两金标 PASS、CPU↔CUDA 机器精度一致（max|Δcentre|≤1.1e-15 Å）、complex128 eigh 正常——DESIGN §9-3 前半句缺口关闭。

**核心缺口（B-1，Major）**：**全部生产与基准体系都需要 disentanglement**（Li2MnO3 356/288、VirCef 180/142、Ce 332/~250、v1 数据集全部 nb>nw），而 Phase 1 只支持 nb==nw——**生产覆盖率为 0**，且 Amdahl 账显示只 GPU 化 wannierise 的整帧收益上限约 2.4×，报告的 5-15× 依赖尚未动工的 Phase 2。DESIGN 披露了阶段边界但无一处明说"生产覆盖为零"。配套隐患：`wannierise()`/`load_case` 无 nb==nw 守卫，喂纠缠数据会**静默**产出貌似合理的非 MLWF 结果（B-8）——下一位使用者最可能犯的就是这个错。

其余值得注意：
- **B-4**：li2mno3 离线 A/B 的 spin2 解纠缠**未收敛**（打满 200 迭代），103/288 WF 偏差>0.1 Å、2/96 Mn 位点级翻转；判定级（元素平均）等价成立，但文档"收敛通道逐 WF 同解"的措辞掩盖了另一通道未收敛。离线工作流投产用于位点级分析（极化子）前必须调高 dis 迭代上限重跑。
- **B-2/F-11**：warm start 一词两义——RESULTS 实测无增益的是 guiding-centres 注入（+3.2%），可行性报告标 3-10× 的是帧间 U 复用（未实测）；两处互不引用，读者极易误以为 3-10× 有实测支撑。
- **B-6**：DESIGN 显存账漏计每迭代 ~8 个同尺寸 torch 中间量，真实峰值 ~3-4 GB（差一个量级）；16GB 单帧结论不翻，但"单卡放几十帧"的直觉要修正为 3-4 帧。
- **S7-4**：热循环每迭代至少 2 次 device→host 同步（float() 取回），生产尺度会封顶 GPU 利用率——Phase 2 落地时应改为设备端收敛判断。
- **SCDM-only 证据链**（B-5）：stray≈288、67×、0/12 全部在集群原始数据核数成立；"完整迭代侧 0-1"实为 0-4。escalate 闸门缺 STRAY_DIST=1.6 距离判据且未读 v3.1 FLOW_N_STRAY 头（B-3）——"中心离位但 spread 紧凑"的坏 WF 会漏检。

**优化路线评估**：项目自己的判断（先离线 MPI k 并行 + 自旋并发的零代码路线，GPU 化作为独立后处理工具、不碰上游）经复核是正确的。真正的关键路径是 **Phase 2 disentangle**——没有它，GPU 路线对生产负载是空转；有了它才谈得上 5-15×。w90_offline 的 li2mno3/vircef .mmn/.amn/.eig 就是现成的纠缠金标。

---

## 6. 算法科学性

**站得住的**：β=10 在实测平台 [8,12] 腹地（β=10↔12 零翻转）；仲裁默认 (1.0,0.20) 在 SAFE 平台腹地、四邻全 SAFE；γ=0 时 tier2=linear 严格无副作用；OS_BOUNDS 15+ 元素与文献相符、DEFAULT_BOUNDS 与电负性兜底在 74 案例集**零触发**；幽灵偏置机理有干净的 A/B 实证（vircef 共享 num_wann ghost=2→Ce+4 vs per-spin ghost=0→Ce+3）；R 指标的实现（±1 完整重解）与 §3.9 定义一致且被暴力枚举验证。

**要修的**：
- **C-2（Major）**："β∈10-15 平台稳定无副作用"被证伪——β=15 翻掉 Cr(CO)6（dev 6.00）、β=20 再翻 tmao；该句是 36 案例时代结论，74 案例后未重验。实际平台 [8,12]（β=8 还多过一例 rhcl）。**C-3**：β 选型扫描无原始数据/脚本存档，本审的 8 点重扫表（audit_evidence/findings/C_bench_b*.tsv）应固化入 data/。
- **C-1（Major）**：库默认 β=3.0（见 4.1）。
- **仲裁统计正当性**（C-5/C-6/E-2）：安全平台矩形内含 unsafe 点 (1.5,0.10)（入库 tsv 自身数据），矩形表述应改为"TAU≤1.0 × PM∈[0.10,0.30] 全 SAFE；TAU=1.5 需 PM≥0.2"；选参与验参共用同一 44 案例、无留出集——首版 (1.0,0.05) 事故正是样本内选参会误导的现成实证，文档应明写此限制；SAFE 判据不含非氢翻动（(0.25,0.02) 带 15 个非氢翻动仍 SAFE，C-16）；"31/37 vs 31/38"分母混用（F-6）。
- **R 指标文档**（C-8/C-9）：§5.4 的白话解释描述的是 wf_margin（单弧下界）而实现是完整重解——同一文档 §3.9 是对的、§5.4 掉队；且 R 作为失效探测器的假阴性未披露（10 个偏差案例中 5 个 R_elem≥54，wcl6 高达 83——R 量"解的唯一性"不量"正确性"，必须配 stray/ghost 使用，否则审稿人会拿 wcl6 反驳）。
- **诊断阈值溯源**（C-10）：STRAY_SPREAD=3.0、ABORT=0.35、batch dNe 1e-3 均为无据经验值（对照 STRAY_DIST=1.6 有书面推导）；分布检验显示 dNe 类不临界，但 stray 计数有 3 案例随 ±10% 扰动翻"干净/带警告"。
- **归因修正**（E-4/C-15）：vircef 残留 H(0) 的"幽灵旧账"归因被 v3.1 自身验收数据**证伪**（ghost 归零后 H(0) 仍在，stray 反而恶化到 36.1 Å²）——根因指向 Wannier 化质量/stray，RESULTS §3.8、arbiter_scan.py docstring、index.md 三处表述需更新。
- **费用函数论证缺口**（C-4）：d²（而非 d）的选择、β 项量纲（β=10 隐含"电负性满差 ≈ 3.16 Å"）、χmax 按案例内归一（同一元素对在不同体系费用差近 2 倍、且与 ARBITER_TAU 绝对阈值耦合）——三处在全部文档中无一句书面论证，是审稿人最可能追问的方法学问题。
- 小项：Ir 上界 8 排除已报道的 Ir(IX)（当前案例集无 Ir，休眠）；eneg() 对缺表 32 元素（Tc、Pr-Lu、Tl、锕系等）静默返 0.0 → 吃满 β 罚且无警告（A2-4）——对当前基准零影响，对论文扩展集是暗雷。

---

## 7. 工程化与可复现性

**正面**：full_bench 确定性双跑逐字节一致，且与入库基线表可复现——**天然 CI 金标**；两份 sbatch 模板对照 sai-pod 站点规范全项合规（QOS/≥4 GPU/无 --cpus\*/--mem\*/mps_mapping 守卫）；pyleos 测试 12/12 过；HANDOFF_V31 交付清单 7 项全部落地勾对；分析侧新人可零依赖走通（cli case / regress 实测）。

**Major 缺陷（全部实证复现）**：
- **D-1**：`RETRY_LADDER["FAIL_W90"]` 唯一键 NBANDS_SCALE 被 retry() 过滤——INCAR 逐字节不变、计数+1、原样重投，随后梯子耗尽；docstring 声称的 EOS_SCDM_ONLY=.FALSE. 也未实现；重试日志还打印 "staged retry #1 (NBANDS_SCALE)" 主动误导。每个 FAIL_W90 白烧 4 GPU × walltime 一次。
- **D-2**：triage() 只查 EOSCAR 存在即 OK——v3.1 abort 诊断 EOSCAR（真实样本 |dNe|=5.49）被判 OK 且永不重试，FAIL_LEOS 分类与其重试梯成死代码；summary.tsv 会虚报成功。此为 v3.1 引入的真回归（v3.0 abort 不写 EOSCAR）。当前生产二进制还是 v3.0 所以尚未生效——一旦晋升 v3.1 立即暴露。
- **D-3**：VASP_BIN 硬编码个人路径且指向 v3 二进制；v3.1 构建未晋升（无 vasp_gam）；"leos_v2 目录装 v3 代码"的命名陷阱。现在发的任何新作业都跑旧物理。
- **D-4**：squeue_ids() 不查 returncode、USER 空时静默空集（集群实测 `squeue -u ""` 返回 0 行 rc=0）→ 运行中作业被提前 triage 成 FAIL、INCAR 被改写、可重复投递；QUEUED/RUNNING 分支不回写状态，一次瞬断永久毒化账本。
- **D-5**：200 个许可 POTCAR（25 个唯一 blob）在 GitHub 私有仓库被跟踪；公开前必须历史重写。

**Minor 面**（选列）：集群 batch.py 副本落后仓库 5 行（bvector→FAIL_W90 修复未部署，D-6）；harvest keep 列表含死条目 vasp.out、缺 KPOINTS/CONTCAR/*.win（D-7）；triage 漏斗把 walltime/OOM 全判 FAIL_DFT 吃错药、DONE/RETRY_DONE/ESCALATED 死状态（D-8）；retry/triage/submit/harvest/parsers/cli 零测试、无 CI（D-9，测试空白与 D-1/D-2/D-4 缺陷位置精确重合）；`git add gpuw90/` 会被 .gitignore 的 \*.amn/\*.mmn 静默拦下 4 个金标输入（D-10）；无 requirements（D-11）；full_bench/arbiter_scan 的 --out 默认值直指已提交数据（本审已两次踩到设防，D-12）；根 README 的 patch 导航指向 v2 旧补丁（D-13）；reference_os.json 40/73 案例无逐条出处注（D-14）。

---

## 8. 文档对齐矩阵（60 条声明核验汇总）

绝大多数硬数字当场命中（R 中位数、Li2O2/CeO2/M5-B 全部数值、Ce 三像 dNe、OOM 取证、gpuw90 金标、仲裁扫描 30 格逐字节）。失配清单：

| 声明位置 | 文档说 | 实际 | 级别 |
|---|---|---|---|
| RESULTS §1/§3.1、PLAN:182 | 运输 61/73 (84%)、wave1 9/12 | **62 OK、wave1 10/12；权威口径 62/72=86.1%** | Major E-1/F-1 |
| RESULTS §3.1 | 贪心 ~42 (58%) | 同口径重算 **49/72=68%**；58% 无出处 | Major F-18 |
| RESULTS:458 | β∈10-15 平台无副作用 | β=15 翻 Cr(CO)6；平台实为 [8,12] | Major C-2 |
| index.md 多处 | fortran_review 待补 / 五项未编译 / 4GPU 未测 / patch 1349 行 | 同提交内均已完成 / 2112 行——**中间时刻快照未更新** | Major F-4 |
| RESULTS:253-255 | 安全平台矩形零新增氢错 | 矩形内 (1.5,0.10) unsafe（自家 tsv） | Minor C-5/E-2/F-5 |
| index.md ×3 | 31 网格点 | 30（6×5）+1 基线 | Minor F-3 |
| RESULTS §3.8 | 基准 31/37 与 31/38 混用 | 38 含 1 个 run-only，需注明 | Minor F-6 |
| comas_gap/summary.tsv | 4 个 FAIL_DFT | 四目录健康 EOSCAR、重算全过——陈旧快照 | Minor E-3/F-7 |
| RESULTS:559/PLAN:255 | 队列已清空 | scf_recipe 三作业在跑、G3/w90-export 双 FAILED、零记载 | Info S-2/F-16 |
| patch/README:47 vs RESULTS §7 | 41 GB=66% vs 49.2 GB=79% | sacct 实证 79% 为峰值口径，两处未互注 | Minor F-17 |
| RESULTS:195 | 完整迭代 stray 0-1 | 实为 0-4 | Info B-5 |
| RESULTS:206 | warm start +3.5% | 实算 +3.2% | Info F-10 |
| RESULTS:404 | legacy 七目录在 data/ | 本地六目录（缺 reo3 镜像） | Info F-15 |
| patch/README | v3.1 六项/编译三轮 | 漏记运输回移（OS(flow) 列/EOS_FLOW 开关）与第四次构建 612888 | Minor F-8 / Major S-1 |
| 三份报告 | eos.F 行号 | 基于 787 行 v3 版，HEAD 1539 行全漂移，未注明基线 | Info F-12 |

---

## 9. 集群侧纪要（sai-pod，全程只读 + 一个轻量验证作业）

- **S1 版本裁决**：见第 2 节。
- **S2 M5-B**：180 EOSCAR 精确对账（9×20）；本地镜像抽查 md5 全等；**m5b_dft 122G + scdm 15G 为集群单副本**——论文主图底层数据无异地备份，重算需 ~50 GPU 日（S-8：建议至少把 180 个 EOSCAR+INCAR+POSCAR 文本级备份）。
- **S3 补丁闭环**：独立复现成立（正面）。
- **S4 作业史**：sacct 与三份报告叙事逐字吻合。
- **S5/S6 环境与脚本**：模块/库/二进制在位；集群副本无私改未回流（仅滞后，parsers=c5ab54d 版、transport=90999a4 版——在集群直接复跑分析会用无 flow_os 的旧解析器）。
- **S7 CUDA 首测**：16V100 节点作业 613251（4 卡分配 1 卡计算、flood-gpu、22 秒）——两金标 PASS、CPU↔CUDA 机器精度一致、fd ratio +1.0000；金标尺度 GPU 无加速（nw=4 不可外推）。**8A100M40 唯一节点宕机（unk\*）**，1-2 GPU 快速验证通道当前不可用。远端审计副本在 `sai-pod:~/tmp_gpuw90_audit/`（复查后可删）。
- 07-10 家目录出现新散件（首轮 mpirun rc=213 失败日志等），归拢纪律回退（S-7）。

---

## 10. 发现总登记表

**Major 15 条（全部红队 CONFIRMED）**：

| ID | 一句话 | 位置 |
|---|---|---|
| A2-1 | #WP 多 rank 归约丢跨带组交叉项（PLAUSIBLE，需编译定谳） | eos.F:986-1040 |
| A2-2 | 评审 M1 未修但报告注记称已修；轨迹占据增长场景不可用 | eos.F:277-286 / fortran_review.md |
| B-1 | gpuw90 生产覆盖率 0（全部生产体系需 disentangle）；只 GPU 化 wannierise 上限 ~2.4× | gpuw90/ + w90_offline |
| C-1 | 库默认 β=3.0 是已知会翻错 X-H 键的取值 | transport.py:365,445 |
| C-2 | "β∈10-15 平台"被证伪，实为 [8,12] | RESULTS:458 |
| D-1 | FAIL_W90 重试空操作 + docstring 虚假承诺 | batch.py:292,315 |
| D-2 | abort 诊断 EOSCAR 被 triage 判 OK（v3.1 真回归） | batch.py:253 |
| D-3 | VASP_BIN 硬编码且指向旧版；v3.1 未晋升、无 gam | batch.py:40 |
| D-4 | squeue 失效静默空集 → 毒化账本/改写运行中作业 | batch.py:221-227 |
| D-5 | 200 个许可 POTCAR 被 git 跟踪（公开前合规债） | data/v1/ |
| E-1/F-1 | 主表 61/73 应为 62/72=86% | RESULTS §1/§3.1 |
| F-18 | 贪心 58% 无出处，同口径实为 68% | RESULTS §3.1 |
| F-4 | index.md 与同提交文档自相矛盾（快照未更新） | 20260710_index.md |
| S-1 | regress_v31 两件出自未存档中间构建；漏记构建 612888 | patch/README:41 |

Minor 39 条、Info 45 条按维度分布：A1 7 / A2 9 / B 12 / C 15 / D 18 / E 7 / F 18 / S 8 / S7 5（共 99 条，全文见 audit_evidence/digest_findings.md 与各 findings_\*.md）。

---

## 11. 改进建议（分级，仅记录不实施）

### A. 论文前必办（正确性/可信度直接相关）
1. **重定主数字**：62/72=86%（并决定 As4S4 "~" 的口径归属）；贪心侧用同一发现器重跑评分并入库（本审 49/72=68% 的脚本在 audit_evidence/verify/F18_greedy_recount.py 可直接复用）；同步 RESULTS §1/§3.1、PLAN 里程碑表。
2. **A2-1 定谳**：单 rank vs 4 rank 同体系 #WP 逐行 diff；坐实则修归约（CW 在 COMM_INTER 上 M_sum_z 后再取模方）并重产仲裁扫描数据基础；论文引用 31/38 前此项是闸门。
3. **溯源修复**：用构建 612888 重产 ce_fs_4gpu 与 mumnox_abort2（并给 mumnox 归档补 POSCAR）；patch/README 补记 612888 与"数据↔作业↔构建"对应表 + 运输回移条目。
4. **定点稳健性重跑**：ISPIN=2 且 ghost>0 的 ~20 个体系 + Ce 三像在 v3.1 下 A/B，作为附录；同时验证 Ce 三像幽灵是否归零、CH4 活化叙事逐原子保持。
5. **合规与备份**：公开前 git filter-repo 清 POTCAR 历史（以 SHA256+版本标识替代）；m5b 180 帧 EOSCAR 文本级异地备份。
6. **文档措辞修正包**（一次提交清完）：β 平台 [8,12] + β 扫描表入库；安全平台矩形逐点化；R §5.4 重写对齐实现 + 补假阴性口径；vircef H(0) 归因改 stray（三处）；61/37/38、31 网格点、stray 0-4、+3.2%、49.2 GB、index.md 快照注记、comas summary.tsv 重生成、"队列已清空"限定时点并补记 scf_recipe。
7. **方法节补论证**：d² 选择、β 量纲换算、χmax 案例内归一的取舍（或用全局固定 χmax 重扫一遍验证不变性）；披露 74 案例输入横跨 v1/v3；仲裁参数"选参即验参、无留出集"限制明写。

### B. 恢复生产前必办（工程）
1. 晋升 v3.1 二进制到 ~/apps 带版本号目录 + current 符号链接（补建 vasp_gam）；VASP_BIN 参数化。
2. 修 batch.py 三缺陷：实现 NBANDS_SCALE（或按 docstring 落实 EOS_SCDM_ONLY 降级）；triage 首查 `# LEOS_ABORTED`；squeue_ids 用 getpass 兜底 + returncode 检查 + 失败时保守哨兵；顺手清 triage 漏斗（sacct 终态感知，资源类失败单列）与死状态。
3. 同步集群 leos_analysis/tools 副本到 HEAD 并在 README.txt 标注提交号。
4. `git add -f gpuw90/testdata`（或 .gitignore 加 `!gpuw90/testdata/**`）后把 gpuw90 入库；escalate 闸门补 STRAY_DIST 判据或读 FLOW_N_STRAY 头。
5. 最小 CI：GitHub Actions 单 job 跑 unittest + full_bench --no-reliability 与入库表四列 diff（确定性已验证，直接可用作金标）；为 retry/triage/abort-EOSCAR 补三个单元测试（夹具都是现成的仓库文件）。
6. β 默认 3.0→10.0；eneg() 缺表元素打警告；--out 默认改带日期文件名；requirements-cluster.txt / requirements-gpuw90.txt。

### C. 可延后（性能与路线）
1. **gpuw90 Phase 2（disentangle）是 GPU 路线的先决里程碑**——用 li2mno3_off/vircef_off 的 .mmn/.amn/.eig 做纠缠金标；落地前先给 wannierise()/load_case 加 nb==nw 守卫（一行 NotImplementedError）。
2. 生产尺度 GPU bench（nw≈288 真实帧）回答加速比问题；顺手修热循环 device→host 同步（设备端收敛判断）。
3. 离线 w90 工作流投产前调高 spin2 dis 迭代上限重跑（位点级分析前提）；w90_offline 目录补 README 说明 EOSCAR 出处（B-11/S-4）。
4. warm start 术语拆分（guiding-centres 已证伪 vs U 复用未实测）并互引。
5. harvest keep 列表修正；reference_os.json 补 h2s3/h2s5 注 + 统一下划线语义；README 导航更新（m5b、patch 入口）；报告行号基线注记。

---

## 附录 A：证据归档

- 会话 scratchpad：`.../scratchpad/audit_evidence/`
  - `findings/`（77 个文件）：九路维度的完整证据日志（A1_transport.md、A2_eosF.md、B_wannier.md、C_science.md、D_engineering.md、E_recompute.md、F_docs.md、S_cluster.md、S7_cuda.md）、全部复算 TSV（reliability_74 三份、arbiter 三份、β 扫描 8 点）、对抗测试脚本（A1_adversarial.py、A1_reliability_brute2.py 等）、交叉对拍产物（a41_\*）、S7 作业日志。
  - `verify/`（49 个文件）：15 条 Major 的逐条红队复核报告、Minor 分维度批量复核、版本裁决全文（version_adjudication.md）、三方复算复现（replication.md）、贪心重算脚本（F18_greedy_recount.py）。
  - `digest_{dimensions,findings,verdicts}.md`：结构化汇总三件套。
- 远端：`sai-pod:~/tmp_gpuw90_audit/`（gpuw90 CUDA 首测现场，复查后可删）。
- 基线快照：`.../scratchpad/baseline/`（审查前后工作区哈希，证明仓库零污染）。

*报告完。审查执行：Claude（9 发现代理 + 24 核实代理 + 主会话综合），2026-07-10。*
