# LEOS 方法学论文 — 工程实现计划

（2026-07-07 与 Claude 逐项确认后的决策记录与任务分解）

## 1. 项目定位

LEOS（Lattice Effective Oxidation States）：VASP 6.4.3 补丁 + 独立分析包，从
Wannier 中心/投影自动判定固体与分子晶体的 IUPAC 氧化态。目标：一篇方法学
期刊论文（JCTC / npj Comput. Mater. 量级），年底前投稿。

**核心主张（A）**：把氧化态的电子归属表述为带容量约束的最优运输 / 最小费用流
问题——统一几何、电负性（或第一性投影权重）与物理容量约束，全局最优、可证；
熵正则化（Sinkhorn）软分配给出共价度连续描述符。
**兼顾（B）**：稳健自动化工作流（诊断-自愈闭环、WF 质量检查、全量基准）。
**兼顾（C）**：去经验化——PAW 投影权重替代电负性表（作为 A 框架内的对照实验）。

## 2. 关键文献坐标（新颖性边界）

- **Comas-Vilà et al., JCTC 2025, 21, 7075**（工作区有 PDF）：EOS 名字出处。
  QTAIM+AOM+EFO Aufbau，QE+critic2+APOST-3D 重后处理链。40 固体测试集 =
  我们的主基准（参考值直接引其 SI，不自建链路复现）。我们的差异化：
  免积分、单程序原位、轨迹能力；他们有 R(%) 可靠性指标、fragment、简并平分，
  我们必须给出对应物。
- **Gimferrer, Comas-Vilà & Salvador, Molecules 2020, 25, 234**：已系统论证
  "定域轨道中心 + 最近原子判据"的失效（O-H/N-H 错判氢化物、依赖定域化方案）。
  → 不能声称首次发现失效；贡献定位为"给出修复失效的全局最优框架"。
  其失效案例类型纳入验证（数据集已有 SrBe(OH)4、BH6N）。
- **Marzari 组, npj Comput. Mater. 2025（arXiv:2412.01652）**：DFT+U+V 占据矩阵
  的 MD 氧化态追踪 + redox-aware MLIP。→ "首次轨迹 OS"不能写；差异化打点：
  全元素、无 U/V 参数、能标氧价态（阴离子氧化还原是对方结构性盲区）。
- Jiang-Levchenko-Rappe PRL 2012（Wannier 中心计数开山）、2013 Bader/MLWF/投影
  对比（ResearchGate 258751163）：需引用。
- 运输/最小费用流做电子归属：两轮检索无先例，为本工作新组件。

## 3. 已确认决策（grill 记录）

| # | 决策点 | 结论 |
|---|--------|------|
| 1 | 交付形态 | 期刊论文（学位论文模板已删）；9 月底方法定型、年底投稿 |
| 2 | 核心主张 | A 运输形式化为主，B 稳健工具、C 去经验化兼顾 |
| 3 | 软件形态 | **双路径 C**：Fortran 原位快速路径（贪心修复版，服务 MD）+ Python 包（运输、诊断、基准）；算法在 Python 定型后再评估回移 Fortran |
| 4 | 代价函数 | **档 2**：α·d² − β·p(w,a)，p = WF 在原子 PAW 球投影权重（CPROJ×U 旋转导出）；电负性表降级为回退项 |
| 5 | 容量约束 | 默认**宽容量**（价壳物理界：主族 8 / TM 18，按 PAW 价态配置生成）；紧容量（已知 OS 库）做可选诊断模式，两模式分歧输出警告（= 新奇价态探测器） |
| 6 | 稳健性分层 | VASP 内只诊断（(k,band) 级电子丢失、per-WF 全量导出）；自动补救在 Python 工作流层（参数矩阵重跑）；不上 site_symmetry |
| 7 | 效率取舍 | SCDM-only 基准 + MD warm start 进本轮；Γ-only、库内存直传 → future work |
| 8 | 应用演示 | **A+B+C 全做**，按风险排序：A Li2O2 氧化分解 AIMD（阴离子氧化还原）→ C CeO2 氧空位 Ce3+/4+（与 f 电子稳健性专项协同）→ B 富锂 Li2MnO3（MACE 微调+抽帧，设 2 周时间盒，超期降级为展望） |

## 4. 架构

### Fortran 侧（patch v2，目标文件 patch/leos.patch）
- 修复：电负性表换 Allen（可选 Pauling）；`MOD(V+60,1)` → `V-FLOOR(V)`；
  删除匈牙利分支（容量模型结构性错误：ATOM_CAP 按中性原子设槽，阴离子
  永远拿不到超中性电子）；main.F 缩进污染清理。
- 诊断：占据带落 lwindow 外的电子丢失，定位到 (k点, 带) 写入输出；
  EOS_NELECT_ABORT 语义待数据支撑。
- 导出 EOSCAR v2：per-WF 中心/spread/占据/近邻原子表/**PAW 投影权重 p(w,a)**
  （CPROJ 用 U_matrix 旋转到 WF 基，LORBIT 类基础设施）。
- warm start：EOS_EVERY_N 轨迹模式复用上帧 U 矩阵；EOS_SCDM_ONLY 快速模式开关。
- 自动窗口：GAP_TOL=1e-3 修正；带交叉时按能量直方图找 gap。

### Python 包（新建，名称待定，暂定 pyleos）
- 解析器：EOSCAR v2 / wannier90.wout / POSCAR / EIGENVAL。
- 运输分配器：宽容量 min-cost flow（整数解）+ Sinkhorn 软分配（共价度/置信度）；
  等价类容量池（简并混价，Fe+2.5 类）；fragment 汇节点。
  参考实现：script/ot_assign_proto.py（纯标准库 SSP，用于对拍）。
- 稳健性工作流：坏 WF 检测（spread 阈值+中心合理性）→ INCAR 参数矩阵自动重跑；
  失败三级归因（DFT / Wannier / LEOS）。
- 基准工具链：74+40 案例自动对比 Bader/文献值/四配置（贪心、档1、档2、SCDM-only），
  参数敏感性扫描报告。

## 5. 里程碑

- **M1 运输分配器原型**（Python，现有数据即可开发）：包骨架、档 1 代价、
  等价类/fragment、45 分子集回归。验收：As4S4 容量失控修复（已验证 ✅）、
  与贪心全面对比。
  **进度 2026-07-07**：pyleos 包落地（parsers/capacity/transport/cli/tests，
  纯标准库；单元测试 7/7 通过；As4S4 冒烟与原型一致）。V1 未含：等价类
  容量池、fragment、紧容量模式、per-WF 占据（等 EOSCAR v2）。
  **M1 验收达成（同日）**：36 案例回归 33 OK（β=10、Allen），运输 vs 贪心
  6 优 1 劣（劣的 CSe2 是输出截断的坏数据）。调试中修复三个求解器 bug：
  ①SSP 势更新遗漏不可达节点→指数退化；②BIG_NEG=1e5 的 ULP 与 1e-12 容差
  同阶→Dijkstra 假改进挂死（诱因是 F 的 n_min=7 首次触发负边，CH2F2 定位）；
  ③满占据假定注入幽灵电子（li2o2 的 num_wann=9>7）→ 解法：以 EOSCAR 的
  NELECT(Wannier) 为总流量上限的 min-cost k-flow，幽灵 WF 自动排除，
  CO 的电子丢失也因此免疫。**方法学发现**：β 翻转点在 6-10（X-H 键几何
  偏移 vs 电负性競争），β∈[10,15] 平台稳定；金属 ReO3、过氧 li2o2 运输
  全对而贪心 dev 2.0/1.0。剩余 3 案例（As4S4 ~、CSe2 FAIL、feo 曾波动）
  全部归因 stray WF/坏数据 → 稳健性闭环的目标清单。
- **M2 Fortran patch v2**（与 M1 部分并行）：上述修复+导出+warm start。
  **进度 2026-07-07**：patch v2 已重写（909→746 行）：删匈牙利分支、
  Allen 电负性默认（EOS_ENEG_PAULING 回退）、FLOOR 回卷修复、GAP_TOL
  1e-3→5e-2、main.F 缩进清理、(spin,k,band) 级电子丢失诊断（OUTCAR+
  EOSCAR 头 N_LOST_*）、EOSCAR #WF 明细块（EOS_WF_DETAIL）。静态平衡
  检查通过，**未编译验证**（本机无 VASP 源码树，M3 在 sai-pod 上验证）。
  推迟到 M3：PAW 投影权重导出（CPROJ×U）、EOS_SCDM_ONLY、warm start、
  能量直方图找 gap。
- **M3 sai-pod 编译与批跑**：官方工具链（vasp/6.4.3-nvhpc25.7-ompi5.0.10 或
  自编 CPU 链）+ wannier90；GPU VASP ↔ CPU W90 库兼容性实测（退路：CPU 编译）；
  Python 批跑驱动替代 bash（16V100 分区 4-GPU 粒度，一节点 4 算例并行）。
  **进度 2026-07-07**：发现生产树 vasp.6.4.3_ccqn 已含 LEOS v1 + W90 链接
  （历史 198 算例即出自它）。新建隔离树 ~/apps/vasp.6.4.3_leos_v2（v2 eos.F
  覆盖 + GAP_TOL 修改 + makefile.include 去 ccqn 硬路径），
  `make DEPS=1 -j8` 并行编译通过（nvhpc 24.3 GPU 链，作业 593747）。
  CO2 验证（作业 593776，4×V100）：EOSCAR v2 全功能正常——N_LOST_* 诊断头、
  #WF 明细块（占据/中心/spread/归属）；物理结果 C+4/O−2 正确，dNe=9.6e-11；
  Allen 标度生效（8 个 C-O 键 WF 全判给 O）。pyleos 已支持 v2 解析并端到端
  对账通过（#WF 占据总和=NELECT，运输结果与贪心/理想三方一致）。
  测试算例位于 ~/leos_v2_test/CO2。**M3 剩余**：Python 批跑驱动、
  投影权重导出（档 2 数据）、EOS_SCDM_ONLY、warm start（M2 递延项，
  现在有编译环境可开发）。
- **M4 全量基准**：Comas-Vilà 40 固体输入准备（COD/CCDC 结构）+ 74 现有案例，
  四配置全跑 + 敏感性扫描。
- **M5 应用演示**：A → C → B（见决策 8）。

## 6. 可行性验证记录（2026-07-07，As4S4，档 1）

- 现状贪心（EOSCAR）：S 从 −8 到 0（均值 −4.12），As +3..+5 —— 灾难性。
- 运输 + 宽容量（纯距离）：S −2..+2（均值 −0.88），As −1..+3 —— 容量失控修复。
- 电负性项 β 扫描（1→10）：均值饱和于 ±1.00，收益有限。
- 根因诊断：88 个 WF 中 19 个 spread>2.5 Ų（最大 22.6，"游离"型 15 个），
  Wannier 化本身未收敛（num_iter=100 硬停）。
- **结论**：容量约束修复灾难 ✅；分配层救不了坏 WF 输入 → 稳健性闭环是
  运输框架的前置条件；档 2 投影权重对坏 WF 仍有信息（中心无意义时投影
  权重依然可用），是鲁棒兜底。CO 案例（dNe=0.354，带交叉部分占据）留待
  档 2 的 per-WF 占据导出后验证。

## 7. 风险清单

1. GPU VASP + W90 CPU 库兼容性（M3 实测，退路 CPU 编译）。
2. Ce 4f 窗口顽固（C 演示降级为限制讨论）。
3. Comas-Vilà 结构获取：CCDC 需订阅（学校代理），COD 免费，高压结构在原文 SI。
4. B 演示 MLFF 链路重（时间盒 2 周）。
5. 审稿人圈 = Salvador/Contreras-García 系，对比实验与引用必须做足。

## 8. 服务器要点（sai-pod）

- 16V100 分区：每作业 ≥4 GPU（QOS huge-gpu 30d / rush-gpu 1d）；禁用
  --cpus-per-task/--mem*（站点钩子自动杀）；mps_mapping source 行保留。
- 用户已有 wannier90+VASP 编译链：~/apps/vasp.6.4.3_ccqn/build_c0.sbatch
  （fftw/3.3.10 + nvhpc/24.3 + ompi/5.0.8），LEOS patch 可套用同链路，
  注意与 CCQN patch 的源码树隔离。
- MLFF 环境现成：mace0.3.15 / deepmd / gpumd 模块（B 演示用）。
