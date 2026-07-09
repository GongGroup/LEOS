# v3.1 技术负债交接（2026-07-09）

本文写给接手 v3.1 的 agent，目标是零上下文可开工：先读"背景五分钟"，
再逐条做三个任务。每条给出现状、根因证据、代码位置、改法方向、复现
体系与验收标准。项目总纲见 [PLAN.md](PLAN.md)（第 8 节是当前阶段），
全部结果见 [RESULTS.md](RESULTS.md)。

**建议顺序：先做任务三**（纯本地 Python，零编译零排队），再做任务一、
二（要改 Fortran 补丁并在服务器编译回归，可以一次编译带两个改动）。

## 背景五分钟

- 项目一句话：LEOS 从 DFT 计算自动判定每个原子的氧化态——VASP 补丁
  把电子波变换成局域 Wannier 电子云写进 EOSCAR 文件，pyleos 把"每团
  电子归哪个原子"解成带容量约束的最小费用流（运输问题）。
- 两本账，别弄混：EOSCAR 里的逐原子表是 VASP 原位**贪心**（对照
  基线，可能出现荒谬值，那是故意保留的失效样本）；**正式结论**是
  pyleos 运输判定（`python3 pyleos/cli.py case <目录> -v`，纯标准库
  无依赖）。
- 代码布局：`patch/leos_v3.patch` 是单一 diff（触及 .objects、
  main.F、mlwf.F、eos.F 四个文件）；`patch/v3/eos.F` 是新增源文件
  eos.F 的当前全文，读代码从它入手最方便，**改动时两处要同步**。
  Python 侧在 `pyleos/`。
- 服务器：`ssh sai-pod`。项目树 `~/LEOS/`（索引见服务器
  `~/LEOS/README.txt`）；补丁的 apply → 编译 → 回归全流程记录在
  [../patch/README.md](../patch/README.md)，现成工作区
  `~/LEOS/leos_merge/{vanilla,rebuild,regress}`（原版树 / 重建区 /
  回归区）。16V100 分区单作业至少 4 GPU，长作业用 `--qos=huge-gpu`，
  GPU 作业不要写 `--cpus-per-task`/`--mem*`（站点钩子会自动杀）。
- 改动约定：改补丁 → 干净原版树 apply → 编译 → 回归对比（正常案例
  输出逐字节不变，目标案例按验收标准变好）→ 文档销账（见文末交付
  清单）→ 提交推送（提交信息风格看 `git log` 最近几条，结尾带
  Co-Authored-By）。

## 任务一：NELECT-mismatch 先写诊断再跳过（Fortran）

**现状**：体系丢电子超过阈值（`NELECT_ERR > EOS_NELECT_ABORT`，
默认 0.35）时，eos.F 直接放弃分析并跳过后续输出——结果是**恰恰
最需要诊断的体系拿不到逐 (spin, k点, 能带) 的 N_LOST 电子丢失
定位块**。轻度丢失（超 WARN 不超 ABORT）不受影响。

**根因证据**（RESULTS 3.14）：师弟旧目录里凡是手写 wannier90 输入
缺冻结窗口的体系全部丢电子触发中止，当时只报一个丢失总数、没有
定位，诊断只能靠外部取证绕路。根修复（AUTO_WINDOW）已上线，但
"丢了电子却不告诉你丢在哪"的输出顺序问题还在。

**位置**：`grep -n EOS_NELECT_ABORT patch/v3/eos.F` 找到 abort
分支；N_LOST 定位块的写出逻辑在正常路径里（EOSCAR 头部已有
`N_LOST_STATES` / `N_LOST_E` 字段，逐层明细跟在其后）。

**改法**：把 N_LOST 定位块的写出移到 abort 判断之前——先把"丢在
哪个自旋、哪个 k 点、哪条带"完整写进 EOSCAR/print-out，再决定
跳不跳过分析。保持现有字段格式不变（pyleos 解析器已认这个格式）。

**复现与验收**：取服务器 `~/LEOS/legacy_rerun/mumnox/` 的输入，
去掉冻结窗口（关 `LWANNIER90_AUTO_WINDOW`，win 里不写
dis_froz_max）重跑，故意触发 abort——改后应能在输出里看到完整的
逐层丢失定位；同时跑一个正常案例（CO2，输入在
`~/LEOS/leos_merge/regress/`）确认输出与改前逐字节一致。

## 任务二：num_wann 按自旋通道分别设定（Fortran）

**现状**：`NUM_WANN` 是单一标量，两个自旋通道共用（自动设定时取
两通道占据需求的较大者）。磁性体系里少数自旋通道的占据数更少，
共用值必然超出——**每个超出的 Wannier 函数就是一个幽灵**（名义
容量比真实电子多）。

**根因证据**（RESULTS 3.14、第 7 节）：VIRCEF 的 spin2 通道恒定
+2 个幽灵；Ce18O36Co4 三像 ghost 约 6-10。幽灵会被 pyleos 的总量
封顶自动排除，结论通常不受影响，但 RESULTS 5.2 记录了它的系统性
偏置风险（上游丢失叠加幽灵时，电正性金属可能被高估）——根治要在
上游让幽灵不产生。

**位置**：`patch/v3/eos.F` 约 175-191 行（NUM_WANN 的自动设定
`NUM_WANN = NREQ` 与校验），以及 `WF_OCC(MLWF%NUM_WANN, MLWF%ISPIN)`
的分配处；wannier90 库对每个自旋通道本来就是独立调用的（mlwf.F
段的自旋循环），技术上支持每通道用自己的值。

**改法方向**：把单一 NUM_WANN 改为按自旋各设（数组，或自旋循环内
的局部值），每通道用自己的占据需求 NREQ(ISP)。EOSCAR 的 #WF 块
本身带 spin 列，格式不用动；注意用户在 INCAR 显式给 NUM_WANN 时
的兼容行为（可保持"显式值 = 两通道上限"的旧语义并在输出里注明）。

**复现与验收**：VIRCEF 重跑（输入在 `~/LEOS/legacy_rerun/vircef/`），
spin2 幽灵应从 +2 降到 0；Ce18O36Co4_IS（同目录 ce_is）ghost 从 6
明显下降；非磁性案例（CO2、As4S4）回归输出不变。幽灵消除后，
顺手复核 LASWEH/VIRCEF 的金属价态是否仍正确（5.2 偏置应同步缓解）。
注意 ce 系体系大，作业要 8 GPU（此前 4 GPU 出过 OOM）。

## 任务三：仲裁阈值在含 C-H 基准上重校准（纯 Python，先做这个）

**现状**：档 2 仲裁模式的两个阈值 `ARBITER_TAU = 1.0`、
`ARBITER_PMARGIN = 0.05`（[../pyleos/transport.py](../pyleos/transport.py)
185-186 行）是首版拍定的，未系统扫描。RESULTS 3.8 有警示：几何
近平局的 C-H 键上，投影成分可能把氢的电子误判（氢的投影范数
天然小，成分裁决对它不公平）。

**机制**（transport.py 268-279 行）：档 1 几何代价的前两名候选
相差小于 TAU 时才启用仲裁；第二名的投影成分领先第一名超过
PMARGIN 才允许翻案。TAU 越大仲裁越常介入，PMARGIN 越小翻案越容易。

**改法**：在含 C-H 键的案例（MUMNOX、VIRCEF、LASWEH——EOSCAR 在
`data/v3/legacy_rerun/`——加上 m4 两波里的有机分子案例）上扫
TAU × PMARGIN 网格（建议 TAU 0.5/1.0/1.5/2.0 × PMARGIN
0.05/0.10/0.20），同时用 `script/full_bench.py --tier2 arbiter`
跑全量回归。选择标准：全量 dev 不劣于当前值，且 C-H 键零错判
（氢保持 +1，不被翻案成负值）。把选值依据和扫描表写进 RESULTS 3.8。

**验收**：74+ 案例全量回归通过数持平或更好；legacy 三案例的 C-H
判定全对；新阈值写进 transport.py 并附一行注释说明来源。

## 别动的东西

服务器上：`~/zack`（师弟原始数据，是 3.14 取证原料）、`~/srzhan`
（他人）、`~/apps`、`~/potcars`、`~/models`、`~/env`（通用资源，
脚本按绝对路径引用）。生产二进制 `~/apps/vasp.6.4.3_leos_v2` 在
回归通过前不要覆盖。本地：`data/v1/`（历史基线，只读）。
`.gitignore` 已拦截 WAVECAR/CHGCAR/OUTCAR 等重文件——新回归数据
只收 EOSCAR/输入级小文本进 git。

## 交付清单（做完逐项勾）

1. 三条改动落进 `patch/leos_v3.patch`（与 `patch/v3/eos.F` 同步）
   与 `pyleos/transport.py`；
2. 回归证据入库：目标案例改善 + 正常案例不变，数据放
   `data/v3/`（小文本），运输判定输出照
   `data/v3/legacy_rerun/pyleos_rerun_20260709.txt` 的先例存档；
3. 文档销账：RESULTS 4.2 表各加一行改进记录、第 7 节三条负债标记
   已解决、3.8 补仲裁扫描结论；PLAN 8.1 第 3 条更新状态；
4. 提交并推送（惯例见 git log）。
