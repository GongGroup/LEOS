# pyleos

pyleos 读 LEOS VASP 补丁输出的 EOSCAR，把每个 Wannier 函数（WF）的电子分配给
原子，给出整数氧化态、质量警告和可靠性指标。分配写成带容量约束的最小费用流
（运输问题），全局最优、结果为整数。只用 Python 标准库，没有第三方依赖。

方法说明见 [docs/PLAN.md](../docs/PLAN.md)，输出的读法见其中"怎么读结果"一节。

## 用法

```bash
# 单个算例。目录中需要 POSCAR、EOSCAR，以及 wannier90.wout
# （自旋极化时为 wannier90.1.wout 和 wannier90.2.wout）
python3 -m pyleos case cases/as4s4/v3 -v

# 历史回归集：对 36 个 v1 运行比较运输与贪心，并对照参考答案
python3 -m pyleos regress

# 基准（见下文"基准"）
python3 -m pyleos bench full --no-reliability --out /tmp/bench.tsv
```

`python3 pyleos/cli.py ...` 与 `python3 -m pyleos ...` 等价。

`case` 的选项：

| 选项 | 默认 | 作用 |
|---|---|---|
| `--beta B` | 10.0 | 电负性项的权重。小于约 8 时 X–H 键会被判错且没有警告，见 RESULTS.md"参数选择" |
| `--eneg allen\|pauling` | allen | 电负性标度 |
| `-v` | | 输出逐原子结果 |
| `--fragment 2,3,4` | | 把这些原子（从 1 开始编号）合成一个整体接收电子，可重复使用 |
| `--drop-stray` | | 排除未定域的 WF（stray）后再分配，用于判断错误来自输入还是分配。排除后总电子数变少，元素平均值不能再与参考答案比较 |
| `--gamma G` | 0.0 | 线性叠加投影权重：代价减去 G 乘以归一化的投影权重。这种方式对参数敏感、偏向重原子，不推荐；推荐的阈值规则见下文"作为库使用" |

`regress` 的选项：`--beta`、`--eneg`、`--reliability`（同时计算 R）。

## 基准

`python3 -m pyleos bench` 有三个子命令，案例目录都从 cases/INDEX.tsv 查得：

| 子命令 | 内容 |
|---|---|
| `full --out FILE` | 74 案例基准，输出与 results/bench74_golden.tsv 同格式（`--out` 必填，防止误覆盖金标）；`--greedy` 改为给贪心列打分，`--no-reliability` 跳过 R |
| `clean` | clean 集：元素平均 dev 加逐原子正确率；`--method` 指定分配插件（默认 methods/baseline.py） |
| `test --set regression\|homonuclear\|both --method PLUGIN` | 预注册的留出测试集；`--selftest` 只用开发集自检打分路径 |

分配插件是定义了 `assign(case)` 的 .py 文件，放在 methods/；methods/frozen/
下是已冻结、不再修改的候选（sha256 记在 sets/test_set.json）。插件看不到其他
方法的结果和文件系统，规则见 bench.py 开头。

基准定义在 sets/：reference_os.json（74 案例的参考氧化态，按 dataset、oxides、
wave1、wave2 分组）、clean_set.json（逐原子标注）、test_set.json（留出测试集）。

## 作为库使用

```python
from pathlib import Path
from pyleos.parsers import load_case
from pyleos.transport import assign_case, reliability_margins

case = load_case(Path("cases/as4s4/v3"))
res = assign_case(case)                      # 默认：beta=10，Allen，几何 + 电负性
res = assign_case(case, tier2="arbiter")     # 投影辅助判定，需要 EOSCAR 的 #WP 块
res = reliability_margins(case, res=res)     # 计算 R_site / R_elem
```

投影辅助判定的阈值是 transport.py 中的 ARBITER_TAU（1.0）和
ARBITER_PMARGIN（0.20）。

## 模块

| 文件 | 内容 |
|---|---|
| parsers.py | 读 POSCAR、EOSCAR（含 #WF、#WP 块）和 wannier90 .wout，组装成 CaseData |
| capacity.py | 各元素的容量上下限 OS_BOUNDS，Allen 和 Pauling 电负性表 |
| transport.py | 最小费用流分配、stray 和 ghost 诊断、可靠性指标 R |
| cases.py | 读 cases/INDEX.tsv（案例登记表） |
| bench.py | 74 案例、clean 集和留出测试集三个基准 |
| cli.py | 命令行入口（`python3 -m pyleos`） |
| batch.py | 集群批跑驱动（见下） |
| xsd2poscar.py | 把 Materials Studio 的 .xsd 转成 POSCAR |

EOSCAR 解析的一个兼容处理：官方 Xe 赝势的 VRHFIN 字段写作 "X"，会使 EOSCAR
中的元素名被截断，parsers.py 按前缀匹配元素名。

## 模型

- 图：源 → WF 节点 → 原子节点 → 汇。闭壳层（ISPIN = 1）每个 WF 提供 2 个
  电子；自旋极化（ISPIN = 2）每个通道每个 WF 提供 1 个电子，两个通道在同一张
  图里共享原子容量。
- 代价：`d²(WF 中心, 原子) + β·(1 − χ/χmax)`，χ 默认为 Allen 电负性。
- 原子容量 [n_min, n_max] 来自 capacity.py 的 OS_BOUNDS，按价层物理极限设定。
  半芯壳层必须记满，用负代价边实现 n_min。
- 总流量不超过 NELECT，多余的 WF（ghost）由求解器丢掉代价最高的部分。
- 电子数和容量都是整数，所以最优流是整数，氧化态不需要舍入。
- stray 判据：spread > 3 Å²，或中心离所有原子都超过 1.6 Å。

求解器实现中需要注意的三点：

1. 最短路的势（potential）更新必须覆盖没有被访问到的节点，否则会出现负的
   约化代价，算法退化到指数时间。
2. 大数值奖励与浮点容差的数量级要拉开（1e3 与 1e-9），否则舍入残差会被当成
   改进，造成死循环。
3. 迭代次数设上限，异常时抛出明确错误，而不是无限循环。

## batch.py：集群批跑

batch.py 只依赖标准库，可以单独复制到集群上运行；在仓库内也可用
`python3 -m pyleos batch ...` 调用。

```bash
python3 batch.py stage   --dest ~/leos_batch1 --cases cases.json
python3 batch.py submit  --dest ~/leos_batch1
python3 batch.py status  --dest ~/leos_batch1
python3 batch.py retry   --dest ~/leos_batch1
python3 batch.py harvest --dest ~/leos_batch1 --out ~/leos_batch1_results
```

cases.json 的格式：`[{"name": "co2", "src": "<含 POSCAR、POTCAR、KPOINTS 的目录>",
"incar": {"ENCUT": 500}}, ...]`。每个算例目录会得到 LEOS INCAR 模板（合并
算例自己的设置）、一个 4 GPU 的 sbatch 脚本，状态记录在
`<dest>/batch_state.json`。

失败分类（status 和 harvest 输出）：

| 状态 | 含义 | retry 的处理 |
|---|---|---|
| FAIL_DFT | SCF 崩溃（EDDDAV、ZHEGV、"I REFUSE"）或 OUTCAR 不完整 | 先 ALGO = Normal，再 ALGO = All 且 EDIFF = 1e-5 |
| FAIL_W90 | DFT 正常，但 wannier90 没有输出 Final State | 关闭 EOS_SCDM_ONLY，NBANDS 放大到 1.3 倍 |
| FAIL_LEOS | Wannier 化正常，但没有可用的 EOSCAR（缺失，或带 LEOS_ABORTED 标记） | EOS_NELECT_ABORT 放宽到 0.5 重跑，用于诊断 |
| OK | EOSCAR 存在且没有 LEOS_ABORTED 标记 | |

wannier90 报 "Not enough bvectors" 时归为 FAIL_W90（是 k 网格问题，不是 SCF
问题），处理方法见 patch/README.md"常见失败"。

两种附加模式：

- **分级筛查**（`stage --mode tiered`，之后运行 `escalate`）：不含 4d、5d、
  f 区元素的算例先用 SCDM-only 快速模式，含这些元素的直接用完整 MLWF。运行
  结束后，`escalate` 把 EOSCAR 带质量警告（stray、dNe 偏差、ghost 过多）的
  快速模式算例改为完整 MLWF 重新排队。
- **离线 wannier90**（`stage --offline-w90`，之后反复运行 `w90`）：GPU 作业
  只导出 .mmn/.amn（`LWRITE_MMN_AMN = .TRUE.`，`LWANNIER90_RUN = .FALSE.`），
  wannier90 在 DSPRHBM CPU 分区上运行。`w90` 命令推进每个算例的状态：导出完成
  → 生成每个自旋的 .win（tools/offline_w90_prep.py，单文件部署时与 batch.py 放在一起）并提交一个 CPU 作业
  （两个自旋并发）→ 两个 .wout 都有 Final State 后复制回算例目录，供 pyleos
  读取。只导出的运行不产生可用的 EOSCAR，离线模式按导出文件和 .wout 判断状态。

**batch.py 开头的 VASP_BIN、OFFLINE_VASP_BIN、W90_BIN 三个常量写死了 SAI
集群家目录迁移前的路径（/home/sjtu-gongxueqing/zhangxianyang/...），其中
VASP_BIN 还指向旧的 vasp.6.4.3_leos_v2。使用前要改成当前路径，见
patch/README.md"集群部署"。**

## 测试

```bash
python3 -m unittest discover pyleos/tests -v   # 21 个测试
```

| 文件 | 内容 |
|---|---|
| test_mcf.py | 最小费用流求解器与分配 |
| test_batch.py、test_batch_tier.py、test_batch_offline.py | batch.py 的失败分类、重试、分级筛查和离线模式 |

## 已知限制

- 等价位点之间（同核键两端、簇内等价原子）的归属不唯一：元素平均正确，单个
  原子的值取决于求解器落在哪个最优解上。用 `--fragment` 处理。
- 容量只按物理极限设定；按化学常见价态收紧容量的诊断模式没有实现。
