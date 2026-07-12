# LEOS：Automatically determining oxidation states in solids from Wannier functions


| 位置 | 内容 | 入口 |
|---|---|---|
| [docs/PLAN.md](docs/PLAN.md) | 项目计划：背景、核心概念、决策记录、当前状态与待办 | — |
| [docs/RESULTS.md](docs/RESULTS.md) | 全部结果：测试设计、主要数据、逐案例解读、术语表 | 看成果从这里 |
| [pyleos/](pyleos/) | Python 分析包：运输分配器、解析器、批跑驱动（纯标准库，21 个单元测试） | `pyleos/README.md` |
| [patch/](patch/) | VASP 6.4.3 补丁，正式发布物为合并版 **`leos_v3.patch`**（含 v3.1 全部修复；从零安装与验收记录见 README） | `patch/README.md` |
| [gpuw90/](gpuw90/) | GPU Wannier 化独立工具（PyTorch，CPU/CUDA 同码）：局域化 + 解纠缠两阶段，金标对照与 V100 基准 | `gpuw90/README.md` |
| [data/](data/) | 全部数据：`structures/`（48 个测试结构库）、`v1/`（历史基线）、`v3/`（当前版本全部批次：74 案例金标、对照集、结构修复、复验等）、`m5b/`（富锂正极 180 帧轨迹分析） | `data/v3/reliability_74.tsv` |
| [script/](script/) | 独立脚本：全量基准 `full_bench.py`、β 扫描、仲裁参数扫描、离线 w90 输入生成、MACE 轨迹链路等 | 各脚本 docstring |
| docs/ 其余 | 测试集来源登记（`testset_registry.tsv` + `testset_findings.md`）、修正结构包（`structures_fix/`）、参考文献（`refs/`，版权文件不入 git） | — |

## 快速上手

```bash
# 分析一个已归档的算例（无需任何依赖，纯标准库）
python3 pyleos/cli.py case data/v3/m4_wave2/as4s4 -v

# 对全部历史回归集跑一遍运输 vs 贪心对比
python3 pyleos/cli.py regress          # 默认扫 data/v1

# 74 案例全量基准（约 20 秒；加 --greedy 得到贪心基线口径）
python3 script/full_bench.py --no-reliability --out /tmp/bench.tsv
```

