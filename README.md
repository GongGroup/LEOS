# LEOS：从 Wannier 函数判定固体中的氧化态

LEOS（Lattice Effective Oxidation States）根据 VASP 计算得到的 Wannier 函数，
自动给出晶体中每个原子的整数氧化态。它把"每份电子归哪个原子"写成带容量约束的
最小费用流问题，所有电子一起做全局最优分配，结果自动是整数。

仓库包含三部分：

- **VASP 补丁**（[patch/](patch/)）：在 VASP 6.4.3 内完成 Wannier 化和氧化态
  分配，输出 EOSCAR 文件。
- **pyleos**（[pyleos/](pyleos/)）：读 EOSCAR 做分配、诊断和可靠性评估，只用
  Python 标准库。
- **gpuw90**（[gpuw90/](gpuw90/)）：用 PyTorch 实现的 wannier90 解纠缠和局域化，
  可在 GPU 上运行。

## 主要结果

| 测试 | 结果 |
|---|---|
| 74 案例基准（72 个有参考答案） | 运输 64/72（89%），逐个 WF 独立判定的贪心方法 51/72（71%） |
| Comas-Vilà 2025 的 40 个固体 | 全部给出判定 |
| Li2O2 脱锂系列 | O 的平均氧化态 −1.00 / −0.94 / −0.88，与电荷计数一致 |
| 富锂正极 180 帧 | 有序体系过渡金属平均价态等于 4 + 2x，50 ps 内无漂移 |

数据来源、未通过案例的原因和其余结果见 [docs/RESULTS.md](docs/RESULTS.md)。

## 快速开始

pyleos 只需要 Python 3：

```bash
# 分析一个算例（目录中需要 POSCAR、EOSCAR 和 wannier90 的 .wout）
python3 pyleos/cli.py case data/v3/m4_wave2/as4s4 -v

# 74 案例基准，约 25 秒；一定要给 --out，否则会覆盖基准参考表
python3 script/full_bench.py --no-reliability --out /tmp/bench.tsv

# 单元测试
python3 -m unittest discover pyleos/tests -v
```

编译和运行补丁版 VASP 见 [patch/README.md](patch/README.md)；gpuw90 需要
torch 和 numpy，见 [gpuw90/README.md](gpuw90/README.md)。

## 文档

| 文件 | 内容 |
|---|---|
| [docs/PLAN.md](docs/PLAN.md) | 问题、方法、怎么读输出、主要设计决定、文献、工程约定 |
| [docs/RESULTS.md](docs/RESULTS.md) | 全部测试与应用结果，每个数字注明数据文件 |
| [CONTEXT.md](CONTEXT.md) | 术语表 |
| [docs/testset_findings.md](docs/testset_findings.md) | 测试集结构来源核查，逐案例记录在 testset_registry.tsv |
| [docs/comas_vila_2025_notes.md](docs/comas_vila_2025_notes.md) | 参考方法 Comas-Vilà 2025 的笔记与对比 |
| [patch/README.md](patch/README.md) | 补丁安装、INCAR 参数、运行注意、版本记录 |
| [pyleos/README.md](pyleos/README.md) | pyleos 用法、模块、批跑脚本 |
| [gpuw90/README.md](gpuw90/README.md)、[gpuw90/DESIGN.md](gpuw90/DESIGN.md) | gpuw90 用法与设计 |

## 目录

| 目录 | 内容 |
|---|---|
| data/structures/ | 48 个测试结构及其 Bader 分析 |
| data/v1/ | 早期 v1 程序的历史数据，74 案例基准中历史回归集的来源（见 [data/v1/README.md](data/v1/README.md)） |
| data/v3/ | 当前补丁产生的数据：m4_wave1、m4_wave2（扩展基准）、comas_gap（Comas-Vilà 对照集）、m5a_li2o2、m5c_ceo2、legacy_rerun 等，以及基准参考表 reliability_74.tsv |
| data/m5b/ | 富锂正极三体系的逐帧价态序列 |
| script/ | 基准、参数扫描、轨迹分析、离线 wannier90 输入生成等脚本 |
| docs/structures_fix/ | 替换后的正确结构及来源记录 |

data/ 只保存小文本文件（EOSCAR、输入、.wout、汇总表）；波函数等大文件在集群上，
位置见 docs/PLAN.md"工程约定"。

## 待办

https://github.com/GongGroup/EOS/issues
