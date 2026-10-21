# Lattice Effective Oxidation States (LEOS)

LEOS 根据 VASP 计算得到的 Wannier 函数，
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

pyleos 只需要 Python 3，所有功能都从 `python3 -m pyleos` 进入：

```bash
# 分析一个算例（目录中需要 POSCAR、EOSCAR 和 wannier90 的 .wout）
python3 -m pyleos case cases/as4s4/v3 -v

# 74 案例基准，约 25 秒；--out 必填，金标表不会被误覆盖
python3 -m pyleos bench full --no-reliability --out /tmp/bench.tsv

# 单元测试
python3 -m unittest discover pyleos/tests -v
```

编译和运行补丁版 VASP 见 [patch/README.md](patch/README.md)；gpuw90 需要
torch 和 numpy，见 [gpuw90/README.md](gpuw90/README.md)。

## 我想……

| 想做的事 | 去哪里 |
|---|---|
| 分析一个体系 | `python3 -m pyleos case cases/<体系>/<运行>`；体系和运行都列在 [cases/INDEX.tsv](cases/INDEX.tsv) |
| 找某个体系的全部计算 | `cases/<体系>/`，一个化学体系一个目录，每个子目录是一次计算 |
| 跑基准、对比金标 | `python3 -m pyleos bench full\|clean\|test`，见 [pyleos/README.md](pyleos/README.md)"基准" |
| 查一个成绩数字的出处 | [docs/RESULTS.md](docs/RESULTS.md)，每个数字都标了 results/ 或 cases/ 下的文件 |
| 在集群上批量算 | `python3 -m pyleos batch ...`，见 [pyleos/README.md](pyleos/README.md)"batch.py" |
| 理解方法和设计取舍 | [docs/PLAN.md](docs/PLAN.md)，术语见 [CONTEXT.md](CONTEXT.md) |
| 重画论文插图 | [paper/](paper/)，每张图一个脚本加图注 |

## 目录

| 目录 | 内容 |
|---|---|
| cases/ | 全部计算数据：`cases/<体系>/<运行>/`，登记表 INDEX.tsv，命名规则见 [cases/README.md](cases/README.md) |
| results/ | 全部汇总表：金标 bench74_golden.tsv、贪心对照、参数扫描、对照集与测试集结果、富锂正极时间序列（cathode_md/） |
| pyleos/ | Python 包：分配、基准（bench.py）、案例登记（cases.py）、批跑（batch.py）；基准定义在 sets/，分配插件在 methods/ |
| patch/ | VASP 补丁 |
| gpuw90/ | GPU Wannier 化工具 |
| tools/ | 研究用脚本：参数扫描、MD 与轨迹、Bader、离线 wannier90 输入生成 |
| docs/ | 方法、结果、测试集核查、参考方法笔记 |
| paper/ | 论文插图脚本与图注 |

cases/ 只保存小文本文件（EOSCAR、输入、.wout）；波函数等大文件在集群上，
位置见 docs/PLAN.md"工程约定"。

## 待办

https://github.com/GongGroup/EOS/issues
