# LEOS：从 Wannier 函数自动判定固体中的氧化态

本仓库包含 LEOS 方法（VASP 补丁 + Python 分析包）的全部代码、测试数据和文档。
项目目标是一篇方法学期刊论文，核心主张是把"电子归属哪个原子"这个问题
表述为一个有全局最优解的运输问题。

## 目录导航

| 位置 | 内容 | 从哪看起 |
|---|---|---|
| [docs/PLAN.md](docs/PLAN.md) | 项目计划：背景、核心概念、决策记录、里程碑 | **新读者从这里开始** |
| [docs/RESULTS.md](docs/RESULTS.md) | 全部结果：测试设计、主要数据、逐案例解读、术语表 | 看成果从这里 |
| [pyleos/](pyleos/) | Python 分析包：运输分配器、解析器、批跑驱动 | `pyleos/README.md` |
| [patch/](patch/) | VASP 6.4.3 补丁：`leos.patch`（v2 基线）+ `v3/`（当前版本增量） | `patch/v3/README.md` |
| [data/](data/) | 全部数据：`structures/`（48 个测试结构库：XSD + Bader 归档）、`v1/`（历史基线运行目录）、`v3/`（当前版本四批精简结果） | `data/v1/RESULTS_INDEX.md` |
| [script/](script/) | 独立小工具（Bader 统计、运输原型） | — |
| docs/refs/ | 参考文献 PDF（版权文件，不入 git） | — |

## 快速上手

```bash
# 分析一个已归档的算例（无需任何依赖，纯标准库）
python3 pyleos/cli.py case data/v3/m4_wave2/as4s4 -v

# 对全部历史回归集跑一遍运输 vs 贪心对比
python3 pyleos/cli.py regress          # 默认扫 data/v1
```

## 服务器（sai-pod）

- 生产二进制：`~/apps/vasp.6.4.3_leos_v2/bin/{vasp_std,vasp_gam}`（v3 特性全开）
- 工作树：`~/LEOS/`（原始大文件仓：inputs / runs / results / tools）。本地 `data/v3/` 是其 results 的精简镜像（EOSCAR/POSCAR 级），大文件（wout/波函数类）只在服务器保留
- 批跑：`python3 ~/LEOS/tools/leos_batch.py --help`
