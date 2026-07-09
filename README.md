# LEOS：从 Wannier 函数自动判定固体中的氧化态

本仓库包含 LEOS 方法（VASP 补丁 + Python 分析包）的全部代码、测试数据和文档。
项目目标是一篇方法学期刊论文，核心主张是把"电子归属哪个原子"这个问题
表述为一个有全局最优解的运输问题。

（2026-07-09 状态：两个阶段的计算与验证已全部完成并推送远端，项目
处于收尾与论文写作阶段——收尾清单与写作路线见 docs/PLAN.md 第 8 节。）

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

## 结果怎么读（两本账，别弄混）

- 每个算例目录里 EOSCAR 的逐原子表是 VASP **原位贪心**写的——它是
  对照基线（老方法），直接翻它可能看到超出物理容量的荒谬值（比如
  一个氢名下 −2），那正是论文要展示的失效模式，不是最终结论。
- **正式结论是运输判定**，由 pyleos 后处理产生：
  `python3 pyleos/cli.py case <目录> -v`（`--fragment i,j,k` 做片段
  合并记账）。重要体系的运输输出已固化为存档文本，例如
  `data/v3/legacy_rerun/pyleos_rerun_20260709.txt`（Ce 三像逐原子账本）。
- 质量三警告（stray / ghost / dNe）与可靠性 R（改判裕度）的朴素
  读法见 RESULTS 第 5 节。

## 服务器（sai-pod）

- 生产二进制：`~/apps/vasp.6.4.3_leos_v2/bin/{vasp_std,vasp_gam}`（v3 特性全开）
- 项目文件统一入口：`~/LEOS/`（2026-07-09 归拢——第一阶段工作树
  inputs / runs / results / tools，加上第二阶段全部工作目录
  legacy_rerun / gap_* / m5b_* / polaron_li2o2 / fe3o4_* / ws_bench 等；
  散脚本在 tools/、日志在 logs/、作业脚本在 sbatch/；索引见服务器
  `~/LEOS/README.txt`）。本地 `data/v3/` 是集群结果的精简镜像
  （EOSCAR/POSCAR 级小文本），波函数类大文件只在服务器保留，
  .gitignore 已拦截此类重文件误入 git
- 通用资源（不在项目树内，勿动）：`~/apps`、`~/potcars`、`~/models`、`~/env`
- 批跑：`python3 ~/LEOS/tools/batch.py --help`（第二阶段批次驱动；
  第一阶段的 leos_batch.py 在同目录）
