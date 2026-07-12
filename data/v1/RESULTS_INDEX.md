# LEOS 结果索引（v1 历史基线）

> 迁移注记（2026-07-12）：本索引写于 v1 时代，原指向的本机工作区
> `~/Desktop/EOS4VASP/leos_work/` 已不存在——数据早已全部迁入本仓库：
> 汇总表在 `data/v1/oxide_leos_summary.tsv` 与
> `data/v1/eos_dataset/leos_vs_bader.tsv`，运行目录在 `data/v1/oxides/`
> 与 `data/v1/eos_dataset/`。下文正文保留当时的统计口径供追溯；
> 文中旧路径一律按上述仓库内路径读。v1 之后的全部结果见
> `data/v3/` 与 docs/RESULTS.md。

## 氧化物套件（历史输出）
- 汇总 TSV：`data/v1/oxide_leos_summary.tsv`（tag=default 为默认回归集合；tag=stress 为压力/实验用例）
- 原始输出目录：`data/v1/oxides/leos-oxide-*-*`（同一 case 可能多次运行；TSV 仅选取最新一次）

### 氧化物套件概况（以汇总 TSV 为准）
- 总用例：29（default=24, stress=5）
- OK：28；非 OK：1
- 非 OK 列表：
  - reo2 (NO_EOSCAR)

## EOS 数据集（LEOS vs Bader）
- 对比汇总 TSV：`data/v1/eos_dataset/leos_vs_bader.tsv`（按“结构×元素”逐行对齐 LEOS 与 Bader）
- 原始输出目录：`data/v1/eos_dataset/<structure>`（结构名会被净化，例如 `(Cu5I7)2-`→`_Cu5I7_2-`）

### 当前进度（以对比 TSV 为准）
- EOS 数据集总结构数：45
- 已产出汇总记录的结构数：10（OK：9；非 OK：1）

### 已完成/已记录结构一览
| 结构 | 状态 | Bader 总净电荷 | $\Delta N_{\mathrm{e}}$ | LEOS 氧化态范围（按元素汇总） |
| --- | --- | --- | --- | --- |
| (Cu5I7)2- | OK | -2.0000 | 0.000000 | Cu(+1)x5, I(-1)x7 |
| As4S4 | OK | +0.0040 | 0.000000 | S(-8~0)x16, As(+3~+5)x16 |
| AuClB(C6H5)3{P(C6H5)3}2 | OK | -0.0030 | 0.000000 | H(+1)x32, B(+3)x1, C(-3~0)x42, P(+5)x2, Cl(-1)x1, Au(-1)x1 |
| B6H10 | OK | -0.0040 | 0.000000 | H(-1)x40, B(0~+2)x24 |
| BH6N | OK | -0.0040 | 0.000000 | H(-1~+1)x24, B(+3)x4, N(-3)x4 |
| CH2F2 | OK | +0.0000 | 0.000000 | H(+1)x2, C(0)x1, F(-1)x2 |
| CO | OK | +0.0000 | 0.354438 | C(+1~+2)x4, O(-2~-1)x4 |
| CO2 | OK | -0.0000 | 0.000000 | C(+4)x4, O(-2)x8 |
| CSe2 | FAIL_NO_EIGENVAL | +0.0010 | - | - |
| ReO3 | OK | +0.0000 | 0.000000 | O(-2~-1)x3, Re(+4)x1 |

### 需要重点复核的结构（基于简单规则自动标记）
- 规则：状态非 OK / $\Delta N_{\mathrm{e}}>0.10$ / 存在元素氧化态落在 [-4, +6] 之外。
- As4S4（OK, dNe=0.000000）
- CO（OK, dNe=0.354438）
- CSe2（FAIL_NO_EIGENVAL, dNe=-）

### 本地输出但未写入对比 TSV（需重跑或手动补记）
- Ca3AuN（NO_EIGENVAL；-）

### 继续全量批跑（历史命令，环境已不存在，仅存档）
```bash
orb -m Ubuntu-x64 bash -lc 'NP=14 WORK_ROOT=/Users/mastreina/Desktop/EOS4VASP/leos_work CLEANUP=1 bash /Users/mastreina/Desktop/EOS4VASP/leos_eos_dataset_tests.sh'
```
（v1 批跑环境为本机虚拟机，已随 EOS4VASP 工作区退役；现行批跑见 pyleos/batch.py 与集群 `~/LEOS/tools/`。）
