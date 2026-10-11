# data/v1：v1 程序的历史数据

这批数据由项目早期的 v1 程序产生，分配方法是贪心。它是 74 案例基准中
"历史回归集"（36 个）的来源。下面表格中的氧化态都是 v1 贪心结果；用当前
方法分析同一份数据：

```bash
python3 pyleos/cli.py case data/v1/eos_dataset/CO2 -v
```

## 目录内容

| 路径 | 内容 |
|---|---|
| oxides/leos-oxide-<体系>-<随机后缀>/ | 氧化物和卤化物的运行目录，共 198 个；同一体系可能有多次运行 |
| oxide_leos_summary.tsv | 氧化物汇总，每个体系取最新一次运行。tag=default 为常规回归，tag=stress 为压力测试 |
| eos_dataset/<结构名>/ | 分子晶体的运行目录，含 Bader 分析 bader_summary.md。结构名中的特殊字符已替换，例如 (Cu5I7)2- 写作 _Cu5I7_2- |
| eos_dataset/leos_vs_bader.tsv | v1 贪心结果与 Bader 电荷的逐元素对比 |

两个汇总表的 workdir 列是原工作区路径
（/Users/mastreina/Desktop/EOS4VASP/leos_work/...），该目录已不存在，对应
数据就在本目录下。

## 记录时的状态

下表的"状态"只表示程序是否跑完、电子数是否守恒，不表示氧化态是否正确。
例如 As4S4 记为 OK，但其中 S 的氧化态低到 −8。

氧化物汇总（oxide_leos_summary.tsv）：29 个体系（default 24，stress 5），
28 个 OK；reo2 没有产出 EOSCAR。

分子晶体（leos_vs_bader.tsv 等）：

| 结构 | 状态 | Bader 总净电荷 | dNe | v1 贪心氧化态（按元素） |
|---|---|---|---|---|
| (Cu5I7)2- | OK | −2.0000 | 0.000000 | Cu(+1)×5，I(−1)×7 |
| As4S4 | OK | +0.0040 | 0.000000 | S(−8~0)×16，As(+3~+5)×16 |
| AuClB(C6H5)3{P(C6H5)3}2 | OK | −0.0030 | 0.000000 | H(+1)×32，B(+3)×1，C(−3~0)×42，P(+5)×2，Cl(−1)×1，Au(−1)×1 |
| B6H10 | OK | −0.0040 | 0.000000 | H(−1)×40，B(0~+2)×24 |
| BH6N | OK | −0.0040 | 0.000000 | H(−1~+1)×24，B(+3)×4，N(−3)×4 |
| CH2F2 | OK | +0.0000 | 0.000000 | H(+1)×2，C(0)×1，F(−1)×2 |
| CO | OK | +0.0000 | 0.354438 | C(+1~+2)×4，O(−2~−1)×4 |
| CO2 | OK | −0.0000 | 0.000000 | C(+4)×4，O(−2)×8 |
| CSe2 | FAIL_NO_EIGENVAL | +0.0010 | — | — |
| ReO3 | OK | +0.0000 | 0.000000 | O(−2~−1)×3，Re(+4)×1 |

Ca3AuN 当时缺 EIGENVAL，未写入对比表。

"CO"的胞内实际是两个乙烯二酮分子（O=C–C=O），不是一氧化碳，见
docs/testset_findings.md。

v1 的批跑环境是本机虚拟机，已经停用。现在的批跑见 pyleos/batch.py。
