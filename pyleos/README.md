# pyleos

LEOS 的离线分析包：把 Wannier 电子对原子的归属表述为**容量约束最小费用流**
（运输问题），全局最优、整数解，替代 Fortran 原位路径的最近邻贪心。
纯标准库，无第三方依赖。

## 用法

```bash
# 单案例（case 目录需含 POSCAR、EOSCAR、wannier90[.1/.2].wout）
python3 pyleos/cli.py case leos_work/leos_eos_dataset/As4S4 -v

# 全量回归（对照 tests/reference_os.json 的理想氧化态）
python3 pyleos/cli.py regress --root leos_work

# 参数：--beta 电负性权重（默认 10.0；X-H 键的正确判定需要 β≥10，36 案例回归显示 β=10 无副作用）；--eneg allen|pauling（默认 allen）
```

## 模型

- 图：源 → WF 节点（ISPIN=1 每 WF 供 2e；ISPIN=2 每通道每 WF 供 1e，
  双通道同图共享原子容量）→ 原子节点 → 汇。
- 边代价（档 1）：`d²(WF中心, 原子) + β·(1 − χ/χmax)`，χ 为 Allen 电负性。
  档 2（PAW 投影权重代价）等 Fortran 侧 EOSCAR v2 的 `#WF` 明细落地后启用。
- 原子容量 `[n_min, n_max]` 由 `capacity.py` 的宽物理界生成
  （`OS_BOUNDS`，宁松勿紧；半芯壳层通过 n_min 强制满占据，用负代价边段实现）。
- 整数供给 + 整数容量 ⇒ 最优流自动是整数，氧化态无需舍入。
- 诊断：spread > 3 Å² 或中心离所有原子 > 1.2 Å 的 WF 标记为 stray
  （Wannier 化质量问题，分配层无法修复，应触发工作流层重跑）。

## 测试

```bash
python3 pyleos/tests/test_mcf.py
```

## 已知限制（V1）

- 对称等价类容量池、fragment 汇节点未实现（简并键对 2+0/1+1 分裂时
  解不唯一，元素级汇总不受影响）。
- per-WF 占据未知时按满占据近似（部分占据体系如 CO 需 EOSCAR v2）。
- 紧容量（化学先验）诊断模式未实现。
