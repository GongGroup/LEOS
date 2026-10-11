# v2 → v3 增量材料（存档）

本目录是 v3 开发时相对 v2（patch/leos.patch）的增量材料，内容已全部并入
patch/leos_v3.patch。安装请直接用 leos_v3.patch，见 [../README.md](../README.md)。

v3 在 SAI 集群的生产树 `vasp.6.4.3_leos_v2` 上开发和验证（编译作业 593892，
验证作业 593893、593894，2026-07-07）。

## 文件

| 文件 | 内容 |
|---|---|
| eos.F | 全量替换 v2 的 eos.F |
| apply_v3_mlwf.py | 向 mlwf.F 注入 LEOS 变量，以及 .win 写入段中 SCDM-only 和引导中心的处理 |
| apply_v3_main.py | 向 main.F 注入轨迹模式接线：按 EOS_DO_THIS_STEP 逐步做 Wannier 化和分析（只开 LWANNIER90 时仍只在第一步运行），并去掉每步无条件的 MLWF_FREE |

eos.F 相对 v2 新增的功能：

- `EOS_PROJ_WEIGHTS`：每个 WF 在各原子上的投影范数（#WP 块），供投影辅助判定
  使用。曾试过按 QION 增广差加权，软赝势通道的对比度只有 1.35:1，改用投影
  范数后约 5:1，因此放弃前者。
- `EOS_SCDM_ONLY`：跳过 wannier90 的迭代局域化。CO2 上氧化态不变，wannier90
  耗时 2 秒 → 0.023 秒。
- `EOS_WARM_START`（默认开，EOS_EVERY_N > 0 时生效）：把上一步的 WF 中心作为
  引导中心传给 wannier90。

## 开发时的观察

- CO2 在 92 K 下 3 步 MD：贪心结果在第 3 步漂移（C 由 +4 变为 +3 或 +2）。
- 小体系上引导中心的开销大于收益（0.58 秒 → 1.3 秒）。
- MD 中部分占据使 num_wann 增大（32 → 40）；pyleos 用 NELECT 总量上限排除
  多出的 WF。
