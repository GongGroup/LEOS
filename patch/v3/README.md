# LEOS v3 增量（相对 patch/leos.patch = v2）

v3 在 sai-pod 生产树 `~/apps/vasp.6.4.3_leos_v2` 上开发并验证
（编译作业 593892，验证作业 593893/593894，2026-07-07）。

## 内容

1. **eos.F**（本目录，全量替换 v2 的 eos.F）
   - `EOS_PROJ_WEIGHTS`：每 WF 在各原子的投影仪范数指纹（#WP 块），
     喂 pyleos 档 2 运输代价。QION 增广差加权已试验并放弃
     （软赝化通道对比度仅 1.35:1；投影范数恢复到 ~5:1）。
   - `EOS_SCDM_ONLY`：跳过 W90 最小化（CO2 验证：OS 不变，W90 时间
     2s → 0.023s）。
   - `EOS_WARM_START`（默认开，EOS_EVERY_N>0 生效）：把上一步 WF 中心
     交给 mlwf 作为 guiding centres。
2. **mlwf.F 注入**（apply_v3_mlwf.py）：LEOS 钩子变量 + win 写入段的
   SCDM-only/guiding 注入。
3. **main.F 注入**（apply_v3_main.py）：轨迹模式接线——老集成写死
   `IF (NSTEP == 1)`，EOS_EVERY_N 特性在 v1/v2 中从未真正执行；
   v3 起按 EOS_DO_THIS_STEP 逐步 Wannier 化+分析（纯 LWANNIER90
   用户保持首步一次的原行为），并修复每步无条件 MLWF_FREE 的
   悬空引用隐患。

## 已知观察（进论文素材）

- 3 步 92K CO2 MD：贪心 OS 在第 3 步漂移（C +4 → +3/+2），朴素贪心
  在轨迹上不稳定 → 运输后处理的价值点。
- 小体系上 guiding 开销 > 收益（0.58s → 1.3s）；warm start 的收益
  场景是大体系/难收敛（As4S4 类），待基准。
- MD 中 num_wann 随部分占据虚胖（32 → 40）；pyleos 的 NELECT 上限
  流自动排除幽灵电子。
