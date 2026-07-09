# 第二阶段收尾交接（2026-07-07）

> **状态（2026-07-09）：本文四个事项已全部闭环**——极化子 Fe3O4（RESULTS
> 3.10）、测试集补齐 39/40（3.13）、M5-B 全链（3.12）、效率基准（3.7），
> 追加的旧版体系验证也已 6/6（3.14）。本文转为历史存档，当前状态以
> PLAN.md 第 5 节与 RESULTS.md 第 7 节为准。
>
> **2026-07-09 服务器目录迁移**：sai-pod 家目录已归拢，本文出现的所有
> `~/<name>` 路径现位于 `~/LEOS/<name>`（目录平移保名；散脚本在
> `~/LEOS/tools/`，日志在 `~/LEOS/logs/`，作业脚本在 `~/LEOS/sbatch/`）。
> 映射详情见服务器 `~/LEOS/README.txt`。通用资源（`~/apps`、`~/potcars`、
> `~/models`、`~/env`）与师弟原始目录 `~/zack` 未移动。

本文件写给接手集群作业与收割分析的 agent。已完成并入库的工作见
[PLAN.md](PLAN.md) 第 7 节与 [RESULTS.md](RESULTS.md) 3.8-3.11 节及
git 历史（e5b3f79..2cdd3f4+），此处只列**未闭环的四件事**和它们的
收割方法。集群均为 sai-pod；两集群资源位置的说明在用户级记忆
`sai-cluster-resources` 里也有一份。

## 1. 极化子：只剩 Fe3O4 Verwey 电荷序

- 作业 594555，目录 `~/fe3o4_verwey/`（弛豫 + LEOS 单点两步，
  INCAR.leos 会在第二步自动换入）。
- 判读要点：B 位（八面体）16 个 Fe 应出现 ne=6/5 交替（Fe²⁺/Fe³⁺
  电荷序）；A 位（四面体）8 个 Fe 全 +3；元素平均 +2.667 是主指标。
  等价 B 位之间预期 R_site 低（位点简并，同 RESULTS 3.10 的方法学
  边界），用键长/磁矩佐证位点归属。
- 背景发现：历史 fe3o4 案例（data/v3/m4_wave1/fe3o4）的 MAGMOM 顺序
  写反了（磁矩给了 O 不是 Fe）——历史"电荷序摊平"部分源于此，写
  论文时值得注记。
- Li2O2 侧全链已闭环入库：`data/v3/m5a_li2o2/polaron/`
  （u0/u4/u6/u8/sp_pbe/sp_hse），RESULTS 3.10 已成文。

## 2. 测试集补齐：35 案例批次在跑

- 批次目录 `~/gap_batch`（29 个 Comas-Vilà 缺口 + 6 个 scdm_* 重跑）。
  QOS 限额由 nohup 守护自动补交（日志 `~/gap_submit.log`，
  batch_state.json 无 STAGED 后守护自杀）。
- 收割：`python3 ~/batch.py status --dest ~/gap_batch`，然后
  `harvest --out ...`；失败案例走 `retry`（梯子在 batch.py 头注释）。
- 评分前置：`pyleos/tests/reference_os.json` 需新增 gap 段参考值。
  参照 Comas-Vilà 论文 SI 的 S 表（他们的 EOS 结论）；正文解析产物在
  本会话 scratchpad 的 comas 目录（会话清理后需重新用 mineru-pdf 解析
  docs/refs/ 的 PDF）。
- scdm_* 六案例的目的是 escalate 触发实测：跑完 `batch.py escalate`
  验证 WCl4/5/6+MoCl2 被元素规则预升级、GaSe/CSe2 靠质量警告触发；
  快筛输出以 `*.scdm` 归档供对比。
- 尾巴（等用户放文件到 docs/refs/）：Wiley SI
  （10.1002/chem.202102679，Table S1 → MgH4/SrH4/CeH4/ThH4 的精确
  I4/mmm 高压相）与 ACS SI（10.1021/jacs.5b08162 → MgXe）。到位后
  转录 CIF 补算这 5 个，替换 4 个 APPROX 近似相。收集清单：本会话
  scratchpad `comas40_structures/manifest.tsv`（25 OK / 4 APPROX /
  1 MISSING），结构已同步到集群 `~/comas40_structures` 并转好
  `~/gap_src`。

## 3. M5-B 生产：轨迹在跑，管线已验证

- 生产 MD 作业 594556/594557/594558（每体系串 3 个脱锂点，50 ps，
  产物 `~/m5b_md/<sys>/<x>/prod/`：traj.extxyz + md.log +
  verdict.json）。验收 9/9 已通过并入库（data/m5b/acceptance）。
- 轨迹到手后铺抽帧 DFT：集群 `~/frames_to_dft.py`（mace env 里跑；
  `--system li2mno3|li2ruo3|nmc`；POTCAR 在 `~/potcars/`；默认取
  30 帧；生成 chain.sbatch 一个队列位串行跑完所有帧）。
- SCDM-only 对比组：同一轨迹再跑一份 `--scdm-only`（任务 4 的效率
  数据源，见下）。
- 时间序列汇总：仓库 `script/timeseries.py`（默认 tier2=arbiter，
  输出运输与贪心并排的元素平均曲线 + 稳定性摘要 + R_elem）。
- 冒烟链 594565（`~/m5b_dft_smoke/li2mno3_x025`，3 帧）已验证
  POSCAR/POTCAR/INCAR 物种顺序对齐。
- 预期管理（PLAN 已记）：O 空穴位点级可能 PBE 离域——元素平均是
  主指标；Li2RuO3 若暴露 4d 困难，如实写进适用性边界。

## 4. 效率基准：warm start 对照在跑

- 作业 594559（on）/ 594560（off），目录 `~/ws_bench/{on,off}`：
  96 原子 Li2MnO3 x0.25 的 VASP-AIMD 40 步，EOS_EVERY_N=5。
- 分析：对比两目录 wannier90*.wout 的迭代步数与耗时、print-out 的
  逐步 LOOP 计时；机制验证点——on 组从第 2 个分析帧起 win 文件里
  应出现 guiding_centres 块。
- SCDM-only vs 完整 MLWF 的大体系提速/一致率：用第 3 节的抽帧链
  `--scdm-only` 变体同帧双跑得到。
- 兑底方案（Li2O2 大超胞）未启用；M5-B 主路线通畅就不需要。

## 全局收尾清单（数据齐后）

1. 终版全量表：补齐案例并入后，以 `--tier2 arbiter` 默认重跑
   `script/full_bench.py`（含 reliability），替换
   data/v3/reliability_74.tsv 为 100+ 案例版。
2. RESULTS 未决事项节滚动更新，PLAN 里程碑表补 M5-B/补齐的完成态。
3. 两份 SI 到位后的 5 案例补算（见第 2 节尾巴）。
