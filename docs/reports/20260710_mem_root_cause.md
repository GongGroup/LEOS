> **归档注记** — 生成日期：2026-07-10；生成者：主会话派出的 OOM 取证 sub-agent；原路径：`/private/tmp/claude-501/-Users-mastreina-Desktop-EOS/04bfbd96-cc02-4d31-8159-643385b8af3d/scratchpad/mem_root_cause_report.md`（原件保留未动）。由成果收集会话拷贝入库，总索引见 [20260710_index.md](20260710_index.md)。

# LEOS/VASP 6.4.3 内存根因报告 — Ce18O36Co4CH4 每 rank ~10 GB 的构成与 OOM 根因

分析对象:VASP 6.4.3(LEOS v3 补丁)OpenACC GPU 版,sai-pod 16×V100 节点。
源码引用:`vasp.6.4.3.constrain/src/`(主干与原版 6.4.3 一致)、`patch/v3/eos.F`、`patch/leos_v3.patch`、sai-pod `~/apps/wannier90-3.1.0/src/`。
实测证据:sai-pod `~/LEOS/legacy_rerun/ce_{is,fs,ts}/` 的 OUTCAR/INCAR/wout/sbatch 与 sacct 记录(本报告已逐项核对)。
单位约定:GB = 10⁹ 字节;VASP 自报 kBytes = 10³ 字节;sacct MaxRSS 的 K = KiB。简报中 "56.2 GB/106.7 GB" 即 sacct 的 56,155,564 K / 106,703,532 K。

---

## ① 一页结论

**每 rank ~10 GB 的主体不是 bug,是三份全尺寸波函数数组 + 全 k 点投影算符,且这两块在 k 方向对 KPAR 完全不摊薄、在带方向只按 Nranks/KPAR 摊薄。**

1. **ALGO=ALL(IALGO≥50 → `INFO%LONESW=.TRUE.`,main.F:1998)让 main.F 在初始化时一次性分配 W、W_F、W_G 三份波函数**(main.F:2551-2554),每份 `CPTWFP(NRPLWV, NB_TOT/NB_PAR, NKPTS, ISPIN)`(wave.F:413)。k 维度是**全部 20 个 IBZ k 点**——`WDES%NKPTS=KPOINTS%NKPTS` 与 KPAR 无关(wave.F:123-136);KPAR 分组靠 `KPAR_SYNC_ALL` 事后广播(wave_mpi.F:1562),**KPAR 越大、每 rank 带数越多、波函数内存反而按 KPAR 倍增**。
2. 自报 10,165,693 kB 正是 ce_fs 成功重跑(**8 rank、KPAR=2**,NBANDS=332)OUTCAR 里 `DUMP_ALLOCATE`(main.F:2560,初始化完成点)打印的注册总量:wavefun 8.39 GB(=3×W)+ nonl-proj 1.37 GB + grid 0.30 GB + 杂项 0.11 GB。**其中 wavefun 与 nonl-proj(共 9.76 GB,96%)是每-rank 冗余或仅按 Nranks/KPAR 摊薄的块**,这就是 RSS 随 rank 数近线性增长(4r→8r:56→107 G-KiB)的原因——8 rank 作业同时把 KPAR 从 1 提到 2,带分组数没变(都是 4 组),每 rank 基底不变。
3. **56/107 GB 一句话**:每 rank 基底 ≈ 注册 10.2-10.3 GB + M_matrix 等 wannier 常驻 ~1.1 GB + 运行时 ~0.4-0.7 GB ≈ 11.5-12 GB;4 rank ≈ 46 GB、8 rank ≈ 92-98 GB,再加 rank0 在 SCDM 段的 **~11.5 GB 瞬态**(PSI + PSI_D_SCDM 两份全软网格×全带矩阵 + ZGEQP3 工作区,mlwf.F:3883/3910/3917)即得实测峰值 57.5 GB(4r)与 109.3 GB(8r,均为十进制 GB)。
4. **FS/TS 4-rank OOM 根因:INCAR 里 KPAR=2(IS 是 KPAR=1)**。4 rank + KPAR=2 → NB_PAR=2 → 每 rank 波函数翻倍:注册量 ≈ 18.9 GB/rank,4 rank 需 ~77-78 GB ≫ 62.5 GiB(67.1 GB)配额 → 在 main.F:2551-2554 三份 ALLOCW 清零(首次触页)期间被 cgroup 杀掉,发生在启动后 16-21 s(sacct:602730/602732/602742/602743/602749/602750 全部 OUT_OF_MEMORY,16-21 s)。**WAVECAR 假说证伪**:三个 INCAR 均 LWAVE=.FALSE.,ce_fs 目录的 WAVECAR 是 0 字节空壳;INWAV_FAST 读缓冲只有单带级(~MB,fileio.F:268/357),不存在全量读入缓冲;干净页缓存可回收、不会触发 cgroup OOM。
5. 次要发现:FS/TS 的 INCAR 少了 MAGMOM 行 → 默认 NBANDS 从 308(IS)涨到 332(推断机制,数值为 OUTCAR 实测),各块再放大 ~8%(M_matrix ~16%)。早先作业 594619(62.5G、18.5 分钟才 OOM)与 "4 rank KPAR=1 跑 332 带、死在 SCDM 段" 的模型峰值(~59-60 GB,贴着配额)吻合(作业身份为推断)。

---

## ② 内存预算表(Ce 参数代入)

参数(OUTCAR 实测):NB_TOT=332(FS/TS)/308(IS),NKPTS=20(IBZ),max PW=51,660(NRPLWV≈51,660-51,668,NPROD≈960),ISPIN=2,LMDIM=32,NTYP=5,NIONS=63,软网格 NGX×NGY×NGZ=90×90×126=1,020,600,细网格 180×180×252,NCORE=1(GPU 版强制),NB_PAR=Nranks/KPAR;wannier90:全网格 nk=32,num_wann=253/自旋,nnb=8,num_bands=332。复数 16 B。

### A. 初始化段(启动 ~20 s 内即分配,决定 OOM)

| # | 数组/段 | 尺寸公式(字节) | 4r,KPAR=1 (332带) | 8r,KPAR=2 (332带) | 4r,KPAR=2 (332带,OOM配置) | 冗余属性 | 源码位置 |
|---|---|---|---|---|---|---|---|
| 1 | W/W_F/W_G 的 CPTWFP ×3 | 3·16·NRPLWV·(NB_TOT/NB_PAR)·NKPTS·ISPIN | 8.23 GB | 8.23 GB | 16.47 GB | k 维全量冗余;带维按 NB_PAR=N/KPAR 摊 | wave.F:413;main.F:2551-2554;LONESW: main.F:1998,2552 |
| 2 | 同上 CPROJ ×3 | 3·16·NPROD·(NB_TOT/NB_PAR)·NKPTS·ISPIN | 0.15 GB | 0.15 GB | 0.31 GB | 同上 | wave.F:424;注册 wave.F:436 |
| | ↑ 合计 = OUTCAR "wavefun" | | 8,386,158 kB✓ | 8,386,158 kB✓ | 16,772,316 kB(推算) | | |
| 3 | NONL_S%QPROJ | 8·NGDIM·LMDIM·NTYPD·NKPTS = 8·51660·32·5·20 | 1.322 GB | 1.322 GB | 1.322 GB | **每 rank 全量冗余**(全 k、全类型;仅 NCORE>1 才切) | nonl.F:125-129 |
| 4 | NONL_S%CREXP | 16·NGDIM·NIONS | 0.052 GB | 0.052 GB | 0.052 GB | 每 rank 冗余 | nonl.F:124 |
| | ↑ 合计 = OUTCAR "nonl-proj" | | 1,374,569 kB✓(逐字节吻合) | 同左 | 同左 | | 注册 nonl.F:128 |
| 5 | grid(CHTOT/CHTOTL/DENCOR/CVTOT/CSTRF/CHDEN/SV) | 16·(11.5·MPLWV_C + 4·MPLWV_soft),按 COMM_KIN=N/KPAR 摊 | 0.296 GB | 0.296 GB | 0.591 GB | 分布(分母=N/KPAR) | main.F:2460-2467;注册 2494 |
| 6 | CHAM+CHF(LONESW 子空间矩阵) | 2·16·NB_TOT²·NKPTS·ISPIN | 0.141 GB | 0.141 GB | 0.141 GB | 每 rank 冗余,**未注册** | main.F:2555-2556 |
| 7 | one-center(CDIJ/CQIJ/CRHODE) | 3·16·LMDIM²·NIONS·NCDIJ | 6.2 MB | 6.2 MB | 6.2 MB | 冗余(63 离子全存) | main.F:2496-2498 |
| 8 | fftplans + base | 常数 | 0.103 GB | 0.103 GB | ~0.10-0.18 GB | 常数 | fftmpi_map.F:69-88;main.F:491 |
| 9 | CELTOT/FERTOT/AUXTOT 等 | ~NB_TOT·NKPTS·ISPIN·32 | ~1 MB | ~1 MB | ~1 MB | 冗余(可忽略) | wave.F:439-443 |
| 10 | INWAV_FAST 读缓冲(ISTART=1) | ~5·MAXVAL(NPLWKP_TOT)·8-16 | 3-4 MB | 3-4 MB | 3-4 MB | 瞬态;**无全量缓冲** | fileio.F:268,357 |
| **Σ 注册值(=自报)** | | | **9.56 GB (308带,实测) / 10.17 GB (332带,推算)** | **10.17 GB✓实测** | **≈18.9 GB(推算)** | | DUMP_ALLOCATE main.F:2560 |

### B. SCDM / wannier90 / LEOS 段(SCF 收敛后)

| # | 数组/段 | 尺寸公式 | 数值 (332带) | 位置/属性 | 源码 |
|---|---|---|---|---|---|
| 11 | MLWF%M_matrix | 16·NB²·nnb·nk_full·ISPIN = 16·332²·8·32·2 | **0.90 GB**(308带:0.78) | **每 rank 冗余**,常驻至运行结束 | mlwf.F:1224 |
| 12 | A_MATRIX_TMP → A_matrix | 16·NB_TOT·num_wann·nk_full·ISPIN ×2(瞬态并存) | 86 MB ×2 → 86 MB 常驻 | 每 rank | mlwf.F:3838-3839, 573 |
| 13 | U_matrix + U_matrix_opt + lwindow | 16·(NW²+NB·NW)·32·2 | 65.5+86 MB | 每 rank(EOS_ANALYZE 需要) | mlwf.F:1863-1865 |
| 14 | **SCDM: PSI** | 16·(NGX·NGY·NGZ)·NB_TOT = 16·1,020,600·332 | **5.42 GB**(308带:5.03) | **仅 rank0(NODE_ME==1)**,瞬态,跨两个自旋循环存活 | mlwf.F:3842,3883 |
| 15 | **SCDM: PSI_D_SCDM(转置副本)** | 同上 | **5.42 GB** | rank0,QRCP 期间与 PSI 并存 | mlwf.F:3910-3914 |
| 16 | SCDM: ZGEQP3 工作区 | LWORK≈(N+1)·NB_blk·16,N=1,020,600 + RWORK 2N·8 | 0.5-1.0 GB + 16 MB(推断,依 LAPACK 块大小) | rank0,瞬态;**串行 CPU LAPACK ZGEQP3** | mathtools.F:266-302(RRQR=xGEQP3);调用 mlwf.F:3887,3917 |
| 17 | SCDM: SVD(锚点小矩阵) | 332×253 'S' 模式 | ~3 MB | rank0 | mathtools.F:105-193;mlwf.F:3891,3959 |
| 18 | wannier_run(w90 3.1.0 库) | m_matrix_orig + m_matrix_orig_local 各 16·NB²·8·32;m_matrix(NW) 16·253²·8·32 | ~0.90+0.90+0.26 ≈ **2.0-2.3 GB** | **仅 rank0,串行**;自旋两通道顺序两次调用 | w90 overlap.F90:65,68;mlwf.F:1875-1904 |
| 19 | MMN 期 PEAD overlap 瞬态 | S(NB²) + 波函数旋转/收集缓冲 | ~0.5 GB(推断) | 每 rank,瞬态 | mlwf.F:1221,1250(PEAD_CALC_OVERLAP) |
| 20 | GET_WAVE_FUNCTIONS 缓冲 | WVFN 全软网格 + NB_PAR 条带 CR | 16 MB + ~65 MB | 每 rank,瞬态 | mlwf.F:3141-3142,3136-3139 |
| 21 | **eos.F(LEOS 本体)** | UTOT 16·332·253;CW 16·NPRO(≈960)·253;PW 8·253·63·2;WF_* | **1.3 + 3.9 + 0.26 MB + ~40 KB = <6 MB** | 每 rank,用后即 DEALLOCATE ✓无辜 | eos.F:843-844(UTOT/CW)、497(PW)、416-421、764-765 |
| 22 | 未注册运行时(CUDA ctx/UCX/heap/libs) | — | ~0.4-0.7 GB/rank(对账余项,推断) | 每 rank | — |

k 点维度说明:A/B 表中 SCF 数组用 NKPTS=20(IBZ);wannier 系数组用 nk_full=32(`MLWF%NKPTS=KPOINTS_FULL%NKPTS`,mlwf.F:964-966,NEQUIV 映射回 IBZ)。MMN 的 k 循环仅按 COMM_KINTER(=KPAR 组)切分(mlwf.F:1241):KPAR=1 时每 rank 冗余计算全部 32×8 个重叠(时间冗余,非内存)。

---

## ③ 与实测对账

| 实测 | 值 | 模型 | 误差 |
|---|---|---|---|
| ce_fs OUTCAR 自报(8r,KPAR=2,332带) | 10,165,693 kB | wavefun 8,386,158(=3×W,逐 kB 吻合)+ nonl 1,374,569(**逐字节**:8·51660·32·5·20 + 16·51660·63 = 1,374,569,280 B)+ grid 295,617 + fft 73,156 + base 30,000 + 1c 6,193 = 10,165,693 | 0 |
| ce_is OUTCAR 自报(4r,KPAR=1,308带) | 9,559,465 kB | 同构,NBANDS_loc=77:wavefun 7,779,930 | 0 |
| lg-ce_is MaxRSS(602731,62.5G 过) | 56,155,564 KiB = 57.5 GB | 4×(9.56+0.14+0.78+0.24+~0.5 运行时) + rank0 SCDM(5.03×2+0.6)+w90 段 ≈ 56-58 GB;峰值落在 SCDM 窗口 | <2% |
| lg-ce_fs MaxRSS(604671,8r,125G 过) | 106,703,532 KiB = 109.3 GB | 8×(10.17+0.14+0.90+0.24+~0.5) + rank0 SCDM(5.42×2+0.6) ≈ 107-109 GB | <2% |
| lg-ce_fs/ts 4r OOM(602730/32/42/43/49/50) | 62.5 GiB 配额,**16-21 s** 全灭 | KPAR=2→NB_PAR=2:注册 ≈18.9 GB/rank,4 rank 需 ~77-78 GB;死于三份 ALLOCW 清零(wave.F:446 `W%CPTWFP=0` 立即触页)——16-21 s = MPI/读 POTCAR/FFT 规划/SETUP_FOCK 后到达 main.F:2551 的时刻 | 定性+定量吻合 |
| IS 过 / FS-TS 挂 | — | **KPAR=1 vs KPAR=2**(INCAR diff 实证)+ NBANDS 308 vs 332(FS/TS 缺 MAGMOM 行,推断机制) | — |
| WAVECAR 假说 | 证伪 | ce_fs WAVECAR=0 B(LWAVE=F);INWAV_FAST 仅单带缓冲(fileio.F:357);页缓存可回收不致 cgroup OOM | — |
| 594619(62.5G,18.5 min OOM,采样 35.9 GB) | — | 推断为 4r/KPAR=1/332带 的 FS 早期尝试:init ~44 GB 通过,SCF ~15 min 后进 MMN+SCDM,峰 ~59-60 GB 贴配额被杀 | 推断 |

RSS 线性增长的解释:4r→8r 时把 KPAR 同步从 1 提到 2,NB_PAR 恒为 4,于是"每 rank 注册基底"不变(10.2 GB),总量随 rank 数近线性;若 8r 保持 KPAR=1,每 rank wavefun 应减半(预测总量 ~76 GB 而非 107)。**KPAR 是伪装成并行参数的内存倍增器**。

---

## ④ SCDM / wannier90 / LEOS 逐段分析

**SCDM(mlwf.F:3799-4006,stock;LEOS 经 `IF (LEOS) LSCDM=.TRUE.` 强制开启,patch main.F 钩子)**
- 仅 `COMM%NODE_ME==1` 干活(3842),其余 rank 只陪跑 FFT/gather 并等广播(3865-3880,4005)。
- rank0 分配 `PSI(NP·NRSPINORS, NB_TOT)`,NP=`GRID%RL%NP`;NCORE=1 时 RL%NP=全软网格 1,020,600 点(GET_WAVE_FUNCTIONS 内 `WVFN(NGX·NGY·NGZ)` 经 `M_sum_g(COMM_INB)` 得全网格带,mlwf.F:3141-3186)→ **PSI=5.42 GB**;再显式转置成 `PSI_D_SCDM(NB_TOT, NP)` 又一份 5.42 GB(3910-3914),两者在 QRCP 期间并存。
- QRCP=LAPACK **ZGEQP3**(mathtools.F RRQR_DECOMPOSER,查询式 LWORK≈(N+1)·块宽),**单 rank 串行 CPU** 跑 332×1,020,600 矩阵;每个自旋通道各一次(锚点=Γ);非锚点 k 只做小矩阵(NB_TOT×NUM_WANN)SVD(3950-4000)。
- 内存后果:rank0 瞬态 +11.3-11.9 GB;时间后果:全部 2 自旋×32 k×332 带的实空间重建(CPU FFT)+ 一次大 QRCP,串行段拉长 wannier 阶段。

**wannier90 输入与库调用(mlwf.F)**
- M_mn:`M_matrix(332,332,nnb=8,nk=32,ISPIN=2)`=0.90 GB **每 rank 分配**(1224),k 循环只按 KPAR 组切(1241),用 `M_sum_z8(COMM_KINTER)` 归并(1272)→ KPAR=1 时计算与存储双冗余。overlap 经 PEAD_CALC_OVERLAP(1250)。
- A_mn:SCDM 路径得 `T(NB_TOT,253,32,2)` 每 rank 86 MB,广播自 rank0(3838,4005);复制进 `A_matrix` 后释放 TMP(573-584)。
- eig:每 rank 小数组;.mmn/.amn/.eig 文件由 IONODE 写(WRITE_MMN 1299)。
- `wannier_run` **仅 rank0**(`IF (IO%IU6>=0)`,1875),ISPIN=2 时**顺序两次**(1886-1902,seedname wannier90.1/wannier90.2);U/U_opt/centres 广播回所有 rank(1908-1913)。w90 3.1.0 库内部再持有 m_matrix_orig(_local) 2×0.90 GB + num_wann 级 m_matrix 0.26 GB(overlap.F90:65,68)→ rank0 w90 窗口 +2.0-2.3 GB(小于 SCDM 窗口,峰值由 SCDM 决定)。
- `LWANNIER90_AUTO_WINDOW`(stock ESTIMATE_AUTO_WINDOW,mlwf.F:884;patch 增强了带隙中点选择)只扫 CELTOT,内存可忽略。

**LEOS(eos.F,787 行)**
- EOS_ANALYZE/EOS_PROJ_WEIGHTS 全部数组:NE_ASSIGNED/OS_*(63)、WF_OCC/WF_I/WF_E(253×2)、PW(253×63×2)、UTOT(332×253 复)1.3 MB、CW(NPRO≈960×253 复)3.9 MB —— **合计 <6 MB/rank,用后 DEALLOCATE(eos.F:764-765),量级确认:无辜**。
- 所有 rank 重复同样计算(U 收缩 253²×32×2 ≈ 0.4 GFLOP,秒级),结尾 `M_sum_d(COMM, PW)`(916)要求全员参与——设计一致,无内存问题。
- 患处仅在于:LEOS 默认强制 LSCDM + LWANNIER90(EOS_READER,eos.F:297-301;main.F 钩子 `IF (LEOS) LSCDM=.TRUE.`),即**每个 LEOS 作业都自动背上 ④ 的 rank0 SCDM 峰值与每 rank 0.9 GB 的 M_matrix**。

---

## ⑤ 修复候选

### (a) 补丁范围内(eos.F、patch 已触及的 mlwf.F/main.F 段)

| 候选 | 做法 | 预计收益 | 风险 |
|---|---|---|---|
| A1. SCDM 免双缓冲 | CALC_WANNIER_TRAFO_SCDM 不先存 PSI 再转置:逐带取回 WVFN 直接写入 PSI_D_SCDM(IB,:)(共轭+cutoff 就地做);锚点 pivot 行事后按需重建(仅 NUM_WANN 行) | **rank0 −5.0/−5.4 GB**(308/332 带);IS 4r 峰值 57.5→~52 GB,给 62.5G 配额留出安全边 | 中:改 stock 例程(补丁本就 patch mlwf.F);数值不变,需一次回归(EOSCAR 对比) |
| A2. M_matrix 归约到 root | mlwf.F:1224 只在 IONODE 分配全量;各 rank 算完本地 k 后 reduce-to-root(替代现 allreduce 语义 M_sum_z8:1272);w90 只在 rank0 用它 | 非 root **−0.90 GB/rank**(8r 省 6.3 GB) | 低-中:WRITE_MMN 与 wannier_run 都在 IONODE,天然兼容;注意 EOS/重启路径没人再读 M |
| A3. 用后及早释放 | wannier_run 返回并广播 U 后立刻 `DEALLOCATE(MLWF%M_matrix, MLWF%A_matrix)`(EOS_ANALYZE 只用 U/U_opt/lwindow/centres/spreads;EOS_EVERY_N 重跑前 MLWF_FREE 已会重建) | 每 rank 常驻 −1.0 GB(w90 之后直至作业结束) | 低 |
| A4. QRCP 网格子采样(SCDM 文献常规) | 锚点 QRCP 在按 |ψ|² 加权抽样的 ~10-20% 网格点上做(pivots 只需空间代表点) | rank0 SCDM 矩阵 5.4→0.5-1.1 GB,QRCP 时间 −80%+ | 中-高:科研代码改动,需验证 pivots/初猜质量(可用 EOSCAR/spread 回归) |
| A5. eos.F | 无需改(MB 级已确认);可在文档注明 LEOS 默认强制 LSCDM 的内存代价 | — | — |

### (b) VASP 主干 / wannier90 上游(只诊断 + 配置绕行)

| 问题(上游) | 诊断 | 配置级绕行 | 预计收益/风险 |
|---|---|---|---|
| B1. KPAR 不切 k 维存储 | wave.F:413 全 NKPTS 分配 + KPAR_SYNC_ALL(wave_mpi.F:1562);每 rank wavefun ∝ KPAR/Nranks | **铁律:rank 数不足时 KPAR=1**。本体系 4 GPU/62.5G:必须 KPAR=1;8 GPU/125G:KPAR=2 可用(现行成功配置) | 4r 下 KPAR 2→1:每 rank −8.5 GB(注册 18.9→10.2),OOM→可跑 |
| B2. ALGO=ALL 三份 W(LONESW) | main.F:2551-2556 | 改 ALGO=Normal(Davidson)或 Fast:只剩 1 份 W+CPROJ | **−5.6 GB/rank**(8r/KPAR2)、−11.2(4r/KPAR2);风险:LDAU+真空表面收敛性,需在 IS 上试收敛等价性(NELM/AMIX 已调过,中风险) |
| B3. NBANDS 膨胀 | FS/TS 缺 MAGMOM → 默认 NBANDS 332(IS 308)(机制推断,数值实测) | INCAR 显式 `NBANDS=308`(或补 MAGMOM);num_wann=253 仍有 55 空带 | −0.6 GB/rank(wavefun)−0.13(M_matrix);零风险 |
| B4. QPROJ 全 k 常驻 | nonl.F:125(NKPTS 维,NCORE=1 不切) | 唯一杠杆是 k 网格密度(如 4×4×1 或 Γ-偏移减 NKPTS);LREAL=Auto 可换成 RPROJ 实空间(63 原子体系 LREAL=.FALSE. 是精度选择) | k 减半:QPROJ −0.66 GB、wavefun −4.2 GB/rank;风险:收敛性/能量一致性需验证 |
| B5. w90 库串行 rank0 + SCDM rank0 峰 | mlwf.F:1875;3842 | 给 rank0 留余量:峰值预算 = 每rank基底 + **12-13 GB**(SCDM+w90);62.5G/4r 场景等价于"其余 3 rank 只能用 ~50 GB"。或 `LWANNIER90_RUN=.FALSE.` 产出 mmn/amn 后独立跑 wannier90.x(错峰,不省 SCDM) | 规划级,无代码风险 |
| B6. sai-pod 配额与 GPU 数联动 | 62.5G=4 GPU、125G=8 GPU(sacct ReqMem 实证) | 内存吃紧时宁可申请 8 GPU 只为 125G(现行做法);或 `--mem` 显式提额(若站点允许) | 调度成本换稳定性 |
| B7. PREC 软网格 | SCDM PSI ∝ NGX·NGY·NGZ(90×90×126) | PREC=Normal 软网格约 −40-55% → PSI 5.4→~2.5-3.2 GB | 精度权衡(力/应力),表面反应能需复核,不作首选 |

**推荐组合(立即可用,零代码)**:4 GPU 档 → `KPAR=1 + NBANDS=308 + MAGMOM 补回`(峰值 ~52-58 GB,配 A1 后更稳);8 GPU 档维持 `KPAR=2 + 125G`。**代码上最划算**:A1(rank0 −5 GB)+ A2/A3(每 rank −1 GB 常驻),合计把 4r/62.5G 档从"贴线"变"从容"。

---

### 附:证据索引
- OUTCAR 内存块:`sai-pod:~/LEOS/legacy_rerun/ce_is/OUTCAR`(9,559,465 kB)与 `ce_fs/OUTCAR`(10,165,693 kB;`distrk: each k-point on 4 cores, 2 groups`)。
- sacct:602730/32/42/43/49/50(lg-ce_fs/ts,62.5G,16-21 s,OUT_OF_MEMORY);602731(lg-ce_is,56,155,564K,COMPLETED);604671/604672(125G,106,703,532K/100,193,604K,COMPLETED);594619(62.5G,18:28,OOM,MaxRSS 采样 35.1 GB)。
- INCAR diff(仓库 `data/v3/legacy_rerun/`):ce_is `KPAR=1`+MAGMOM;ce_fs/ce_ts `KPAR=2`、无 MAGMOM。
- wannier90.1.wout(ce_fs):4×4×2=32 点、num_wann=253、input Bloch=332、b_k 8 支。
- w90 3.1.0 源:`sai-pod:~/apps/wannier90-3.1.0/src/overlap.F90:65,68`。
