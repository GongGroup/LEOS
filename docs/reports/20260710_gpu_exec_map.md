> **归档注记** — 生成日期：2026-07-10；生成者：主会话派出的 CPU/GPU 执行位置地图 sub-agent；原路径：`/private/tmp/claude-501/-Users-mastreina-Desktop-EOS/04bfbd96-cc02-4d31-8159-643385b8af3d/scratchpad/gpu_exec_map.md`（原件保留未动）。由成果收集会话拷贝入库，总索引见 [20260710_index.md](20260710_index.md)。

# LEOS 单点任务 CPU/GPU 执行位置地图(源码 + sai-pod 实测)

**范围**:VASP 6.4.3 + LEOS 补丁,OpenACC GPU 版(nvhpc 24.3 + openmpi 5.0.8,`-D_OPENACC -DVASP2WANNIER90`),sai-pod 16×V100 节点。
**主例**:`~/LEOS/legacy_rerun/ce_fs/`(Ce18O36Co4CH4 类体系,NBANDS=332,ISPIN=2,IBZ 20 k 点 / 全 BZ 32 k 点,波函数网格 90×90×126 ≈ 1.02M 点,num_wann=253,ALGO=ALL,8 rank × 8 GPU,KPAR=2)。
**副例**:`~/LEOS/ws_bench/{on,off}`(Li2MnO3 MD 40 步,EOS_EVERY_N=5,NBANDS=num_wann=352,4 rank × 4 GPU)。
**源码基线**:本地 `vasp.6.4.3.constrain/src`(与生产树 `~/apps/vasp.6.4.3_leos_v2/src` 交叉核对:两边 `mlwf.F`、`eos.F` 的 `!$ACC` 均为 0);LEOS 补丁 `patch/v3/eos.F` + `leos_v3.patch`;wannier90-3.1.0 源码在 sai-pod `~/apps/wannier90-3.1.0/`。
标注约定:**[实证]** = 源码行或服务器日志直接支撑;**[推断]** = 由机制推导,未单独计时。

---

## ① 一页结论

1. **已在 GPU 的只有 DFT SCF 本身**(ALGO=ALL 全带 CG:哈密顿作用、FFT、非局域投影、子空间旋转、最终对角化)。Ce 单点实测 SCF = **13 120 s,占总时 77.2%**(554 步 × 23.7 s/步)。
2. **SCF 结束的一瞬间,GPU 的参与就归零**:`electron_all.F` 出口把波函数一次性拷回 host 并**删除 device 副本**(`UPDATE SELF` + `ACC_DELETE_TYPED_VAR(W)`)。此后 SCDM、M_mn、A_mn、wannier90 库、LEOS 分析、全部收尾 I/O——**全部 CPU**,期间 8 张 V100 完全闲置(Ce 单点闲置 3 775 s ≈ 8.4 GPU·h)。
3. **Wannier 管线(CPU 段)= 3 775 s,占 22.2%**,内部分解:
   - MMN(重叠矩阵)≈ 482 s(8 个 CPU 核并行);
   - SCDM(实空间波函数组装 + QRCP + SVD)≈ 470 s(组装全员参与,QRCP/SVD 只在 1 个 rank 的 1 个核上);
   - **wannier90 库 = 2 817 s(占管线 75%)**:单核、串行、两自旋顺序跑;其中自旋↓的 disentangle 打满 1000 次迭代耗 1 743 s,是单点里最大的单块 CPU 时间;
   - LEOS 分析段(eos.F)≈ 4 s,可忽略。
4. **wannier90 是"最差配置"在跑**:库模式只有 master rank 调用,链接的是 nvhpc 自带**参考版 LAPACK/BLAS**(`-llapack -lblas`,无线程、无优化 GEMM),`FCOPTS=-O2` 无并行。也就是说 128 核 + 16 GPU 的节点上,这 2 817 s 只用了 1 个 CPU 核。
5. **轨迹/逐帧模式下主次颠倒**(m5b 实测):分析帧上 wannier90 占帧内电子步耗时 **~81%**(479 s / 590 s),全轨迹占 **63-64%**;SCF 只占 27%。若目标是 MD 逐帧氧化态,GPU 化/并行化 Wannier 管线的收益远大于继续优化 SCF。
6. host↔device 大宗搬运只有两处:ELMIN 入口整套 W 上卡、出口整套 W 下卡(初始 ORTHCH 自带一次小的上/下卡)。Wannier 段没有任何 H↔D 流量,但有大量 **host 侧 MPI 流量**(每带全网格广播/规约 ×21 248 次、903 MB 的 M_matrix allreduce)。

---

## ② 阶段地图表(Ce 单点,总墙钟 17 000 s)

时间边界来自 OUTCAR 计时 + 输出文件 mtime + wout 内部计时,三者交叉一致 **[实证]**。

| # | 阶段 | 源码位置(file:line) | 执行位置 | Ce 实测秒数(占比) | host↔device / 通信搬运点 |
|---|------|------------------|----------|------------------|--------------------------|
| 0 | 启动:读 POTCAR/INCAR、WAVECAR(本次 0 字节→未读)、随机初猜 WFINIT、PROALL、初始 ORTHCH、初始电荷 | main.F:2862-2971,3019;wave.F:1310/1334;choleski2.F:66 | CPU 为主;**ORTHCH 是自包含 GPU 区域** | ~53 s(0.3%) | ORTHCH 内:H→D `ACC_COPYIN_TYPED_VAR(W)`(choleski2.F:70)→ GPU Cholesky → D→H `ACC_COPYOUT_TYPED_VAR(W)`(choleski2.F:372) |
| 1 | SCF,ALGO=ALL 全带 CG,554 步 | electron_all.F:241(PUSH GPU)+ rot.F/hamil.F/nonl(r).F/fftmpiw.F/subrot.F 等 | **GPU**(波函数常驻 device) | **13 120 s(77.2%)**,23.7 s/步 | 入口 H→D:electron_all.F:246 `ACC_COPYIN_TYPED_VAR(W)` 及 NONL/GRID/P 等(241-258);迭代内小量 D↔H(本征值/占据/混合)[推断] |
| 2 | 最终对角化(final diagonalization)+ 波函数落地 | electron_all.F:630-633(UPDATE SELF→打印)、subrot.F:163 | GPU(EDDIAG 自开) | 含在 LOOP+ 缝隙,~秒级 | **D→H 全量**:electron_all.F:630、906 `UPDATE SELF(W%CW,W%CPROJ,W%CELTOT)`;**device 副本删除**:electron_all.F:929-930 `ACC_DELETE_TYPED_VAR(W/WDES)`;POP GPU:933。**此后再无 GPU 参与** |
| 3 | wannier_setup(库模式,k 网格/近邻表) | mlwf.F:1064-1072(仅 IU6≥0 的 master rank) | CPU,单 rank | <1 s(wannier90.wout 06:46:01) | 结果 bcast 给各 rank |
| 4 | M_mn 重叠矩阵计算 | mlwf.F:1194-1275 → pead.F:1530 → CALC_OVERLAP_ pead.F:4151-4356(逐带 W1_ROTATE + FFTWAV_W1 + ORTH2/zgemm) | **CPU**,k 点按 KPAR=2 分组、组内 4 rank 带并行(mlwf.F:1240-1242)→ 共 8 核 | ~465 s(2.7%) | 末尾 `M_sum_z8` 对 M_matrix(332²×8×32×2 复数 ≈ **903 MB**)全归约(mlwf.F:1272) |
| 5 | 写 .mmn(1.04 GB×2)+ .eig | mlwf.F:1280-1344、1727-1771 | CPU,rank0 格式化写 | ~17 s | 磁盘 |
| 6 | SCDM:实空间波函数组装 + QRCP + SVD → A_mn | mlwf.F:3799-4006(QRCP 3917=xGEQP3;SVD 3959/3995=xGESDD;A=U·V† 3964/4000);组装 GET_WAVE_FUNCTIONS mlwf.F:3073 → W1_GATHER_W1 wave_high.F:4000-4022(逐带 FFTWAV_W1) | **CPU**;组装阶段全 8 rank 参与 FFT/通信,**QRCP/SVD 仅 NODE_ME==1 的单核**(mlwf.F:3842) | ~470 s(2.8%)(含写 .amn 139 MB×2 约 2 s) | 每带:全网格 CR 广播(wave_high.F:4030,~16 MB/次)+ `M_sum_g` 全网格 allreduce(mlwf.F:3186);332 带×32k×2 自旋 = **21 248 次**;结果 T 矩阵 bcast(mlwf.F:4005) |
| 7 | wannier90 库:disentangle + wannierise,自旋↑ | mlwf.F:1875-1904(仅 master;ISPIN=2 时 DO ISP 顺序,mlwf.F:1886);w90 源 disentangle.F90 / wannierise.F90 | **CPU 单核单 rank**(参考 BLAS/LAPACK) | **716 s(4.2%)** = dis 374.2 + wann 342.2(wannier90.1.wout) | 其余 7 rank 阻塞在 `M_sum_z8(U_matrix…)`(mlwf.F:1908-1913);8 GPU 闲置 |
| 8 | 同上,自旋↓ | 同上 | 同上 | **2 101 s(12.4%)** = dis 1 743.2(打满 1000 迭代)+ wann 357.3(wannier90.2.wout) | 同上 |
| 9 | LEOS 分析段(占据统计、中心配对、写 EOSCAR) | patch/v3/eos.F:195+(EOS_ANALYZE);挂钩 main.F(leos_v3.patch @@3329) | **CPU**(eos.F 零 ACC,无 BLAS) | **~4 s(0.02%)**(EOSCAR mtime − w90 结束 = 4 s) | 输入 wann_centres/U_matrix 已被 M_sum 同步到全 rank |
| 10 | 收尾 I/O:PROCAR/DOSCAR/vasprun.xml/AECCAR2 等 | main.F(LOOP+ 计时器内,main.F:3200/4916) | CPU | ~61 s(0.4%) | 磁盘 |

**分段合计**:SCF(GPU)13 120 s;Wannier 管线(CPU,阶段 3-9,06:46:01→07:48:56)= **3 775 s = 22.2%**;其中 w90 库 2 817 s = 管线的 74.6%。
交叉验证:LOOP+ = 16 947 s(main.F:3200-4916 覆盖阶段 1-10)≈ 13 120 + 3 775 + 收尾,自洽 **[实证]**。

---

## ③ 各问题的证据

### Q1 SCF 初始波函数:ISTART=0 与 ISTART=1 各步在哪执行

- **随机初猜(ISTART=0/INIWAV=1)**:`main.F:2862-2864` → `WFINIT`(wave.F:1310)→ `WFINIT_DEFAULT`(wave.F:1334):`RANE()` 逐元素在 **host** 填充 `W%CPTWFP`,该子程序无任何 ACC 指令 **[实证]**。
- **初始正交化**:`main.F:2971 CALL ORTHCH`。ORTHCH 是**自包含 GPU 区域**:choleski2.F:66 `PUSH_ACC_EXEC_ON(.TRUE.)`,:70 把整套 W 拷上 device,GPU 上做 Cholesky/子空间重叠,choleski2.F:372 `ACC_COPYOUT_TYPED_VAR(W)` 拷回并删除 device 副本 **[实证]**。
- **初始电荷(来自波函数)**:`main.F:3019 SET_CHARGE` 的 GPU 包装被 `#ifdef _OPENACC_TST` 包住(main.F:3008-3018);生产 makefile.include 只定义 `-D_OPENACC -DVASP2WANNIER90`,**无 `_OPENACC_TST`** → 该块编译剔除,初始电荷在 **CPU** **[实证]**。
- **初始子空间迭代/延迟步(NELMDL)**:发生在 ELMIN_ALL 内部,已处于 GPU 区域(electron_all.F:241 PUSH TRUE + :246 COPYIN(W))→ **GPU** **[实证]**。
- **ISTART=1(读 WAVECAR)**:`main.F:2585 INWAV_FAST`(fileio.F,`!$ACC`=0)→ **CPU 读盘**;之后同样经 PROALL/ORTHCH,最终在 ELMIN 入口统一上卡 **[实证]**。
- **本次 Ce 实测的真实路径**:INCAR 写 `ISTART=1`,但目录里 WAVECAR 为 0 字节;print-out 打印 `WAVECAR not read` 与 `WARNING: random wavefunctions but no delay for mixing` → **实际走了随机初猜**,554 步 CG 从零收敛 **[实证]**。
- **结论**:"生成初始波函数"在 GPU 版里 = host 随机数(秒级)+ GPU 正交化 + GPU 内延迟迭代;整个初始化 ~53 s(0.3%),不构成性能问题。真正的"还在 CPU"不在这里(见 ④)。

### Q2 SCDM 段(mlwf.F)

- **实空间波函数矩阵怎么来**:`CALC_WANNIER_TRAFO_SCDM`(mlwf.F:3799)对每个自旋、每个 k(锚点 Γ + 其余 31 个)调 `GET_WAVE_FUNCTIONS`(mlwf.F:3073):`W1_GATHER_W1`(wave_high.F:4000)逐带做 `FFTWAV_W1`(wave_high.F:4021)把平面波系数 FFT 到实空间,再把全网格 CR 广播给 COMM_INTER 全体(wave_high.F:4030),重排成 NGX·NGY·NGZ 顺序后 `M_sum_g` 全网格规约(mlwf.F:3186)。`NSTRIP=1`(mlwf.F:3102)→ 一次一条带。
- **FFT 在 GPU 吗?不在。** 机制:VASP 的 ACC 指令都带 `IF(ACC_EXEC_ON)`/`ACC_IS_PRESENT` 守卫;`ACC_EXEC_ON` 是栈式开关,0 层默认 `.FALSE.`(openacc.F:41-42),只有显式 `PUSH_ACC_EXEC_ON(.TRUE.)` 的区域才开。全源码扫描:自开 GPU 的例程清单里没有 mlwf/pead 重叠链上的任何一个(pead.F 仅 PEAD_ELMIN 自开:2131/2508);且 SCF 出口已删除 W 的 device 副本(electron_all.F:929)。因此 MLWF 段里 FFTWAV_W1→FFTWAV 走 **CPU FFTW 分支** **[实证(机制)+推断(未逐条跑运行时验证)]**。既然从未上卡,也就**不存在"算完拷回 host"**——数据从 SCF 出口那一次 UPDATE SELF 之后就一直住在 host。
- **QRCP 用什么、谁执行**:`RANK_REVEAL_QR%CALC`(mlwf.F:3917),类型 `RRQR_DECOMPOSER`(mathtools.F)封装 LAPACK **xGEQP3**(错误信息直接写明,mlwf.F:3889/3919;复数版本即 zgeqp3)。**只在 `W%WDES%COMM%NODE_ME==1` 那一个 rank 执行**(mlwf.F:3842 `DO_WORK`),其余 rank 只陪跑波函数组装然后在 bcast 处等(mlwf.F:3865-3880)。矩阵尺寸 NB_TOT×(NP·NRSPINORS) = 332 × ~1.02M,只在锚点 k 做一次/自旋。执行位置:**CPU 单核**(mathtools 无 ACC)**[实证]**。
- **SCDM 后 A_mn 构造**:选出 pivot 列 → 每个 k 上取 332×253 小矩阵,`SVD%CALC` = **xGESDD**(mlwf.F:3959/3995),A(k)=U·V†(mlwf.F:3964/4000,MATMUL);全部 k 做完后 `M_bcast_z_from` 广播 T(mlwf.F:4005),回 MLWF_MAIN 重排进 `MLWF%A_matrix`(mlwf.F:571-584),rank0 写 .amn(mlwf.F:1674-1699)。全 CPU **[实证]**。
- **实测**:.eig/.mmn 写完(06:54:03)→ .amn 写完(07:01:55)= **~470 s**,其中组装(21 248 次 FFT+广播+规约)与 QRCP(332×1M 列选主元)为主 **[实证(区间)+推断(区间内部拆分)]**。
- 顺带(内存,归内存报告管):rank1 上 `PSI(NP·NRSPINORS, NB_TOT)` 与 `PSI_D_SCDM` 按满网格×全带分配,复数 16 B → 单块 ~5.4 GB 量级。

### Q3 M_mn / A_mn / 重叠构造在哪算、数据何时离开 device

- **M_mn**:`MLWF_WANNIER90_MMN`(mlwf.F:1194):对 32 k × 8 近邻 × 2 自旋逐对调 `PEAD_CALC_OVERLAP`(pead.F:1530,薄包装)→ `CALC_OVERLAP_`(pead.F:4151-4356):逐带 `W1_ROTATE`(IBZ→全 BZ 对称旋转)+ `FFTWAV_W1` + 相位/增强电荷(`OVERL_AND_APPLY_PHASE_SHIFT`,PAW Qij)+ `ORTH2`(dfast.F,zgemm 型重叠)。**全 CPU**(同 Q2 机制;dfast.F 不自开 GPU)。并行度:k 循环按 `COMM_KINTER`(KPAR=2)瓜分(mlwf.F:1240-1242),组内 4 rank 按带分工 → 共 8 个 CPU 核 **[实证]**。
- **数据离开 device 的时刻**:不在 MLWF 内——**早在 ELMIN_ALL 出口**(electron_all.F:906 UPDATE SELF;:929 删除)就完成了。MLWF 全程只读 host 数组 **[实证]**。
- **A_mn**:LSCDM 路径下 A_mn 就是 SCDM 的 T(Q2);非 SCDM 路径才走 `MLWF_CALC_AMN`(mlwf.F:1553)。
- **实测**:MMN 计算 ~465 s + 写盘 ~17 s **[实证]**;末尾 903 MB 的 `M_sum_z8` allreduce 含在其中 **[实证(代码)+推断(耗时占比)]**。

### Q4 wannier90 库调用

- **串行、单 rank、两自旋顺序**:`wannier_run` 只在 `IO%IU6>=0` 的 master rank 上调用(mlwf.F:1875);ISPIN=2 时 `DO ISP=1,2` 顺序执行(mlwf.F:1886-1902),打印顺序与 print-out 一致;其余 rank 阻塞在结果规约(mlwf.F:1908-1913)**[实证]**。库本身以 `libwannier.a` 静态链入(makefile.include `LLIBS += -lwannier`),构建:`F90=nvfortran, FCOPTS=-O2`,`LIBS=-L$(NVROOT)/math_libs/lib64 -llapack -lblas` = **nvhpc 参考版 LAPACK/BLAS,无线程化** **[实证]**。
- **disentangle 主循环**(disentangle.F90 `dis_extract`,迭代循环 ~1703 起):每迭代、每 k:由邻居重叠构建 Z 矩阵——成对 zgemm(1780/1783、1965/1968:`('C','N')` num_wann×ndimwin×ndimwin ≈ 253×332×332,乘 nntot=8 个邻居),然后 **ZHPEVX** 对维度 ndiff(≤332,= ndimwin−ndimfroz)的矩阵全谱对角化(1826)。Ce 自旋↓:1000 迭代 × 32 k ≈ 3.2 万次 ~300 维 ZHPEVX + ~50 万次 253×332×332 zgemm 片段 → 1 743 s;自旋↑ 早收敛(Δ<1e-10 连续 3 迭代)→ 374 s **[实证]**。
- **wannierise 主循环**(wannierise.F90 `wann_main`,:416 `do iter=1,num_iter`):每迭代、每 k:`internal_new_u_and_m` 里 **zheev(num_wann=253)**(:1228,梯度步的酉化)+ ~4 次 `utility_zgemm`(253³)(:1246-1301)+ 每邻居 2 次 zgemm 更新 M(:1299-1301,×8 邻居)。Ce:100 迭代 × 32 k → 342/357 s(两自旋)**[实证]**。
- 收尾一次性 num_wann²×num_kpts×nrpts 的 Fourier(:1002/1036)+ 写 .chk(338 MB)。

### Q5 LEOS 段(eos.F)

- `patch/v3/eos.F`:`!$ACC` = 0,无 GEMM/HEEV 调用(生产树 leos_v2 的 eos.F 亦为 0)→ **CPU** **[实证]**。
- `EOS_ANALYZE`(eos.F:195+):循环量级 NUM_WANN×NKPTS×NB_TOT(占据统计)与 NUM_WANN×NIONS(中心→离子配对),253×32×332 ≈ 2.7M 次标量运算,毫秒-秒级;含 EOSCAR 写盘。
- **实测**:w90 结束 07:48:52 → EOSCAR 07:48:56,**~4 s(0.02%)** **[实证]**。挂钩:`leos_v3.patch` 在 main.F MLWF 块后插 `EOS_ANALYZE`,并由 `EOS_PREPARE_MLWF`(eos.F:153)自动置 NUM_WANN=占据带数、默认 `LWANNIER90_RUN=.TRUE.`(eos.F:148)。

### Q6 实测耗时分解

**Ce 单点(ce_fs,17 000 s 墙钟)** — 见 ② 表。要点:SCF 77.2%;Wannier 管线 22.2%(MMN 2.8% + SCDM 2.8% + w90↑ 4.2% + w90↓ 12.4% + LEOS 0.02%)。ce_is/ce_ts 未逐一取数,同配置同量级 **[推断]**。
细节:rank0 的 LOOP+ CPU 时间(33 207 s)≈ 2× 墙钟,反映 w90 单核运行期间其余 rank 在 MPI 忙等烧 CPU **[推断]**。

**m5b 逐帧基准(ws_bench,Li2MnO3,4 rank×4 GPU,40 MD 步,每 5 步分析一帧,共 9 个分析帧)[实证]**:

| 量 | warm-start ON | warm-start OFF |
|----|--------------|----------------|
| 全轨迹墙钟 | 6 742.6 s | 6 532.4 s |
| SCF(LOOP 总和) | 1 794 s(26.6%) | 1 792 s(27.4%) |
| w90 总和(两自旋 9 帧) | 4 313 s(**64.0%**) | 4 104 s(**62.8%**) |
| 分析帧 LOOP+ 均值 | 590 s | ~580 s |
| 其中 w90 均值/帧 | 479 s(**81%**) | 456 s |

- num_wann=num_bands=352 → 除首帧外无需 disentangle(dis 仅首帧:34 s + 140 s),后续帧成本 ≈ 纯 wannierise ×2 自旋 ≈ 420-500 s/帧,与记忆中的"逐帧 554 s"对应(554 ≈ 分析帧 LOOP+ 均值)。
- EOS_WARM_START on/off 对 w90 耗时无实质改善(ON 反而略慢),说明逐帧耗时卡在 wannierise 迭代本身而非初值 **[实证]**。

---

## ④ 与用户直觉的对照

**直觉 1:"SCDM、解纠缠和幺正变换耗时"——大方向对,但内部权重要修正。**
- 单点里这三样合计 3 770 s(22%),分布极不均匀:**wannier90(解纠缠+幺正/局域化)2 817 s 占管线 75%**,SCDM 本身只有 ~470 s,MMN ~482 s。最大单块是自旋↓ disentangle 的 1 743 s——而且它是**打满 DIS_NUM_ITER=1000 上限未收敛**才停的,属于"算法在空转"而非算力不足 **[实证]**。
- 单点总账里 SCF(77%,GPU)仍是绝对大头;但**逐帧轨迹模式下直觉完全成立**:w90 占 63-81%,是首要优化目标 **[实证]**。
- 可加码的事实:这 2 817 s 是**1 个 CPU 核 + 参考版 BLAS/LAPACK** 跑出来的。上限估算 **[推断]**:换线程化 BLAS(nvpl/OpenBLAS/MKL,zheev/zgemm 多线程)即可数倍;dis/wannierise 的每 k 运算天然可按 k 并行(w90 本身有 MPI 版 postw90/wannier90.x,但库模式 `wannier_run` 是串行入口);GPU 化(cuSOLVER zheevd + cublas zgemm,253-332 维批量)在 32 k 批量下也有一个量级空间——但先修"未收敛打满 1000 迭代"和"参考 BLAS 单核"这两个低垂果实,收益/工作量比最高。
- SCDM 段的 470 s 里,矩阵代数(QRCP/SVD,单核)与数据组装(21 248 次 逐带 FFT+全网格广播+规约,NSTRIP=1)约各占一半量级 **[推断]**;后者是通信/结构问题,GPU 化波函数组装意义有限,改成分块/批量(NSTRIP>1、按带聚合规约)更直接。

**直觉 2:"初始波函数还在 CPU"——字面上部分属实,性能上不成立。**
- 属实的部分:随机数填充(wave.F WFINIT,host)、WAVECAR 读盘(fileio.F,host)、初始电荷 SET_CHARGE(GPU 包装被 `_OPENACC_TST` 编译开关剔除)确为 CPU **[实证]**。
- 不成立的部分:初始正交化 ORTHCH 是自包含 GPU 区域(自己上卡/下卡),初始子空间迭代在 ELMIN 内 GPU;整个初始化 ~53 s、占 0.3%,不值得动 **[实证]**。
- 另一个容易被直觉忽略的事实:本次 Ce 生产跑 INCAR 虽写 ISTART=1,**WAVECAR 是 0 字节,实际全部从随机初猜冷启动**(print-out "WAVECAR not read")。554 步 CG 的 13 120 s 里有相当一部分是冷启动代价——若逐帧/复算场景能真正续用 WAVECAR(LWAVE 目前 =.FALSE.),SCF 段可显著缩短(m5b 暖波函数帧的 SCF 份额即是佐证)**[实证+推断]**。
- 真正"还在 CPU"的是**从 wannier_setup 到 EOSCAR 的整段 Wannier 管线**(22%),并且其 75% 收缩在 1 个核上;这才是与直觉里"初始波函数"对应的、事实上的 CPU 短板。

---

### 附:关键证据索引

- GPU 开关机制:openacc.F:41-42(0 层默认 OFF)、choleski2.F:66/372(ORTHCH 自管)、electron_all.F:241/246/933(SCF 区域)
- SCF→host 移交:electron_all.F:630 与 906(UPDATE SELF)、929-930(删除 device 副本)
- mlwf.F/eos.F 零 ACC:本地与 `~/apps/vasp.6.4.3_leos_v2/src` 双向核实
- SCDM:mlwf.F:3842(单 rank)、3887-3919(xGEQP3)、3959-4000(xGESDD+MATMUL)、3186(全网格规约);wave_high.F:4021/4030(逐带 FFT+广播)
- MMN:mlwf.F:1240-1242(KPAR 分工)、1250(PEAD_CALC_OVERLAP)、1272(903 MB allreduce);pead.F:4151-4356
- w90:mlwf.F:1875-1913(master-only、双自旋顺序、结果规约);make.inc(nvfortran -O2、参考 -llapack -lblas);disentangle.F90:1780-1834、wannierise.F90:416/1228-1301
- 实测:ce_fs OUTCAR(LOOP×554=13 120 s、LOOP+ 16 947 s、Elapsed 17 000 s)、wannier90.{1,2}.wout(374/342、1743/357 s)、目录 mtime 时间线(06:46:01→06:54:03→07:01:55→07:13:51→07:48:52→07:48:56→07:49:57);ws_bench on/off OUTCAR+wout
