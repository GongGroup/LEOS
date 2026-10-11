# LEOS VASP 补丁

`leos_v3.patch` 是对原版 VASP 6.4.3 的单个 diff，改动 4 个文件：
`src/.objects`、`src/eos.F`（新增）、`src/main.F`、`src/mlwf.F`。打上补丁后，
VASP 在计算末尾调用 wannier90 库做 Wannier 化，并输出 EOSCAR 文件。当前
版本为 v3.1 加 SCDM 分块组装（见"版本记录"）。

## 安装

前提：wannier90 3.1 的**串行版**静态库（见"运行注意"第一条）。

```bash
tar xzf vasp.6.4.3.tgz && cd vasp.6.4.3/..   # 顶层目录名必须是 vasp.6.4.3
patch -p1 < leos_v3.patch
cd vasp.6.4.3
# makefile.include 中加入：
#   CPP_OPTIONS += -DVASP2WANNIER90
#   WANNIER90_ROOT ?= /path/to/wannier90-3.1.0
#   LLIBS += -L$(WANNIER90_ROOT)/lib -lwannier
#   INCS  += -I$(WANNIER90_ROOT)/include
make DEPS=1 -j8 std gam
```

SAI 集群上的构建环境：nvhpc 24.3 + openmpi 5.0.8 + fftw 3.3.10（module
`nvhpc/24.3-pgi-tuned`、`openmpi/5.0.8-nvhpc24.3-pgi-auto`、`fftw/3.3.10`），
一次完整编译约 5–11 分钟。

## INCAR 参数

最小设置是 `LEOS = .TRUE.`。打开后，下面四个 VASP 参数如果 INCAR 里没写，
会被自动设为 `.TRUE.`：`LWANNIER90`、`LWANNIER90_RUN`、`LSCDM`、
`LWANNIER90_AUTO_WINDOW`。`NUM_WANN` 不写时按每个自旋的占据带数自动设定。

| 参数 | 默认 | 作用 |
|---|---|---|
| LEOS | .FALSE. | 总开关 |
| EOS_FLOW | .TRUE. | 在 VASP 内做最小费用流分配，结果写入 EOSCAR 的 OS(flow) 列（β = 10，Allen 电负性，与 pyleos 默认设置相同） |
| EOS_WF_DETAIL | .TRUE. | 写出 #WF 块：每个 WF 的中心、spread、占据数。pyleos 依赖它 |
| EOS_PROJ_WEIGHT | .TRUE. | 写出 #WP 块：每个 WF 在各原子 PAW 投影上的权重，供投影辅助判定使用 |
| EOS_EVERY_N | 0 | 0：只在最后一步分析（含单点计算）；1：每步分析；N > 1：每 N 步和最后一步分析 |
| EOS_NELECT_WARN | 0.01 | Wannier 电子数与 NELECT 之差超过此值时警告 |
| EOS_NELECT_ABORT | 0.35 | 差值超过此值时写出诊断版 EOSCAR 并中止 |
| EOS_SCDM_ONLY | .FALSE. | 跳过 wannier90 的迭代局域化，只用 SCDM 初猜的中心（快速筛查用） |
| EOS_WARM_START | .TRUE. | EOS_EVERY_N > 0 时，把上一个分析步的 WF 中心作为引导中心传给 wannier90 |
| EOS_ENEG_TOL | 0.05 | 贪心分配：电负性差小于此值时视为相同 |
| EOS_ONSITE_CUT | 0.35 | 贪心分配：WF 中心离原子小于此距离（Å）时判为该原子的孤对 |
| EOS_O_MAX | 8.0 | 贪心分配：O 最多接收的电子数 |
| EOS_ENEG_PAULING | .FALSE. | 贪心分配改用 Pauling 电负性 |

后四个参数只影响 EOSCAR 中的贪心列 OS(int)，不影响 OS(flow)。

一个可用的单点 INCAR（data/v3/m4_wave2/as4s4/INCAR）：

```
PREC = Accurate
ENCUT = 500
EDIFF = 1e-06
ISMEAR = 0
SIGMA = 0.05
ISYM = 0
NSW = 0
IBRION = -1
LWAVE = .FALSE.
LCHARG = .FALSE.
LEOS = .TRUE.
LSCDM = .TRUE.
LWANNIER90_AUTO_WINDOW = .TRUE.
EOS_WF_DETAIL = .TRUE.
EOS_PROJ_WEIGHT = .TRUE.
```

## 运行注意

### wannier90 库与离线 wannier90

- **库模式必须链接串行版 libwannier.a。** 补丁的调用方式是只有主进程进入
  wannier90 库，其他进程在库外等待。链接 MPI 版库时，库内部会再发起集体
  通信，与库外等待的进程互相等待，在解纠缠第一步死锁：所有进程 100% CPU
  空转，没有 .werr 文件。很多站点预编译的库是 MPI 版（如 Chinahpc
  /data/app 下全部版本）。检查方法：`nm libwannier.a | grep -ci mpi_`，
  结果非零就是 MPI 版，需要用 `COMMS=`（留空）重新编译串行版。
- **离线运行的 wannier90.x 用 MPI + 线程化 MKL 构建。** 上一条限制只针对库
  模式。make.inc 要点：`F90=ifort`、`COMMS=mpi`、`MPIF90=mpiifort`，LIBS
  用 `-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5`。运行时设
  `MKL_NUM_THREADS=<n> MKL_DYNAMIC=FALSE`。两个自旋并发运行时每进程 2 线程
  最快（数据见 docs/RESULTS.md"效率"）。gnu + OpenMPI 5 编译的 wannier90 在多
  进程下会崩溃，要用 ifort + Intel MPI + MKL。SAI 集群上的现成构建在 LEOS
  目录下的 `w90_offline/build/wannier90-3.1.0-impi-thread/`；离线流程的说明
  在 `w90_offline/README.md`。

### 常见失败

- **NELECT mismatch too large。** 原因几乎都是解纠缠缺少冻结窗口。设置
  `LWANNIER90_AUTO_WINDOW = .TRUE.`（配合 LSCDM），不要手写 wannier90.win。
  中止时补丁会写出诊断版 EOSCAR（带 `# LEOS_ABORTED` 标记），其中 `#NLD`
  块列出丢失最多的（自旋，k 点，能带）。如果丢失集中在占据带顶部、分散在
  所有 k 点上，就是缺冻结窗口。
- **`kmesh_get_bvector: Not enough bvectors found`。** 体心或倾斜的小晶胞，
  以及某些分子晶体的大晶胞，在特定 k 网格下会通不过 wannier90 的 b 向量检查。
  依次尝试：换成偶数或与晶胞形状匹配的网格；或在 INCAR 中加
  `WANNIER90_WIN = "\nkmesh_tol = 0.00001\n"`（把默认的 1e-6 放宽到 1e-5；
  这个多行字符串也可以同时写入 dis_num_iter 等其他参数）。不要改用 Gamma-only
  来回避：Gamma 模式下局域化行为不同，会产生大量 stray。
- **#WP 数据。** 2026-07-11 修复（提交 e0ad691）之前，所有多进程运行产出的
  #WP 块都是错的：带并行时先取模方再跨进程求和，丢掉了交叉项。这些数据要用
  修复后的版本重新计算。另外，**GPU（nvhpc）版和 CPU（Intel）版产出的 #WP
  有系统性差异**，每原子聚合值相差约一倍，原因未查明，所以不同平台的 #WP
  不能混用；需要 #WP 的计算统一在一个平台上完成。

### 内存与 SCF

- **KPAR 与内存（GPU 版）。** KPAR 不会把波函数按 k 点分开存放，每个进程
  的内存占用反而与 KPAR / 进程数成正比。进程少时 KPAR > 1 会耗尽主机内存：
  Ce 体系在 4 GPU、KPAR = 2 下启动即内存不足。4 GPU / 62.5 GB 时用
  KPAR = 1；8 GPU / 125 GB 时 KPAR = 2 可用。
- **磁性体系写全 MAGMOM。** 可以避免 NBANDS 被默认规则调高（Ce 体系从 332
  降到 308），也让 SCF 更快收敛。
- **Ce18O36Co4CH4 末态的实测配置。** KPAR = 1、NBANDS = 308、写全 MAGMOM，
  4 GPU 上用时 5 小时 24 分钟，内存峰值 49.2 GB（配额的 79%）。两个自旋的
  解纠缠各 100 步收敛（自旋共享冻结窗口时，第二自旋跑满 1000 步、耗时
  1743 秒仍未收敛）；按自旋设置的 num_wann 为 255 / 246；dNe 由 0.193 降到
  0.0008，ghost 由 6 降到 1。
- **Ce 体系的 SCF 设置（四组对比，Ce 初态）。** ALGO = ALL 是必需的，换成
  Davidson 会在 176 步内发散到 −10⁷ eV 量级。ALGO = ALL 下混合参数影响很小：
  保守设置 545 步、默认混合 526 步、温和 Kerker 511 步，相差 3–6%。有效的
  加速只有写全 MAGMOM 初值（末态从 539 步降到 454 步，少 16%）。像与像之间
  复用波函数尚未测试（issue #18）。

### SLURM 作业

- **相对路径。** `--output=print-out` 这类相对路径以提交时所在目录为准，作业
  的工作目录也是。从家目录提交别人的脚本时，VASP 会在家目录找 INCAR 并立即
  退出，错误信息在 ~/print-out 里。提交前先 cd 到算例目录，或在脚本中使用
  绝对路径并显式 cd。
- **资源参数写在脚本里。** SAI 集群的检查钩子匹配的是提交命令行：在命令行上
  给 `--cpus-per-task` 的作业会在 2 秒内被 root 取消，连输出文件都不生成；
  同样的值写在脚本的 `#SBATCH` 行里则正常。`--export=ALL,VAR=x` 不受影响。
- **重新编译二进制。** 如果有排队中的作业引用同一个二进制路径，在删除旧
  目录到编译完成之间，它们可能被调度运行并立即失败（已在运行的作业不受影响）。
  重新编译前先用 squeue 检查，或部署到新路径。

### 开发注意

- VASP 内 MPI 求和函数 M_sum_z8 / M_sum_d8 的计数参数 n 是 INTEGER(qi8)，
  调用时必须写 `SIZE(x, KIND=qi8)`。写成 `SIZE(x)`（32 位）时高位是栈上的
  垃圾值，会导致堆越界和 SIGKILL，而且串行运行不会出现，很难排查。

## 集群部署

SAI 集群上的补丁版 VASP 在
`/org/sjtu-gongxueqing/zhangxianyang/temp/org/apps/vasp.6.4.3_leos_v31_20260711/bin/`
（vasp_std、vasp_gam），由提交 e0ad691 时的 leos_v3.patch 编译。之后加入的
SCDM 分块组装（提交 54eec66）只改变速度、不改变结果，还没有重新部署。

软链接 `/org/sjtu-gongxueqing/zhangxianyang/temp/org/apps/leos_current` 和
LEOS 目录下的旧作业脚本仍指向家目录迁移前的路径
`/home/sjtu-gongxueqing/zhangxianyang/`，使用前要改成
`/org/sjtu-gongxueqing/zhangxianyang/temp/org/`。

## 版本记录

### v3.1 之后

| 改动 | 内容 |
|---|---|
| SCDM 分块组装（提交 54eec66，issue #17） | SCDM 收集轨道时按约 16 个能带一组组装，用一次集体通信代替逐能带的广播和归约。每个缓冲元素只有一个非零贡献，结果逐位不变 |

### v3 → v3.1

| 改动 | 内容 |
|---|---|
| 中止前写诊断版 EOSCAR | 电子丢失超过 EOS_NELECT_ABORT 时，先写出带 `# LEOS_ABORTED` 标记和完整丢失记录的 EOSCAR 再跳过。pyleos 读到时给出说明；轨迹中个别帧中止不影响其他帧 |
| 弥散型丢失定位 `#NLD` | 缺冻结窗口时，占据带会在所有 k 点上连续泄漏，这是最常见的丢失方式。`#NLD` 块记录丢失位置；MUMNOX 实测 1348 条记录，4.35 个电子定位到"自旋 1 第 117 带、自旋 2 第 124 带，分散在全部 k 点"。窗口型丢失的 `#NL` 块也写入 EOSCAR |
| 按自旋设置 num_wann | 磁性体系的少数自旋不再产生 ghost（VIRCEF 第二自旋的 ghost 归零，Ce 恢复为 +3）。同时按自旋设置自动冻结窗口；若改写 .win 副本失败，退回共享值；INCAR 显式给出的 NUM_WANN 保持原有含义 |
| SCDM 不再双缓冲 | 0 号进程不再同时持有两份全网格矩阵（Ce 体系该段内存增量约 8 GB），结果逐位相同 |
| M/A 矩阵及时释放 | 非主进程归约后即释放 M（每进程 −0.9 GB）；wannier90 结束后释放 M/A（每进程常驻 −1 GB）；A 只在 SCDM 路径下释放，不影响 LOCPROJ、elphon |
| 原位最小费用流：OS(flow) 列 | EOSCAR 逐原子表最后一列为 VASP 内完成的最小费用流分配（β = 10，Allen 电负性，镧系用与 pyleos 相同的电负性表），开关 EOS_FLOW 默认开。与 pyleos 逐原子相同：7 个测试的 Python 对照代价逐位相同，CO2、VIRCEF 实机验证 |
| #WP 带组归约修正（提交 e0ad691） | 多进程带并行时，先把各带组的复数分量跨组求和，再取模方。串行结果逐位不变；修正后串行与多进程的每原子聚合值相差 2×10⁻⁴（收敛噪声量级），两个集群回归通过 |
| 批跑与离线流程 | 见 pyleos/README.md 的 batch.py 一节和 script/offline_w90_prep.py |

### v2 → v3

| 改动 | 内容 |
|---|---|
| 删除匈牙利分配分支 | 设计缺陷（阴离子拿不到多于中性时的电子），从未被使用 |
| 按（自旋，k 点，能带）诊断电子丢失 | 不只报告总数，可以定位 |
| EOSCAR 的 #WF 块 | 每个 WF 的中心、spread、占据数，是 pyleos 的输入 |
| #WP 投影权重块 | 投影辅助判定的数据来源 |
| EOS_SCDM_ONLY | 快速筛查模式；CO2 上 wannier90 耗时 2 秒 → 0.023 秒，氧化态不变 |
| 轨迹模式接线 | 原来主程序写死 `IF (NSTEP == 1)`，EOS_EVERY_N 从未生效；现在按 EOS_DO_THIS_STEP 逐步分析（只开 LWANNIER90 的用户仍只在第一步运行），并修复了每步无条件 MLWF_FREE 的悬空引用 |
| 热启动 | 上一步的 WF 中心作为引导中心 |
| 其他 | 周期回卷公式、自动窗口容差、编译警告 |

## 验证记录

### v3（2026-07-07，SAI 集群）

1. 来源：原版 6.4.3 与生产树 `vasp.6.4.3_leos_v2`（剥离 CCQN 后）的 diff。
   生产树是 M4/M5 全部数据的出处，因此以它为准。
2. 闭环：原版 + leos_v3.patch 与剥离后的生产树逐字节一致；生产树与剥离版
   的差异恰好是 CCQN 的文件集，没有混入 LEOS 内容。
3. 编译：在打过补丁的原版树上全新编译通过（作业 594329，约 11 分钟）。
4. 回归：CO2 和 TiO2 的同一输入分别用生产二进制和新编译二进制计算，EOSCAR
   整数氧化态完全相同，数值字段在 1e-6 内一致（差异是不同编译产物之间的浮点
   末位差别），#WP、#WF 块齐全。

合并时的发现：

- 生产树上曾叠加过另一个项目（CCQN-TS）的补丁，leos_v3.patch 已剥离。
- 生产树中有一处未入库的修改 `IF (LEOS) LSCDM = .TRUE.`，已收入
  leos_v3.patch。
- 仓库中旧版 leos.patch（v2）的 main.F 轨迹部分是开发中的旧版本，缺少只开
  LWANNIER90 时的 NSTEP == 1 保护，与生产实际不符。

### v3.1（2026-07-10，SAI 集群）

1. 编译：原版树 + 本补丁三次编译通过（作业 612186、612454、612747，各约
   5–11 分钟）。改动经过两轮代码评审，记录在 git 历史中。
2. 回归：CO2 整数氧化态完全相同，数值字段差异在 1e-14 量级（判据仍为"整数
   相同、数值 1e-6 内"；GPU 版在不同节点之间有浮点末位差异，做不到逐字节
   相同）。VIRCEF 第二自旋 ghost 由 2 降到 0，Ce 恢复 +3。MUMNOX 在无窗口和
   错误窗口两种设置下分别验证了 `#NLD` 和中止行为。
3. 内存：Ce18O36Co4（63 原子，ISPIN = 2）在 4 GPU、KPAR = 1、NBANDS = 308
   下 VASP 总内存 41 GB（62.5 GB 配额的 66%）。
4. 回归数据：data/v3/regress_v31/。

## 归档文件

- `leos.patch`：v2 补丁，仅存档。其中 main.F 部分是旧版本（见上），不要使用。
- `v3/`：v2 → v3 的增量材料（eos.F 全量替换和 main.F / mlwf.F 注入脚本），
  已并入 leos_v3.patch，见 [v3/README.md](v3/README.md)。
