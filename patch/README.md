# LEOS VASP 补丁

## 正式发布物：leos_v3.patch

单一 diff，从原版 VASP 6.4.3 一步得到 LEOS v3。触及 4 个文件：
`src/.objects`、`src/eos.F`（新增）、`src/main.F`、`src/mlwf.F`。

### 从零安装

```bash
tar xzf vasp.6.4.3.tgz && cd vasp.6.4.3/..   # 顶层目录名须为 vasp.6.4.3
patch -p1 < leos_v3.patch
cd vasp.6.4.3
# makefile.include 需要 wannier90 v3.1 静态库（LEOS 依赖 SCDM 与库模式接口）：
#   CPP_OPTIONS += -DVASP2WANNIER90
#   WANNIER90_ROOT ?= /path/to/wannier90-3.1.0
#   LLIBS += -L$(WANNIER90_ROOT)/lib -lwannier
#   INCS  += -I$(WANNIER90_ROOT)/include
make DEPS=1 -j8 std
```

### 验证记录（2026-07-07，sai-pod）

1. **来源**：diff(原版 6.4.3, 生产树 `vasp.6.4.3_leos_v2` 剥离 CCQN 后)。
   生产树是 M4/M5 全部数据的出处，故以它为准而非仓库补丁材料的叠加。
2. **闭环**：原版 + leos_v3.patch 与剥离后的生产树逐字节一致；
   生产树与剥离版的差异恰为 CCQN 文件集（无 LEOS 内容混入）。
3. **编译**：从打过补丁的原版树全新编译通过（作业 594329，约 11 分钟，
   nvhpc 24.3 + openmpi 5.0.8 + fftw 3.3.10）。
4. **双树回归**：CO2 与 TiO2 同一输入分别经生产二进制与新编译二进制，
   EOSCAR 整数氧化态精确相等、数值字段 1e-6 内一致（残差为编译产物间
   的浮点末位抖动），#WP/#WF 数据块齐全。

### v3.1 验证记录（2026-07-10，sai-pod）

1. **改动**：六项——abort 先写诊断版 EOSCAR（`# LEOS_ABORTED` +
   窗口型 `#NL` + 弥散型 `#NLD` 丢失定位）；per-spin num_wann 幽灵
   修复（含 per-spin 自动冻结窗口与 win 改写回退的硬不变量）；
   SCDM 免双缓冲；M/A 矩阵生命周期收紧。敌对评审两轮
   （敌对评审两轮通过，评审记录见 git 历史 2026-07-10）。
2. **编译**：干净原版树 + 本补丁三轮全过（作业 612186/612454/612747，
   各约 5-11 分钟）。
3. **回归**：CO2 整数 OS 精确相等、数值字段 1e-14 级（判据沿用 v3
   验证口径"整数相等 + 1e-6 内"；GPU 版存在跨节点运行间浮点末位
   抖动，逐字节不可达）；VIRCEF spin2 幽灵 +2→0、Ce 恢复 +3；
   mumnox 无窗口/错窗口两种复现各验证 `#NLD`/abort 行为。
4. **内存**：Ce18O36Co4（63 原子 ISPIN=2）在 4 GPU + KPAR=1 +
   NBANDS=308 下 VASP 总 RSS 41 GB（62.5 GB 配额的 66%）。

### 使用注意（实测积累）

- **wannier90 库必须是串行编译版**（2026-07-11 Chinahpc 实测）：补丁的
  库模式调用约定是"仅主进程进库、其余进程在库外等归约"。链接 MPI 版
  libwannier.a（很多站点的预编译库都是,如 Chinahpc /data/app 全系)会在
  解纠缠第一步全员死锁——库内部再次发起集体通信,与库外等待的进程互等,
  四进程 100% CPU 空转无 .werr。换集群先验库:
  `nm libwannier.a | grep -ci mpi_`,非零即 MPI 版,须用站点源码
  `COMMS=`(空)自编串行版重链。
- **离线 wannier90.x 用"MPI + 线程化 MKL"构建**(与上条不冲突:
  串行限制只针对库模式;离线跑的独立 wannier90.x 可以吃满 CPU)。
  make.inc 要点:`F90=ifort  COMMS=mpi  MPIF90=mpiifort`,LIBS 用
  `-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5`——与
  串行 MKL 版只差 threading layer 一行;运行时
  `MKL_NUM_THREADS=<n> MKL_DYNAMIC=FALSE`。实测(288 函数帧,
  data/v3/w90_threads_scan_20260712.tsv):参考版数学库对线程零
  响应;MKL 版 8 进程 × 4 线程时解纠缠/局域化每迭代对单核参考版
  约 16×/13×;两自旋并发的生产形态每 rank 2 线程即最优(再加线程
  内存带宽打满反而变慢)。判定一致性:收敛通道与参考版逐位一致。
  sai-pod 现成构建:`~/LEOS/w90_offline/build/wannier90-3.1.0-impi-thread/`。
- **#WP 投影数据的两条铁律**(2026-07-11 定谳,提交 e0ad691):
  ① 修复前(v3.1 含)所有多进程运行的 #WP 块全部无效——带并行下先取
  模方后归约丢交叉项,输出是平台相关伪数据;修复后多进程=单进程,
  历史 #WP 须用修复版重产。② GPU(nvhpc)版与 CPU(Intel)版的 #WP 存在
  独立的系统性平台差异(每原子聚合约差一倍,根因待查)——**不同平台产的
  #WP 永不混用**,仲裁/档 2 重产统一在单一平台完成。
- **VASP 内 MPI 归约的计数参数**:M_sum_z8/M_sum_d8 系列的 n 是
  INTEGER(qi8)——调用必须写 `SIZE(x, KIND=qi8)`;裸 `SIZE(x)`(32 位)
  高位吃栈垃圾,堆越界 SIGKILL,且串行不触发、极难排查。
- **bvector 失败**（`kmesh_get_bvector: Not enough bvectors found`）：
  体心/倾斜小胞或某些分子晶体大胞的特定 k 网格会触发 wannier90 的
  b 向量完备性检查失败。两个对症解，按序尝试：换偶数/匹配胞形状的
  网格；或 INCAR 注入 `WANNIER90_WIN = "\nkmesh_tol = 0.00001\n"`
  （默认 1e-6 放宽到 1e-5，多行字符串可同时注入 dis_num_iter 等）。
  不要降到 Gamma-only 回避——Gamma 模式下局域化行为不同，实测会
  产生大面积 stray。
- **NELECT mismatch too large**：根因几乎总是解纠缠缺冻结窗口。开
  `LWANNIER90_AUTO_WINDOW = .TRUE.`（配 LSCDM）即根治；不要手写
  wannier90.win。v3.1 起 abort 也会写出诊断版 EOSCAR：`#NLD` 块
  直接列出丢失最大的 (自旋, k, 能带)——丢失集中在占据流形顶、
  弥散于全部 k 点，就是缺冻结窗口的指纹。
- **KPAR 与内存（GPU 版铁律）**：KPAR 不切分 k 方向的波函数存储，
  每 rank 注册量反而正比于 KPAR/Nranks——rank 数少时 KPAR>1 会
  直接把主机内存顶爆（Ce 系 4 GPU + KPAR=2 启动即 OOM 的根因）。
  4 GPU/62.5 GB 档必须 KPAR=1；8 GPU/125 GB 档 KPAR=2 可用。
  磁性体系补全 MAGMOM 可避免 NBANDS 被默认规则抬高（332→308）。
- **sbatch 相对路径陷阱**：`--output=print-out` 这类相对路径按
  **提交时所在目录**解析，作业工作目录同理——从家目录重提别人的
  脚本，VASP 会在家目录里找 INCAR 秒死，真实报错还藏在 ~/print-out。
  重提一律先 cd 进案例目录，或脚本内用绝对路径 + 显式 cd。
- **资源参数写进脚本，别放 sbatch 命令行**（sai-pod 实测）：站点
  查杀钩子匹配的是提交命令行——命令行给 `--cpus-per-task` 的作业
  连 CPU 分区也会在 2 秒内被 root 取消（输出文件都不生成）；同样
  的值写进脚本 `#SBATCH` 行则放行（32/64 核对照实验双向验证）。
  `--export=ALL,VAR=x` 传环境变量不受影响；资源类开关一律进脚本。
- **重建二进制的空窗**：排队中（PENDING）的作业引用同一二进制
  路径时，`rm` 旧树到编译完成之间的窗口会让它们被调度进来秒死
  （运行中的作业不受影响，二进制已载入内存）。重建前查 squeue
  有无待跑的消费者，或部署到全新路径。
- **Ce 系 SCF 配方定论（2026-07-10 四组对照）**：ALGO=ALL 是必需品
  （换 Davidson 直接发散到 −10⁷ eV 量级）；ALL 模式下混合参数
  无关紧要（极端保守 545 步 vs 默认混合 526 vs 温和 Kerker 511，
  仅差 3-6%）；实证有效的加速是补 MAGMOM 初猜（FS 像 539→454 步，
  −16%）。像间 WAVECAR 接力未测，留为候选。

### 合并时的审计发现

- 生产树上曾叠加过另一个项目（CCQN-TS）的补丁，`leos_v3.patch` 已剥离；
- 生产树含一处未入库的 LEOS 修改 `IF (LEOS) LSCDM = .TRUE.`（LEOS 开启时
  自动启用 SCDM 初猜），已收入 leos_v3.patch；
- 仓库旧版 `leos.patch`（v2）里 main.F 的轨迹集成块是开发中的旧变体
  （缺纯 LWANNIER90 用户的 NSTEP==1 保护），与生产实际不符——v2 补丁
  仅作历史归档，勿再使用。

## 归档材料

- `leos.patch` — v2 历史补丁（见上：main.F 部分为旧变体，仅归档）。
- `v3/` — v2→v3 的增量材料（eos.F 全量替换 + main/mlwf 注入脚本），
  已并入 leos_v3.patch，保留备查。
