> 归档注记：2026-07-10 由会话内 sub-agent（fortran-review2）生成；原路径 scratchpad/fortran_review.md；评审对象为五项改动的 BLOCKER 修复前状态，两个 BLOCKER 与四个 MAJOR 的修复已在同日落码（见 20260710_index.md）。

# LEOS v3.1 五项改动敌对式评审报告（2026-07-10 重跑）

评审对象：`scratchpad/patched_tree/vasp.6.4.3/src/{mlwf.F, eos.F}`（v3 基线 + 五项新改动），
基线快照 `mlwf.F.v3base` / `eos.F.v3base`。增量共 mlwf.F +117/eos.F +76 行。
交叉证据：`vasp_src/vasp.6.4.3.constrain/src/`（主干源码）、wannier90 3.1.0 官方源
（v3.1.0 tag 的 `wannier_lib.F90`/`parameters.F90`/`disentangle.F90`/`io.F90`，已下载到
本次评审 scratchpad 逐行核对）、`pyleos/parsers.py`、`patch/leos_v3.patch`（工作区 vs HEAD）。

**总判定：2 BLOCKER / 4 MAJOR / 9 MINOR。改动 2（幽灵修复）在其自身的目标场景
（自旋极化且两通道占据数不等，即 VIRCEF / Ce 三像）下必然触发 wannier90 库内
硬停或内存越界，不能按现状合入；改动 1/4/5 与改动 3 的 Fortran 侧通过，
改动 3 的 pyleos 侧有一个轨迹场景 MAJOR；补丁完整性全部通过。**

---

## 0 补丁完整性（评审点 6）— 通过

用"实际 apply 后逐字节比对"的方式验证，不依赖 diff 文本目测：

- 工作区 `patch/leos_v3.patch` == 上一进程产物 `scratchpad/leos_v3_new.patch`（逐字节）。
- `patch/v3/eos.F` == `patched_tree/.../eos.F`（逐字节，"两处同步"约定满足）。
- 把工作区补丁与 `git show HEAD:patch/leos_v3.patch` 分别 apply 到干净 `orig_tree`：
  - **`.objects` 与 `main.F` 两个版本 apply 结果逐字节相同**——工作区补丁对这两个文件
    只有 hunk 行号重排（-908→-911 等），语义零变化；
  - 工作区补丁 apply 结果与 `patched_tree` 四个文件全部逐字节相同——**被评审的树
    就是补丁的产物**；
  - HEAD 补丁 apply 出的 mlwf.F/eos.F 与 `.v3base` 快照逐字节相同——基线真实。
- 行宽：新增行最长 mlwf.F 114 列 / eos.F 124 列，均低于两文件既有最长行
  （179 / 136 列，v3base 已在生产工具链编译通过），无编译行宽风险。
- 续行风格：eos.F 新增代码用"行尾 `&` + 续行行首 `     &`"，mlwf.F 用"行尾 `&` +
  行首 `&`/`         &`"，各与所在文件上下文一致；自由格式两种皆合法。

---

## 1 改动 1：A1 SCDM 免双缓冲（mlwf.F）

### 1.1 MPI 集体通信次序（评审点 ⑴）— 通过

- 每自旋调用次数两侧严格相等：非 work 段 anchor×2 + 其余 k×1（mlwf.F:4031-4038）；
  work 段 anchor 第 1 次 4072、第 2 次 4117、k 循环 4145，共 2+(NKPTS-1)。
- `GET_WAVE_FUNCTIONS` 内部集体调用次数只依赖 `WDESK%NBANDS`（各 rank 相同）与
  NSTRIP=1：`W1_GATHER_W1`（3299，COMM_INTER 上的 gather，无浮点归约）与
  `M_sum_g(COMM_INB)`（3335），未传 `EXCLUDE_BAND` → 循环计数全 rank 一致。
- **QRCP 与 pivot 重排确认纯本地**：`RANK_REVEAL_QR%CALC` 是裸 ZGEQP3/DGEQP3
  （mathtools.F:325-332，全文件无 CALLMPI）；spinor pivot 重排（mlwf.F:4088-4109）
  纯数组操作；`CUTOFF_DES_GET_CUTOFF`（4174-4195）纯函数。work rank 在 QRCP 期间
  晚到第 2 次 anchor 的集合点，非 work rank 只是等待，**不构成死锁**。
- 两处 `M_bcast`（4015、4041/4170）两分支都到达。KPAR>1 时其它 k 组的 rank 走
  非 work 分支、集体通信在各自 COMM_INB/COMM_INTER 内自洽（波函数已由
  KPAR_SYNC_ALL 全量同步，与基线同构）。
- `GET_WAVE_FUNCTIONS_PHASE`（3374-3425）纯本地（只算指数，无通信），仅 work rank
  调用，与基线相同。

### 1.2 数值等价（评审点 ⑵）— 通过（附一处已知 ULP 级差异，仅非共线）

逐操作序对照（基线行号取自 v3base，经 delta diff 确认）：

- 锚点：旧 = 先 `NORM=DOT_PRODUCT(PSI(:,IB),PSI(:,IB))`、`PSI=PSI/SQRT(NORM)`，再
  `PSI_D_SCDM(IB,:)=GCONJG(PSI(:,IB))*CUTOFF`；新 = `WVFN_BANDS_T` 直写 + 行内
  `/SQRT(BNORM)` 再 `GCONJG(...)*CUTOFF`（4076-4077）。每元素的浮点操作序均为
  除法→共轭（精确）→乘法，操作数相同 → **NRSPINORS==1 时逐位一致**。
- `BAND_NORMS`（3345-3349）分 spinor 两段 DOT_PRODUCT 相加；旧锚点是整向量一次
  DOT_PRODUCT（v3base 3945-3947），旧 k 循环本来就是分段（v3base 4020-4022）。
  → **NRSPINORS==2（非共线）时锚点范数有求和次序差 = ULP 级**；新代码把锚点统一成
  k 循环的分段公式（基线自身本来就不一致）。共线 ISPIN=1/2 完全不受影响。
  LEOS 目标体系全部共线 → 回归"逐字节不变"预期成立。
- pivot 重访的一致性依赖"同一 (k,ISP) 两次重建结果逐位相同"：重建链为
  gather（无归约）→ 网格重排（3308-3334）→ `M_sum_g`（NCORE=1 时单成员，恒等）
  → 同一 FFT 计划。单次运行内确定性成立；pivot 本身是整数不受扰动。
- 第 2 次 anchor 与 k 循环的抽取式 `GCONJG(PIV_VALS/SQRT(BNORM))*[PHASE]*CUTOFF`
  （4122、4155）与旧式 `GCONJG(PSI(PIVOT,IB))*[PHASE]*CUTOFF`（PSI 已归一）
  操作序相同 → 共线逐位一致。

### 1.3 索引空间匹配（评审点 ⑶，重点怀疑项）— 澄清：目标平台成立，基线同风险

- PIVOT 值域 = QRCP 矩阵列空间 1..NP*NRSPINORS（INIT 于 4051，NP=GRID%RL%NP，4045）；
  WVFN 长 NPW*NRSPINORS（NPW=NGX*NGY*NGZ，3290-3291）。
- **源码依据**：波函数网格的通信域是 COMM_INB（main.F:2352 `GRID%COMM =>COMM_INB`）；
  REAL_STDLAY 把 NGX*NGY 个 (x,y) 柱按 `MOD(IND,NCPU)+1` 轮转分配
  （mgrid.F:833-852），`RL%NP = RL%NCOL*RL%NROW`（mgrid.F:857，NROW=NGZ）。
  **NCPU=1（GPU 版强制 NCORE=1）时全部柱都在本 rank → NP = NGX*NGY*NGZ = NPW**，
  两个索引空间恒等，`WVFN_PIVOT(NWRITTEN,I)=WVFN(PIVOT_IDX(I))`（3350-3354）合法。
- **确认任务书里的猜想**：NCORE>1 时基线本来就不自洽——
  v3base 的 `WVFN_BANDS(:,NWRITTEN) = WVFN` 左边 NP*NRSPINORS 行、右边
  NPW*NRSPINORS 元素，形状不符（非法赋值）；`GET_WAVE_FUNCTIONS_PHASE` 用全局
  IND 写入长度 NP*NRSPINORS 的数组（3398/3413）同样越界。stock SCDM 隐含
  NCORE=1 前提，新代码与基线**同风险、不加剧**。
- MINOR-6：建议在 CALC_WANNIER_TRAFO_SCDM 进入时加一行断言
  `IF (NP /= NGX*NGY*NGZ) CALL vtutor%error(...)`，把隐含前提显式化。

### 1.4 分配尺寸与生命周期（评审点 ⑷）— 通过

- `PIV_VALS(NB_TOT,NUM_WANN)`（4049）vs 调用处 `SIZE(PIVOT_IDX)=NUM_WANN`
  （4118/4146）、行下标 NWRITTEN≤NB_TOT ✓；`BNORM(NB_TOT)`（4048）而
  NWRITTEN 恰好走满 NB_TOT（= NBANDS×NB_PAR，无 EXCLUDE_BAND）✓。
- `CUTOFF_FACTOR` 为 ALLOCATABLE（3989），函数结果赋值走 F2003 自动分配
  （基线同款写法）✓。
- `PSI_D_SCDM`：ISP 循环头 `IF (ALLOCATED) DEALLOCATE`（4070）处理第 2 自旋进入时
  仍是 (NB_TOT,NUM_WANN) 旧形状的情形；4115-4116 QRCP 后重分配；子程序返回时
  allocatable 自动释放 ✓。峰值内存从 PSI+PSI_D_SCDM 两份全网格降为一份，
  与 A1 目标（rank0 −5.4 GB）一致。

### 1.5 本节其余发现

- **MINOR-1（API 隐患）**：`WVFN_PIVOT` 存在而 `PIVOT_IDX` 缺省时
  `SIZE(PIVOT_IDX)`（3351）引用不存在的可选虚参 = 非法。现两个调用点均成对传入，
  属防御性问题；建议加 `IF (.NOT.PRESENT(PIVOT_IDX)) vtutor%bug`。
- **MINOR-2（性能）**：`WVFN_BANDS_T(NWRITTEN,:) = WVFN`（3344）是按行写
  (NB_TOT,NP) 大数组 = 步长 NB_TOT 的 5.4 GB 跨步写，逐带一次；与旧"先连续存
  再显式转置"总流量同量级，可接受，但若后续调优可改分块。

---

## 2 改动 2：幽灵修复 per-spin num_wann（mlwf.F + eos.F）

### 2.1 BLOCKER-1：共享 AUTO_WINDOW 冻结窗 × per-spin num_wann → wannier90 硬停

**触发条件：ISPIN=2 且 NREQ_SPIN(1)≠NREQ_SPIN(2)（即幽灵修复自身的全部目标体系），
LEOS 默认配置（LWANNIER90_AUTO_WINDOW 强制 TRUE，eos.F:150）。单点计算即触发。**

证据链：

1. `ESTIMATE_AUTO_WINDOW` 的 ELIM 对**两个自旋合并**取 min/max（mlwf.F:923-930），
   `DIS_FROZ_MAX = 0.5*(ELIM(2,NUM_WANN)+ELIM(1,NUM_WANN+1))`（gap>0.05 分支，
   mlwf.F:939-943）> ELIM(2,NUM_WANN) ≥ 少数自旋第 NUM_WANN 条带在所有 k 的能量
   → **少数自旋在每个 k 点都有 ≥ NUM_WANN 条带落在冻结窗内**（dis_froz_min 缺省
   = 谱底，w90 parameters.F90:1685）。窗值只有一份、写进共享的 wannier90.win
   （mlwf.F:1868-1869），spin 拷贝原样继承。
2. 幽灵修复把少数自旋的 num_wann 改写为 NW2 = NREQ_SPIN(2) < NUM_WANN
   （mlwf.F:2018-2019 → MLWF_WIN_SET_NUM_WANN 改 .win；wannier_run 实参 2027）。
3. wannier90 3.1.0 disentangle.F90:784-791：
   `if (ndimfroz(nkp) .gt. num_wann) ... io_error('dis_windows: More states in the
   frozen window than target WFs')`——NUM_WANN > NW2 恒成立 → 第一个 k 就报错。
4. `io_error` = 写 .werr 后 `STOP`（串行库；io.F90:307-362；MPI 库则 MPI_abort）。
   master rank 在 wannier_run 内死亡，其余 rank 阻塞在 mlwf.F:2038 的
   `M_sum_z8(MLWF%COMM,U_matrix,...)` → 作业挂死到 MPI runtime 检出/超时，
   EOSCAR 不产出。
5. gap≤0.05 的回退分支（mlwf.F:945，froz_max=ELIM(1,NUM_WANN)-1e-5）同样几乎必炸
   （少数自旋位于该能量下的带数一般仍 > NW2，只是逐 k 依赖）。

**波及面**：VIRCEF（spin2 −2）、Ce18O36Co4 三像（ghost 6-10）——恰是验收案例；
非磁与等占据体系（CO2/As4S4，NW==NUM_WANN 不触发改写）不受影响，回归会给出
"正常案例通过、目标案例全灭"的结果。
**修法方向**：froz_max 必须随通道走——ESTIMATE_AUTO_WINDOW 按自旋各算一个
（ELIM 拆自旋即可），MLWF_WIN_SET_NUM_WANN 在改 num_wann 的同时改写该通道的
dis_froz_max；或对少数自旋通道在 spin 拷贝里降 froz_max 至其自身带隙中点。

### 2.2 BLOCKER-2：投影回退分支造成 .win 与库实参 num_wann 不一致 → w90 库内内存越界

**触发条件：wannier90.win 含 `begin projections` 且 NW<NUM_WANN——即轨迹模式
（EOS_EVERY_N>0）第 ≥2 步的磁性体系（LEOS_GUIDE 于首步末置位，eos.F:686；
第 2 步 MLWF_WANNIER90_WIN 写入投影块，mlwf.F:1881-1893），或 legacy 模式用户
.win 自带投影块。**

证据链：

1. `MLWF_WIN_SET_NUM_WANN` 见 'begin projections' 直接 RETURN 不改写
   （mlwf.F:1953-1957），注释声称"回退到旧的带幽灵行为"；但**调用方不知道回退发生**，
   仍以 NW 调 wannier_run 并把数组切片到 NW（mlwf.F:2018-2031 无状态回传）。
2. wannier90 3.1.0 库模式的真实语义（关键事实，已逐行核对官方源）：
   - `wannier_run` 把 `num_bands_loc` 拷入模块变量（wannier_lib.F90:308），
     **但从不拷 `num_wann_loc`**——它只用于哑元数组的形状声明
     （wannier_lib.F90:265-271）；
   - 随后 `call param_read()`（wannier_lib.F90:324）从 .win 读 num_wann，
     **该读取没有 library 守卫**（parameters.F90:633-636；对照 num_bands 在
     653-654 有 "Ignoring <num_bands> in input file"）。
   → 库内部 num_wann = .win 值 = NUM_WANN，而全部哑元按 NW 声明。
3. 后果（均在 master rank 内）：`a_matrix = a_matrix_loc`
   （wannier_lib.F90:341，左 (nb,NUM_WANN,nk) 右 (nb,NW,nk)，realloc-LHS 把
   内部数组缩成 NW 列，dis_main 仍按 num_wann=NUM_WANN 列读 → 越界读）；
   `u_matrix_loc = u_matrix`（406，往 NW²·nk 的 copy-out 临时数组写
   NUM_WANN²·nk 个元素 → **堆越界写**）；`wann_centres_loc = wannier_centres`
   （422）同理。生产版无边界检查 → 内存踩踏，崩溃或悄悄给错 U。
4. 侧注：投影块本身有 NUM_WANN 行（mlwf.F:1887 循环到 MLWF%NUM_WANN），与
   .win 的 num_wann=NUM_WANN 数目自洽，param_read 不会先报数目错——直接进入 3 的
   越界路径。

**修法方向**：让 MLWF_WIN_SET_NUM_WANN 返回是否改写成功（LOGICAL OUT），失败时
调用方必须把 NW 回退为 MLWF%NUM_WANN（真正的"旧行为回退"）；或投影块行数按通道
写 NW 行并同步改写 num_wann（还需处理 BLOCKER-1 的窗口）。**任何时候必须保证
".win num_wann == wannier_run 实参 num_wann"这一 w90 3.1 库模式的硬不变量。**

### 2.3 wannier_run 非连续切片的 copy-in/copy-out（评审点 ⑴）— 通过

mlwf.F 顶部有显式 INTERFACE（223-252）：`M_matrix_` (nb,nb,nnb,nk) IN、
`A_matrix_` (nb,NW,nk) IN、`U_matrix_` (NW,NW,nk) OUT、`U_matrix_opt_` (nb,NW,nk)
OUT、`wann_centres_`(3,NW)/`wann_spreads_`(NW)/`spread_`(3) OUT。
调用切片形状全部与哑元逐维相等（2029-2031）；非连续段（A 的 1:NW 列、U 的
1:NW×1:NW）由编译器生成连续临时数组，IN 只 copy-in，OUT 保证 copy-out 写回
切片区域；`M_matrix(:,:,:,:,ISP)` 是前导整维切片、天然连续，0.9 GB 无复制。
INTENT(OUT) 临时数组入口未定义 → w90 出口对 u_matrix_loc/centres/spreads 全量
赋值（wannier_lib.F90:406-430），无脏数据。**幽灵尾部安全的根据**：U/U_opt/centres/
spreads 分配后全部清零（mlwf.F:1988-1989），copy-out 只覆盖 1:NW 区域，广播是
M_sum_z8/M_sum_d 求和（2038-2043）而非 master 覆盖——非 master 全 0 + master
尾部 0 → 和仍为 0 ✓。

### 2.4 MLWF_WIN_SET_NUM_WANN 健壮性（评审点 ⑶）

- **MAJOR-3（静默截断改写）**：读取循环在 `N>=MAXLINES` 时把**刚读出的那一行丢弃
  并停止**（mlwf.F:1945-1950），随后 REPLACE 重写只写前 20000 行 → 超长 .win 被
  静默截尾；行长 >255 列同样被 `LINELEN=255` 截断后写回（1935-1936, 1964）。
  LEOS 自动生成的 .win 约数百行、行宽窄，不触发；**legacy 模式（用户手写 .win，
  正是师弟场景）有实际风险**。建议：超限时放弃改写并回退 NW（与 BLOCKER-2 的
  状态回传一并解决），而不是截断重写。
- begin projections 回退分支的 RETURN 时机：只读不写、CLOSE 已执行（1951），
  原文件保持原样 ✓（但见 BLOCKER-2——文件没坏，调用方语义坏了）。
- `IF (IO%IU6<0) RETURN`（1940）：正确但为死代码——唯一调用点在 RUN 的
  IU6>=0 块内（1993→2019）。.win 只在 master 工作目录读写 ✓。
- **MINOR-3**：`INDEX(ADJUSTL(LINES(I)),'num_wann')==1`（1961）会命中任何以
  num_wann 开头的行（含用户注释）；w90 关键词表中无其它 num_wann* 前缀词，
  实害有限。
- **MINOR-4**：硬编码 UNIT=99（1943/1959）与 COPYFILE(99/98)、eig/mmn 写共用；
  时序上串行无重叠（COPYFILE 先 CLOSE 再进入本例程），但新代码建议用 NEWUNIT=。
  改写 OPEN（1959）无 IOSTAT，目录只读等罕见场景直接运行时错误。

### 2.5 eos.F 收窄完备性（评审点 ⑷）— 通过

- `WF_OCC/WF_I1/WF_I2/WF_E1/WF_E2` 分配后全部清零——**v3base 就有**
  （eos.F.v3base:266-270；当前 eos.F:299-303），EOS_O_MAX 两个仍遍历到
  MLWF%NUM_WANN 的循环（506/528）读到的幽灵条目是 0，
  `WF_I2==0 /= I` 与 `WF_E2<=1e-12` 双保险 CYCLE → **无未初始化读取**。
- WF_OCC 幽灵列保持 0 的依据成立：累加循环只到 NWS（335-341）+ U 幽灵列为 0
  （见 2.3）→ `NELECT_WANN=SUM(WF_OCC)` 不被幽灵污染 ✓。
- `EOS_PROJ_WEIGHTS` 仍按 NUM_WANN 全宽计算（776/794/817）：UTOT 幽灵列=0 →
  CW=0 → PW 幽灵行=0，只浪费少量 FLOP，#WP 只打印 1:NWS ✓。
- EOS_GUIDE_CENTRES：首分配清零（679-680），每步只覆盖 1:NWS（683）→ 幽灵槽恒 0。
  这些 0 会作为 `c=0,0,0:s` 投影写进下一步 .win（mlwf.F:1887-1891 循环全宽）——
  该路径已被 BLOCKER-2 覆盖，此处记 **MINOR-5**（修 BLOCKER-2 时应把投影块行数
  与通道 num_wann 对齐）。

### 2.6 轨迹模式校验链（评审点 ⑸）

- **MAJOR-1（增长逻辑死代码 + 陈旧上界）**：`EOS_PREPARE_MLWF` 只在
  `NSTEP == 1` 被调用（main.F:3343），eos.F:199-208 的
  `NUM_WANN_SPIN(ISP)=MAX(旧值,NREQ_SPIN(ISP))`"逐步单调抬升"逻辑在第 ≥2 步
  **永远不会执行**（首步执行时刚赋完值，MAX 是恒等）。后果：轨迹中某通道占据数
  增长（如 CH4 活化中的自旋转移）时，该通道仍按首步 NW 跑 wannier90 →
  占据带覆盖不全 → 丢电子 → abort 分支。v3base 的共享 NUM_WANN 同样冻结在首步
  （同类缺陷），但共享上界更宽松：**per-spin 收紧使这一暴露面变大**，而注释声称
  的保护并不存在。修法：把 PREPARE（或其重算部分）挪进 EOS_DO_THIS_STEP 的
  每步路径；注意 NUM_WANN（数组尺寸）仍需保持首步 max 或同步报错。
- 首步 auto / 后续步语义：模块 NUM_WANN/NUM_WANN_SPIN 是 SAVE，首步赋值后各步
  MLWF_MAIN 经 519-551 重新装进 P_MLWF；MLWF_FREE（main.F:3348）每分析步重建
  数组，尺寸恒为首步 NUM_WANN → 尺寸链自洽 ✓（陈旧性即 MAJOR-1）。
- 显式 INCAR NUM_WANN 分支：NUM_WANN_SPIN 保持 0 → MLWF_MAIN:548 条件不成立 →
  两通道共用显式值（旧语义），提示语打印一次（eos.F:206-207）✓。
- 全极化边界（NREQ_SPIN(2)=0）：`ALL(NUM_WANN_SPIN(1:2)>0)` 为假 → 自动回退共享
  值，不炸 ✓。
- ISPIN=1 与非 LEOS：NW==NUM_WANN → 不改写 .win、切片 1:NW 即全宽，语义与字节级
  行为同基线 ✓（切片可能引入 copy-in/out 的纯性能开销）。

---

## 3 改动 3：任务一诊断版 EOSCAR（eos.F abort 分支 + 正常路径 #NL）

### 3.1 Fortran 侧 — 通过

- **abort 点变量定义完备**：IU/I（229-259 声明）、NLOST（305 置 0、323 累加）、
  LOST_REC/LOST_OCC（258-259 固定 32 条，324-329 在损失扫描循环内填充，
  该循环 305-344 先于 abort 判断 358）、LOST_NE（306/322）、NELECT_WANN/ERR
  （346-347）、EOS_STEP_COUNT/EOS_EVERY_N（模块 33-34）——全部在使用前有值 ✓。
- **IONODE guard**：整个诊断写块在 `IF (IO%IU6 >= 0)`（359-396）内；
  EOS_STEP_COUNT 只在 master 递增（368、594），且它的全部读者都在 IU6 guard 内
  → 各 rank 视图一致 ✓。
- **提前 RETURN 的并行一致性**：所有 rank 用同一份广播数据（U/U_opt/lwindow 来自
  RUN 的 M_sum/M_bcast 2038-2043）串行重算 WF_OCC → NELECT_ERR 逐位一致 →
  分支全员同向；abort RETURN（398）发生在唯一的后续集体调用
  `M_sum_d(COMM,PW)`（841，EOS_PROJ_WEIGHTS 内）之前且全员一起跳过 → 无死锁 ✓。
  PW 在 abort 路径尚未分配，DEALLOCATE 列表（397）不含 PW ✓。
- **REPLACE/APPEND 交错**：abort 与正常分支共用 EOS_STEP_COUNT；首个分析事件
  REPLACE（369/595），其后 APPEND；`# STEP` 行在"多步或 EOS_EVERY_N>0"时写
  （374/600）——abort 步与正常步任意交错时步号连续、文件结构一致 ✓。
  单点(EVERY_N=0)非轨迹 abort → 纯诊断文件、无 STEP 行 ✓。
- 正常路径 #NL 块（613-621）插在 `# EOS_NELECT_ABORT` 之后、逐离子表头（622）
  之前；NLOST=0 时不写 → 正常无损案例 EOSCAR 逐字节不变 ✓。
- **MINOR-7**：APPEND 分支 `STATUS='OLD'` 无 IOSTAT（372/598，后者系基线遗留）——
  运行中 EOSCAR 被外部删除时直接运行时错误；低概率。

### 3.2 pyleos parsers.py（评审点 ⑶ 前缀逐字符核对）

- `'# LEOS_ABORTED'` raise 位置正确：`parse_eoscar` 的 '#' 分支第一个判断
  （parsers.py:150-157），在构建任何数据前抛 ValueError，消息给出 #NL 的含义 ✓。
- `'#NL'` 行：`startswith('# N_LOST_STATES')` 在第 2 字符（'N' vs ' '）即失配；
  `#WP`/`#WF` 前缀不匹配 → 落入通用注释分支 continue ✓。
- `'# N_LOST_SHOWN'`：与 `'# N_LOST_STATES'` 公共前缀到 `'# N_LOST_S'`，第 11 个
  字符 'H'≠'T' → False；与 `'# N_LOST_E'` 在第 10 个字符 'S'≠'E' → False；
  也不含 `NELECT(Wannier)`/`|dNe|` 子串 → 安全落入注释分支 ✓。
- `'# NELECT(input)'` 仅存在于 abort 文件且位于 raise 之后，永不被解析 ✓。
- `'# STEP = n'`：通用注释分支 ✓。
- **MAJOR-2（轨迹文件被单帧 abort 毒化）**：EOS_EVERY_N>0 时 EOSCAR 是多步追加
  文件；只要**任何一帧**走了 abort 分支，文件中就有 '# LEOS_ABORTED' 行，
  `parse_eoscar` 对**整个文件** raise（151-157 不区分步块）——即使末帧正常也无法
  解析。v3base 行为是 abort 帧什么都不写（文件仍可解析）。属新功能与既有解析器的
  边界冲突；修法：parser 按 '# STEP' 分块、只对目标块（默认最后一块）判定
  ABORTED，或 VASP 侧给 abort 块加块级结束标记。（注：多步文件的逐离子行跨步
  累积是 v3 既有的 parser 限制，不在本次增量内，但同一处修复可一并解决。）

---

## 4 改动 4：A2 非 master 提前释放 M_matrix（mlwf.F:1311-1319）— 通过

全文件 + 全源码树消费点核查（`grep '%M_matrix'` 于 patched mlwf.F、eos.F 与
vasp_src 全部 .F）：

- 释放点在 `M_sum_z8(COMM_KINTER)` 归约之后（1311→1316-1319），此前的计算不受影响。
- `MLWF_WANNIER90_WRITE_MMN`：入口 `IF (IO%IU6<0) RETURN`（1346）先于唯一的
  M 引用（1378）→ 非 master 不触碰 ✓（调用点 MLWF_MAIN:565）。
- `MLWF_WANNIER90_RUN`：wannier_run 及其 M 实参整体在 `IF (IO%IU6>=0)` 内
  （1993-2034）；非 master 只执行分配/清零/广播与守卫式 DEALLOCATE（2050）✓。
- SPN（1396-…）、WANPROJ（2132-…）、EIG、UNK/XSF、WAVECAR、MLWF_WRITE_FINAL、
  eos.F（零处引用）、重启/READ_AMN 路径：均无 M_matrix 引用。
- 全树唯一外部引用 elphon.F:3083 只查 `A_matrix`（见改动 5），不查 M ✓。
- 轨迹再进入：MMN 重分配自带 `IF (ALLOCATED) DEALLOCATE`（1261-1263），且每步
  MLWF_FREE（main.F:3348 → mlwf.F:673）守卫式释放 → 重复步安全 ✓。
- master 语义位一致性：IU6>=0 ⇔ COMM 全局 IONODE ⇔ SCDM 的 DO_WORK
  （NODE_ME==1）为同一 rank——释放边界与消费边界重合 ✓。

---

## 5 改动 5：A3 wannier_run 后释放 M/A（mlwf.F:2045-2051）

### 5.1 通过项

- RUN 与 TRAFO 互斥：MLWF_MAIN 的 IF/ELSE（624-640）——SCDM-only / TRAFO 路径
  （`MLWF_TRAFO_RUN`:2061，其 `U_matrix_opt=A_matrix` 在 2122）**不经过 RUN**，
  不会读到已释放的 A ✓。
- RUN 之后的代码：WANPROJ/XSF/WAVECAR/WRITE_FINAL（642-656）与 eos.F 全部只用
  U/U_opt/lwindow/centres/spreads/eigenvalues（eigenvalues 未被 A3 释放）✓。
- 轨迹重分配两处 `ALLOCATE(MLWF%A_matrix)`（595 SCDM 路径、627 PROJ_UNKNOWN
  identity）本身**无 ALLOCATED 守卫**，但每个分析步先走 MLWF_FREE
  （main.F:3348 → mlwf.F:674 守卫式释放）→ 无 already-allocated 崩溃；A3 使
  595 在"上一步走了 RUN 路径"时也安全（双保险）✓。
- main.F:5766 的另一处 MLWF_MAIN 在 GW/RPA 分支（LCHI），与 LEOS 的 ab-initio
  分支互斥，无叠加调用 ✓。

### 5.2 MAJOR-4：ELPH(METHOD_LOCPROJ) + LWANNIER90_RUN 组合回归

elphon.F:3083-3084：`IF (.NOT. ALLOCATED(MLWF_T%A_matrix) ...) vtutor%error
("ELPH_INIT_WANNIER: Projections not set up properly!")`。
`ELPH_INTERPOL%INIT(..., MLWF_GLOBAL, ...)` 在 main.F:3367-3369 紧跟 MLWF 块执行。
v3base：RUN 结束后 A 仍存活 → 该组合可用；A3 后：A 已释放 → 干净报错退出。
LEOS 作业不用 ELPH，但这是 stock 功能组合的行为回归（干净错误、无损坏/挂死，
故 MAJOR 不是 BLOCKER）。修法：A 的释放加条件
`IF (.NOT. ELPH 需要)` 或挪到 main.F LEOS 专属路径；最小改动是仅在 LEOS 时释放 A
（M 的释放无此顾虑，elphon 不读 M）。

---

## 6 改动外但与被审区强相关的存量问题（不计入本次增量计数，建议记入负债）

1. **RRQR pivot 不跨自旋复位（stock）**：`RANK_REVEAL_QR%PIVOT=0` 只在 INIT
   （mathtools.F:266-267），CALC 不清（325-332）。ZGEQP3 的 JPVT 语义是
   "非零=固定列"——ISPIN=2 的第二个自旋进 CALC 时 JPVT 全非零 → **无枢轴选择**，
   出参恒为 1..N → spin2 锚点"pivot"是网格前 NUM_WANN 个点。SVD 仍产出合法
   正交 T、w90 能从劣质初猜收敛，故一直没炸，但 SCDM 质量（尤其 SCDM_ONLY 快速
   路径与 warm start）受损。A1 逐位保留了该行为（bug-for-bug 等价）。
   修一行：每自旋 CALC 前 `RANK_REVEAL_QR%PIVOT=0`（注意会改变 spin2 数值结果，
   需与回归基线协调）。
2. **warm-start 引导中心的自旋分辨失效（v3base:1860）**：
   `IF (INDEX(MLWF%SEED_NAME,'.2')>0) NI=2`——SEED_NAME 是 LEN=9 的 'wannier90'，
   永不含 '.2'（且装不下）→ 投影块永远写 spin1 的中心，spin2 被 spin1 引导。
   v3 已有，与本次增量无关，修 BLOCKER-2 时顺带处理。
3. stock SCDM/UNK 路径对 NCORE>1 的隐含不自洽（见 1.3），上游问题，只诊断。

---

## 7 结论与处置建议

| 级别 | 编号 | 一句话 | 位置 |
|---|---|---|---|
| BLOCKER | B1 | 幽灵修复 + 共享 AUTO_WINDOW：少数自旋冻结态数 > per-spin num_wann，w90 dis_windows io_error 硬停，master 死、其余 rank 挂在 M_sum_z8 —— 一切"目标体系"（VIRCEF/Ce）单点即触发 | mlwf.F:923-946,1868,2018-2031;w90 disentangle.F90:784-791,io.F90:307-362;eos.F:150 |
| BLOCKER | B2 | 投影回退分支不回退 NW：.win num_wann(=max) 与库实参 NW 不一致，而 w90 3.1 库模式的 num_wann 以 .win 为准（parameters.F90:633-636 无 library 守卫）→ 库内非一致形状赋值堆越界（轨迹第 ≥2 步磁性体系 / legacy 投影 .win） | mlwf.F:1953-1957,2018-2031;wannier_lib.F90:308,324,341,406,422 |
| MAJOR | M1 | EOS_PREPARE_MLWF 仅 NSTEP==1 调用 → per-spin"单调抬升"是死代码，轨迹中通道占据增长时按陈旧 NW 跑丢电子（较 v3base 共享上界更易触发） | main.F:3343;eos.F:199-208 |
| MAJOR | M2 | 轨迹多步 EOSCAR 只要有一帧 abort 就整文件被 parse_eoscar raise（不分步块），正常帧也读不出 | eos.F:368-395;parsers.py:150-157 |
| MAJOR | M3 | MLWF_WIN_SET_NUM_WANN 改写路径对 >20000 行 / >255 列的 .win 静默截断重写（legacy 手写 win 风险） | mlwf.F:1934-1950,1959-1966 |
| MAJOR | M4 | A3 释放 A_matrix 使 ELPH(METHOD_LOCPROJ)+LWANNIER90_RUN 组合从可用变报错（stock 功能回归，干净错误） | mlwf.F:2051;elphon.F:3083;main.F:3367 |
| MINOR | 1-9 | WVFN_PIVOT/PIVOT_IDX 配对无防御；NP==NPW 无断言；非共线锚点范数 ULP 变化；UNIT 99 硬编码+改写 OPEN 无 IOSTAT；num_wann 前缀 INDEX 误伤面;幽灵槽 c=0,0,0 投影；EOSCAR APPEND 无 IOSTAT；SET_NUM_WANN 的 IU6 死分支/RUN 内冗余 guard；WVFN_BANDS_T 跨步写性能 | 见各节 |

**通过项汇总**：A1 的 MPI 次序/数值等价（共线逐位）/索引空间（NCORE=1 平台成立并
确认基线同风险）/分配生命周期；A2 全消费点核查通过；A3 除 M4 外通过；改动 3 的
Fortran 侧全通过、parser 前缀逐字符核对通过；补丁完整性以 apply 级逐字节比对通过
（.objects/main.F 语义零变化，仅 mlwf.F/eos.F 扩展；patch/v3/eos.F 同步）。

**合入建议**：A1、A2、任务一（Fortran 侧）可保留；A3 需按 M4 缩小释放条件；
**幽灵修复必须先解决 B1（per-spin dis_froz_max）与 B2（回退语义闭环）再进编译
回归**，否则回归会在 VIRCEF/Ce 上以 wannier90 硬停/越界收场；M1 建议与幽灵修复
同批修（把 PREPARE 的重算挪进每步路径），M2 在 pyleos 侧按 STEP 分块。
