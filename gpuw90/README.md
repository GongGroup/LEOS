# gpuw90

gpuw90 用 PyTorch 实现 wannier90 的两个计算步骤：解纠缠（Souza–Marzari–Vanderbilt）和
Marzari–Vanderbilt 局域化。它读 wannier90 标准格式的输入，输出 Wannier 中心和 spread。
CPU 和 CUDA 共用同一份代码。设计、算法和逐项验证结果见 DESIGN.md。

## 工具简介

gpuw90 在 VASP 之外离线运行，承担 wannier90.x 中解纠缠和局域化两步，不修改 wannier90 和 VASP。
调用顺序取决于带组类型：

| 带组类型 | 判据 | 调用顺序 |
|---|---|---|
| 孤立带 | `num_bands == num_wann` | `wannierise()` |
| 纠缠带（金属或含导带） | `num_bands > num_wann` | `disentangle()` → `wannierise()` |

**纠缠带组必须先做解纠缠。把全带的 M、A 直接交给 `wannierise()`，得到的中心是错的，而且不会报错。**
代码为此设了两道守卫，触发时都抛 `ValueError`：

- `io_w90.load_case()` 遇到 `num_bands > num_wann` 时，要求显式传 `allow_entangled=True`。
- `wannierise()` 拒绝非方阵的 `M0` 或 `U0`。

| 项目 | 内容 |
|---|---|
| 输入 | `seedname.mmn`、`seedname.amn`、`seedname.win`；纠缠带另需 `seedname.eig`；`.win` 没有显式 kpoints 块时另需 `seedname.nnkp` |
| 输出 | Python 对象：每个 WF 的中心（笛卡尔坐标，Å）、spread（Ų）和 Ω 分解。不写文件 |
| 依赖 | `torch`、`numpy` |

代码注释中的 Phase 1 指孤立带局域化，Phase 2 指解纠缠。

## 环境安装

只需要 `torch` 和 `numpy`。本地装 CPU 版 torch 就能跑全部验收，计算与 CUDA 版相同。

```bash
# 方式 A：虚拟环境（.venv/ 已在 .gitignore 中）
python3 -m venv .venv
source .venv/bin/activate
pip install --index-url https://download.pytorch.org/whl/cpu torch numpy

# 方式 B：装到用户目录
pip install --user torch numpy --index-url https://download.pytorch.org/whl/cpu
```

`--index-url .../whl/cpu` 指定 CPU 版，本地不需要 CUDA 版。

在 GPU 上运行时不用改代码：命令行加 `--device cuda`，库调用传 `device="cuda"`。
V100 对 torch 版本有要求，见下文"性能实测"中的"在 SAI 集群的 V100 上运行"。

## 运行验收

验收脚本把 gpuw90 的结果与 wannier90 参考输出（`testdata/<算例>/<seed>_ref.wout`）比较，
最后打印 `VERDICT: PASS` 或 `VERDICT: FAIL`，退出码分别为 0 和 1。
以下命令都在 `gpuw90/` 目录下运行。

### 孤立带：check_against_w90.py

两个孤立带算例都应 PASS。

```bash
python3 check_against_w90.py testdata/gaas --seed gaas          # GaAs，8 个 k 点
python3 check_against_w90.py testdata/diamond --seed diamond    # 金刚石，64 个 k 点
```

脚本依次打印：

1. k 网格的 shell 半径、shell 权重和 B1 残差；
2. 初态（Löwdin 正交化后、未优化）与参考 "Initial State" 的对照；
3. 逐 WF 的中心和 spread 偏差表；
4. Ω 各分量（I、D、OD、总和）的对照；
5. 最大中心偏差和判定。

判据是最大中心偏差 ≤ 1e-4 Å。WF 按最近中心配对，不按序号，所以允许近简并的 WF 交换顺序。
实测最大中心偏差为 GaAs 8.3e-7 Å、Diamond 3.1e-9 Å，完整数字见 DESIGN.md 的"正确性验证结果"。

该脚本只处理孤立带。纠缠带算例在加载时就会被守卫拦下。

| 参数 | 默认值 | 含义 |
|---|---|---|
| `datadir`（位置参数） | `testdata/gaas` | 算例目录 |
| `--seed` | `gaas` | 文件名前缀，必须与目录内的文件一致 |
| `--device` | `cpu` | `cpu` 或 `cuda` |
| `--num-iter` | 2000 | 最大迭代数；良态体系约 10 次迭代收敛 |
| `--step` | 1.0 | 最速下降步长系数 α；1.0 即经典 MV 步长 |
| `--conv-tol` | 1e-10 | 收敛阈值，与 wannier90 默认值相同 |
| `--conv-window` | 3 | 需要连续满足阈值的迭代次数 |
| `--verbose` | 关 | 打印前 3 次及每 50 次迭代的 Ω 和梯度 RMS |

### 纠缠带：check_dis_against_w90.py

纠缠带验收用两套带冻结窗的自旋极化算例，共 4 个自旋通道，都应 PASS：

- `testdata/vircef`：180 带 → 142 个 WF；
- `testdata/li2mno3`：356 带 → 288 个 WF，与 LEOS 生产计算的体系规模相同。

第一步，取回大输入文件。`.win`、`.eig` 和参考 `.wout` 已在 git 中。
`.mmn` 和 `.amn` 每个 10~300 MB，不入库，需要从 SAI 集群拉取一次：

```bash
cd testdata && sh fetch_entangled.sh && cd ..
```

该脚本通过 ssh 访问 SAI 集群。ssh 主机别名和集群上的路径写死在脚本里，运行前按自己的 ssh 配置核对。

第二步，逐个通道运行。不加 `--wannierise` 时，每个通道用时几秒到 1 分钟左右：

```bash
python3 check_dis_against_w90.py testdata/vircef  --seed wannier90.1
python3 check_dis_against_w90.py testdata/vircef  --seed wannier90.2
python3 check_dis_against_w90.py testdata/li2mno3 --seed wannier90.1
python3 check_dis_against_w90.py testdata/li2mno3 --seed wannier90.2
```

外窗、冻结窗和解纠缠参数自动从 `.win` 读取，`.win` 中没有的项用 wannier90 默认值。

判定只看两道门，两道都通过才是 PASS。逐通道的实测数字见 DESIGN.md 的"正确性验证结果"。

| 门 | 比较的量 | 阈值 |
|---|---|---|
| GATE 1 | 解纠缠后的 Ω_I 与参考 `Final Omega_I` 之差 | ≤ 1e-4 Ų |
| GATE 2 | 解纠缠后初态的 WF 中心与参考 "Initial State" 之差 | ≤ 1e-4 Å |

GATE 2 的初态是解纠缠子空间内的 Löwdin 投影规范，不经过任何最小化。
脚本还会打印 `|Omega_I(initial) - Omega_I(dis)|`。这两个量都与规范无关，差值应接近 0。

`--wannierise` 会继续做局域化并打印末态对照，这部分不参与判定。
原因是参考 `.wout` 的局域化段用 CG、固定 100 次迭代、没有收敛判据，四个参考末态都没有收敛。
末态的差异来自两边运行设置不同，不能用来判断实现对错。详见 DESIGN.md 的"正确性验证结果"。

| 参数 | 默认值 | 含义 |
|---|---|---|
| `datadir`（位置参数） | `testdata/vircef` | 算例目录 |
| `--seed` | `wannier90.1` | 文件名前缀；`.1`、`.2` 分别是两个自旋 |
| `--device` | `cpu` | `cpu` 或 `cuda` |
| `--dis-num-iter` | `.win` 中的 `dis_num_iter`，没有则 200 | 解纠缠最大迭代数 |
| `--wannierise` | 关 | 继续做局域化并打印末态对照，仅供参考；288 个 WF 时在 CPU 上较慢 |
| `--num-iter` | 3000 | 局域化最大迭代数，只在 `--wannierise` 时使用 |
| `--step` | 1.0 | 局域化步长系数 α |
| `--omega-tol` | 1e-4 | GATE 1 阈值（Ų） |
| `--tol-centre` | 1e-4 | GATE 2 阈值（Å） |
| `--verbose` | 关 | 打印逐迭代信息 |

### 守卫单测与梯度自检

仓库 CI 不运行 gpuw90 的测试。改动 gpuw90 后手动运行：

```bash
python3 test_guards.py      # 也可以 python3 -m pytest test_guards.py -q
python3 fd_check.py         # 可选参数 [testdata_dir] [seed]，默认 testdata/gaas gaas
```

- `test_guards.py` 有 5 个单测，覆盖纠缠带守卫：`load_case()` 默认拒绝纠缠带，纠缠带缺 `eig_path` 时报错，
  正常加载路径可用，`wannierise()` 拒绝非方阵输入和 nw 不一致的输入。
- `fd_check.py` 在随机反厄米方向 X 上，比较解析方向导数 `Re Tr[g† X]` 与中心差分
  `[Ω(U e^{εX}) − Ω(U e^{−εX})] / 2ε`。ε 取 1e-4、1e-5、1e-6，预期 ratio ≈ +1.0000。
  改动梯度或移植代码后先跑它。原因见 DESIGN.md 的"数值问题与处理"。

## 作为库调用

以下示例假定当前目录是 `gpuw90/`，模块按同目录导入。

### 孤立带（nb == nw）

```python
import torch
import io_w90
import wannierise as wz

case = io_w90.load_case("seed.win", "seed.mmn", "seed.amn", "seed")
# .win 没有显式 kpoints 块时，加 nnkp_path="seed.nnkp"

M0 = torch.from_numpy(case.M)                       # [nk, nnb, nw, nw]
A = torch.from_numpy(case.A).to(torch.complex128)   # [nk, nw, nw]
U0 = wz.lowdin_orthonormalize(A)                    # 初始规范

res = wz.wannierise(
    M0, U0,
    torch.from_numpy(case.bvec),
    torch.from_numpy(case.wb),
    torch.from_numpy(case.kb_idx),
    num_iter=2000, step=1.0, device="cpu",          # GPU 上改为 device="cuda"
)
print(res.centres)   # [nw, 3]，笛卡尔坐标，Å
print(res.spreads)   # [nw]，Å^2
print(res.omega)     # 键为 "I"、"OD"、"D"、"total"
```

### 纠缠带（nb > nw）

```python
import torch
import io_w90
import wannierise as wz
import disentangle as dz

case = io_w90.load_case("seed.win", "seed.mmn", "seed.amn", "seed",
                        eig_path="seed.eig", allow_entangled=True)
# 不传 allow_entangled=True 会抛 ValueError；纠缠带必须给 eig_path

dres = dz.disentangle_case(case, device="cpu")   # 窗口和迭代参数取自 .win，缺省用 wannier90 默认值
print(dres.omega_I, dres.converged, dres.n_iter)  # omega_I 与规范无关，可用于回归比较

res = wz.wannierise(
    dres.M0_reduced, dres.U0,                     # 已降到 nw 维
    torch.from_numpy(case.bvec),
    torch.from_numpy(case.wb),
    torch.from_numpy(case.kb_idx),
    device="cpu",
)
```

需要显式指定窗口、mixing 或收敛参数时，直接调用 `dz.disentangle(...)`。
它的签名和 `.win` 关键字的对应关系见 DESIGN.md 的"算法"。

### 主要函数与返回值

```python
io_w90.load_case(win_path, mmn_path, amn_path, seedname="case",
                 nnkp_path=None, eig_path=None, allow_entangled=False)     # -> W90Case

disentangle.disentangle_case(case, device="cpu", num_iter=None, mix_ratio=None,
                             conv_tol=None, conv_window=None,
                             check_every=1, verbose=False)                  # -> DisentangleResult

wannierise.wannierise(M0, U0, bvec, wb, kb_idx, num_iter=1000, step=1.0,
                      conv_tol=1e-10, conv_window=3, device="cpu",
                      dtype=torch.complex128, verbose=False,
                      check_every=1)                                        # -> WannieriseResult
```

库函数 `wannierise()` 的 `num_iter` 默认是 1000；`check_against_w90.py` 的默认值是 2000。

| 返回类型 | 字段 |
|---|---|
| `W90Case` | `seedname`、`M`、`A`、`bvec`、`wb`、`kb_idx`、`kpoints`、`real_lattice`、`recip_lattice`、`kmesh_info`、`win`（解析后的 `.win` 关键字）、`eig`；以及 `num_bands`、`num_wann`、`num_kpts`、`nntot` |
| `DisentangleResult` | `M0_reduced [nk,nnb,nw,nw]`、`U0 [nk,nw,nw]`、`V [nk,nb,nw]`、`omega_I`（Ų）、`history`、`n_iter`、`converged`、`ndimwin [nk]`、`ndimfroz [nk]` |
| `WannieriseResult` | `U [nk,nw,nw]`、`centres [nw,3]`（Å）、`spreads [nw]`、`omega`、`history`、`n_iter`、`converged` |

`history` 是逐迭代记录的列表：`DisentangleResult` 中每项为 `(iter, Ω_I(i−1), Ω_I(i), delta)`，`WannieriseResult` 中每项为 `(iter, Ω_total, rms_grad)`。

`wannierise()` 和 `disentangle()` 都有参数 `check_every=N`。
收敛计数保存在设备上，每 N 次迭代才向主机读回一次收敛标志。
GPU 上建议取 10。N > 1 时，循环最多比实际收敛多跑 N − 1 次迭代。

## 测试数据

`testdata/` 下有四套验收算例，每套都附 wannier90 参考输出 `<seed>_ref.wout`。

| 目录 | 类型 | nb / nw | k 点 | 输入来源 | 参考输出 |
|---|---|---|---|---|---|
| `gaas/` | 孤立带 | 4 / 4 | 8 | wannier90 3.1 test-suite 的 `testw90_example01` | 该例的 `benchmark.out`，完全收敛 |
| `diamond/` | 孤立带 | 4 / 4 | 64 | wannier90 3.1 test-suite 的 `testw90_example05` | 该例的 `benchmark.out`，完全收敛 |
| `vircef/` | 纠缠带，两自旋 | 180 / 142 | 8 | SAI 集群 `~/LEOS/w90_offline/offline_runs/vircef` 下的 VASP `LWRITE_MMN_AMN` 导出 | 同目录 `run_np1_s{1,2}/` 中 wannier90.x 的串行运行 |
| `li2mno3/` | 纠缠带，两自旋 | 356 / 288 | 8 | SAI 集群 `~/LEOS/w90_offline/bench/li2mno3_x000_f000` 下的 VASP `LWRITE_MMN_AMN` 导出 | 同目录 `run_np1_s{1,2}/` 中 wannier90.x 的串行运行 |

| 目录 | git 中的文件 | 需取回的文件 |
|---|---|---|
| `gaas/` | `gaas.{win,mmn,amn}`、`gaas_ref.wout` | 无 |
| `diamond/` | `diamond.{win,mmn,amn}`、`diamond_ref.wout` | 无 |
| `vircef/` | `wannier90.{1,2}.{win,eig}`、`wannier90.{1,2}_ref.wout` | `wannier90.{1,2}.{mmn,amn}` |
| `li2mno3/` | `wannier90.{1,2}.{win,eig}`、`wannier90.{1,2}_ref.wout` | `wannier90.{1,2}.{mmn,amn}` |

- 纠缠带的参考运行用的是指向同一套文件的符号链接，所以 gpuw90 和参考读的是完全相同的输入。
- 需取回的文件由 `.gitignore` 排除，用 `testdata/fetch_entangled.sh` 下载。

## 性能实测

在 li2mno3 算例上，V100 每次迭代比集群上串行运行的 wannier90.x 快约 14 倍（解纠缠）和约 22 倍（局域化）。
下面的计时都是固定迭代数下的每迭代耗时（s/迭代），可以跨设备直接比较。

### V100（SAI 集群）

测量于 2026-07-11，SAI 集群作业 644120。算例为 li2mno3 生产帧：288 个 WF，356 带，8 个 k 点。

| 阶段 | V100 | 计算节点 torch-CPU | 集群 wannier90.x 串行 | V100 相对串行 |
|---|---|---|---|---|
| disentangle | 0.079 | 0.423 | ~1.10 | ~14× |
| wannierise | 0.041 | 0.303 | ~0.92 | ~22× |

单位：s/迭代。

- CUDA 上的验收同样 PASS：|ΔΩ_I| = 4.4e-9 Ų，初态最大中心偏差 8.3e-7 Å，与 CPU 结果在同一水平。
- 以 docs/RESULTS.md"效率"一节中的 554 秒/帧为基线（wannier90 段），折算后两自旋的 Wannier 化段约半分钟。
  这是按每迭代耗时折算的估计。

### 本地 CPU 基线

测量于 2026-07-11。环境为 Apple M 系列笔记本，torch 2.8.0 CPU 版，全程 FP64。
wannier90.x 串行一列取自参考输出 `testdata/{li2mno3,vircef}/wannier90.1_ref.wout` 中的计时（单进程）。

| 阶段 | 算例 | gpuw90 CPU（s/迭代） | 集群 wannier90.x 串行（s/迭代） | 备注 |
|---|---|---|---|---|
| disentangle | li2mno3，288 个 WF | 0.24（38 次迭代，9.1 s） | 1.10（38 次迭代，41.9 s） | 迭代数相同，4.6× |
| wannierise | li2mno3，288 个 WF | 0.196 | ~0.92（CG，100 次迭代约 92 s） | 每次迭代的工作量同型 |
| disentangle | vircef，142 个 WF | 0.07 | ~0.21（82 次迭代，17.1 s） | |
| wannierise | vircef，142 个 WF | 0.051 | ~0.18 | |
| `.mmn` 读入 | li2mno3（258 MB 文本） | 3.4 s | 4.7 s | 一次性开销 |

本地 CPU 比串行 wannier90.x 快，是因为 batched matmul 用上了多线程 BLAS。
所以衡量 GPU 的收益时要同时对照 torch-CPU 和串行 wannier90.x 两条基线，V100 表中两列都已列出。

### 复现计时

```bash
python3 bench_gpu.py testdata/li2mno3 --seed wannier90.1                 # 本地 CPU
python3 bench_gpu.py testdata/li2mno3 --seed wannier90.1 --device cuda   # GPU
```

`bench_gpu.py` 对一个自旋通道依次跑解纠缠和局域化，打印每个阶段的总时间和 s/迭代。

- 迭代数固定：脚本内部设 `conv_tol=0`，不会提前收敛。
- 脚本不判断正确性。计时前先用 `check_dis_against_w90.py` 通过验收。
- 末尾打印的 Ω_I 是固定迭代数下的值，不是收敛值。

| 参数 | 默认值 | 含义 |
|---|---|---|
| `datadir`（位置参数） | 必填 | 算例目录 |
| `--seed` | `wannier90.1` | 文件名前缀 |
| `--device` | `cpu` | `cpu` 或 `cuda` |
| `--dis-iters` | 50 | 解纠缠的固定迭代数 |
| `--wann-iters` | 100 | 局域化的固定迭代数 |
| `--check-every` | 10 | 热循环中读回收敛标志的间隔 |

### 在 SAI 集群的 V100 上运行

`run_bench_v100.sbatch` 是 V100 计时的作业模板。它按以下顺序执行：

1. 在 CUDA 上对 li2mno3 自旋 1 跑 `check_dis_against_w90.py`，失败即退出；
2. 用 CUDA 计时，`--dis-iters 50 --wann-iters 200`；
3. 在同一计算节点上用 CPU 计时，`--dis-iters 10 --wann-iters 20`。

- 代码路径默认 `$HOME/LEOS/gpuw90`（把本地 `gpuw90/` rsync 过去），数据路径默认
  `$HOME/LEOS/w90_offline/bench/li2mno3_x000_f000`。两者在脚本开头的 `GPUW90`、`DATA` 变量中修改。
- 脚本把 `run_np1_s{1,2}/` 下的参考 `.wout` 软链为数据目录中的 `<seed>_ref.wout`。
- 作业输出写到提交目录下的 `bench-out`。

集群环境的要求如下，模板注释里也有说明：

| 事项 | 要求 |
|---|---|
| 分区与 QOS | `16V100` 分区每个作业至少申请 4 块 GPU。模板申请 4 块、只用 1 块（`CUDA_VISIBLE_DEVICES=0`），QOS 为 `flood-gpu`，时限 4 小时 |
| 资源参数 | 不要加 `--cpus-per-task` 或 `--mem*`。站点钩子会自动取消带这些参数的 GPU 作业 |
| torch 版本 | V100 是 sm_70。站点 mace 环境中的 torch 2.9.1+cu128 不含 sm_70，运行时报 "no kernel image"。需要 cu124 系列的轮子（torch ≤ 2.6）。模板中有 `torch.cuda.get_arch_list()` 断言检查 |
| Python 环境 | mamba 环境 `gpuw90`（python 3.11），一次性安装命令写在模板注释中 |

## 文件清单

| 文件 | 作用 |
|---|---|
| `io_w90.py` | 解析 `.mmn`、`.amn`、`.eig`、`.win`、`.nnkp`；重建 b 向量和 shell 权重；解析参考 `.wout`；`load_case()` 含纠缠带守卫 |
| `wannierise.py` | 批量 MV 局域化：Löwdin 初值、Ω 分解、梯度、矩阵指数、最速下降循环；拒绝非方阵输入 |
| `disentangle.py` | 批量 SMV 解纠缠：窗口与冻结掩码、初始子空间、Z 矩阵不动点迭代；输出 `M0_reduced` 和 `U0` |
| `check_against_w90.py` | 孤立带验收 |
| `check_dis_against_w90.py` | 纠缠带验收（Ω_I 和解纠缠初态两道门） |
| `test_guards.py` | 纠缠带守卫单测 |
| `fd_check.py` | 梯度符号和前因子的有限差分自检 |
| `bench_gpu.py` | 固定迭代数的计时脚本 |
| `run_bench_v100.sbatch` | SAI 集群 V100 计时作业模板 |
| `DESIGN.md` | 设计、算法、数值问题、验证结果、显存估算、与 LEOS 的对接、已知限制 |
| `testdata/` | 四套验收算例，以及取回大输入的 `fetch_entangled.sh` |
