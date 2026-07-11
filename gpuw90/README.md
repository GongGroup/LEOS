# gpuw90 — 运行说明 (Phase 1: wannierise + Phase 2: disentangle)

独立 GPU Wannier-化后处理工具原型。读 wannier90 标准格式(`.mmn`/`.amn`/`.eig`/`.win`),用 PyTorch
批量实现解纠缠(Souza–Marzari–Vanderbilt)与 Marzari–Vanderbilt 幺正最小化,输出 Wannier 中心与
spread。架构与设计见 `DESIGN.md`。

- 孤立带组(`num_bands == num_wann`):直接 `wannierise()`(Phase 1)。
- 纠缠带组(`num_bands > num_wann`,金属/含导带):先 `disentangle()` 再 `wannierise()`(Phase 2)。
  直接把全带 M/A 喂给 wannierise 会得到静默错误结果,`load_case`/`wannierise()` 均有守卫拦截。

CPU/CUDA 同码,本地用 CPU torch 即可复现全部金标。

## 1. 环境

只需 `torch` + `numpy`(CPU 版即可,数学与 CUDA 相同)。

```bash
# 方式 A: venv(推荐,隔离干净)
python3 -m venv venv
source venv/bin/activate
pip install --index-url https://download.pytorch.org/whl/cpu torch numpy

# 方式 B: --user 安装 CPU 版
pip install --user torch numpy --index-url https://download.pytorch.org/whl/cpu
```

> 注意:装 **CPU** 版 torch(`--index-url .../whl/cpu`),不要装 CUDA 版。开发期在本地 CPU 验证正确性;
> 上 sai-pod V100 时同一份代码传 `--device cuda` 即可,无需改代码。

## 2. 跑金标对照(正确性验收)

### 2a. 孤立带(Phase 1)

```bash
# GaAs example01(8 k 点)
python check_against_w90.py testdata/gaas --seed gaas

# Diamond example05(64 k 点,压测批量 k 路径)
python check_against_w90.py testdata/diamond --seed diamond
```

输出:kmesh 权重与 B1 残差、初态对照、逐 WF 中心/spread 偏差表、Ω 分解对照,以及
**最大中心偏差 + PASS/FAIL**(阈值 1e-4 Å)。预期两例均 PASS(实测最大中心偏差 GaAs 8.3e-7 Å、
Diamond 3.1e-9 Å)。

常用参数:
```
--device cpu|cuda     计算设备(默认 cpu)
--num-iter N          最大迭代数(默认 2000;良态体系 ~10 迭代收敛)
--step S              SD 步长 α(默认 1.0 = 经典 MV 步长)
--conv-tol T          收敛阈值(默认 1e-10,对齐 wannier90)
--conv-window W       连续满足次数(默认 3)
--verbose             打印逐迭代 Ω 与梯度 RMS
```

### 2b. 纠缠带(Phase 2,解纠缠)

两套带冻结窗的真实两自旋金标(sai-pod VASP 导出 + 站点 wannier90.x 串行参考):
`testdata/vircef`(180 带 → 142 WF)与 `testdata/li2mno3`(356 带 → 288 WF,LEOS 生产尺度)。
`.win/.eig/参考 .wout` 已入库;10~300 MB 的 `.mmn/.amn` 不入库,先拉一次:

```bash
cd testdata && sh fetch_entangled.sh && cd ..   # 需能 ssh sai-pod
```

再跑四个自旋通道(每条几秒~1 分钟,不含可选的 --wannierise):

```bash
python check_dis_against_w90.py testdata/vircef   --seed wannier90.1
python check_dis_against_w90.py testdata/vircef   --seed wannier90.2
python check_dis_against_w90.py testdata/li2mno3  --seed wannier90.1
python check_dis_against_w90.py testdata/li2mno3  --seed wannier90.2
```

两道硬门(全部 PASS,数字见 `DESIGN.md` §5b):

1. **Ω_I**:解纠缠的 gauge 不变输出,与参考 `Final Omega_I` 之差 ≤ 1e-4 Ų
   (实测 4 通道 3e-9~5e-10 Ų,打印的 8 位小数完全相同,收敛迭代数也逐一相同)。
2. **解纠缠初态**:子空间 + Löwdin 投影初值的中心/spread 对参考 "Initial State"
   ≤ 1e-4 Å(实测 ~8e-7 Å,已到 wout 打印精度地板)。

加 `--wannierise` 会继续跑幺正最小化并对比末态——**仅供参考**:站点参考跑是 CG、100 迭代
封顶且无收敛判据(四个参考末态均未收敛),末态差异是协议差不是正确性信号(详见 DESIGN §5b)。
窗口与 dis 参数自动读自 `.win`(缺省 = wannier90 默认值);`--dis-num-iter` 可覆盖。

## 3. 作为库调用

### 3a. 孤立带(nb == nw)

```python
import torch, io_w90
import wannierise as wz

case = io_w90.load_case("seed.win", "seed.mmn", "seed.amn", "seed")
#   若 .win 无显式 kpoints,可传 nnkp_path="seed.nnkp"

M0 = torch.from_numpy(case.M)          # [nk, nnb, nw, nw]
A  = torch.from_numpy(case.A).to(torch.complex128)
U0 = wz.lowdin_orthonormalize(A)       # 初始规范;帧间 warm start 时改喂上一帧的 U

res = wz.wannierise(
    M0, U0,
    torch.from_numpy(case.bvec),
    torch.from_numpy(case.wb),
    torch.from_numpy(case.kb_idx),
    num_iter=2000, step=1.0, device="cpu",   # V100 上改 device="cuda"
)
print(res.centres)   # [nw, 3] 笛卡尔 Å  —— 供 LEOS 氧化态管线消费
print(res.spreads)   # [nw]
print(res.omega)     # {"I","D","OD","total"}
```

### 3b. 纠缠带(nb > nw):先 disentangle 再 wannierise

```python
import torch, io_w90
import wannierise as wz
import disentangle as dz

case = io_w90.load_case("seed.win", "seed.mmn", "seed.amn", "seed",
                        eig_path="seed.eig", allow_entangled=True)
# 不加 allow_entangled 会直接 ValueError(防静默错误结果);纠缠时 .eig 必传。

dres = dz.disentangle_case(case, device="cpu")   # 窗口/迭代参数取自 .win,w90 默认值兜底
print(dres.omega_I, dres.converged, dres.n_iter)  # gauge 不变 spread,供回归比对

res = wz.wannierise(
    dres.M0_reduced, dres.U0,                     # 降维后的输入,wannierise 不感知 nb
    torch.from_numpy(case.bvec),
    torch.from_numpy(case.wb),
    torch.from_numpy(case.kb_idx),
    device="cpu",
)
```

需要细粒度控制(显式窗口、mixing、收敛参数)时直接调 `dz.disentangle(...)`,签名见
`disentangle.py` 或 DESIGN §7。两个热循环都支持 `check_every=N`(设备端收敛计数,
每 N 迭代才做一次 host 同步;GPU 上建议 10)。

## 4. 测试数据来源

- 孤立带两组:取自 wannier90 3.1 test-suite(`testw90_example01` GaAs、`testw90_example05`
  diamond),`*_ref.wout` 为其 `benchmark.out`,完全收敛金标终态。
- 纠缠带两组(Phase 2):sai-pod `~/LEOS/w90_offline/` 下的 VASP `LWRITE_MMN_AMN` 导出
  (`bench/li2mno3_x000_f000`、`offline_runs/vircef`,各两自旋);`wannier90.{1,2}_ref.wout`
  为同目录 `run_np1_s{1,2}/` 的站点串行 wannier90.x 跑(输入即符号链接同一套文件,同源可溯)。
  大输入(`.mmn/.amn`)不入 git,用 `testdata/fetch_entangled.sh` 重取。

## 4b. 生产尺度计时(CPU 基线 / CUDA 待测)

固定迭代数计时脚本(不做收敛早停,s/iter 可跨设备直比):

```bash
python bench_gpu.py testdata/li2mno3 --seed wannier90.1               # 本地 CPU 基线
python bench_gpu.py testdata/li2mno3 --seed wannier90.1 --device cuda # GPU(集群)
```

本地 CPU 基线(Apple M 系, torch 2.8 CPU, FP64)见 `DESIGN.md` §5c。

**CUDA 实测(等 sai-pod GPU 空档,不阻塞本阶段)**:`run_bench_v100.sbatch` 已备好——先在
V100 上复跑金标门,过了才计时。两个已探明的坑写在模板注释里:16V100 分区最少 4 GPU/作业
(单卡任务申请 4 用 1,flood-gpu 4h 队列);站点 mace 环境的 torch 2.9.1+cu128 **不含 sm_70**,
V100 必须用 cu124 系轮子(torch ≤ 2.6),模板内有 `get_arch_list()` assert 守卫。

## 5. 有限差分梯度自检(可选)

梯度符号/前因子用有限差分方向导数钉死(见 `DESIGN.md` §4.1)。自检脚本 `fd_check.py`
比较解析 `Re Tr[g† X]` 与数值 `[Ω(U e^{εX})−Ω(U e^{−εX})]/2ε`,预期 ratio ≈ +1.0000。
移植/改梯度后建议先跑它:`python fd_check.py`。

## 6. 目录

| 文件 | 作用 |
|---|---|
| `io_w90.py` | `.mmn/.amn/.eig/.win/.nnkp` 解析、b 向量与 shell 权重重建、`.wout` 参考解析、纠缠守卫 |
| `wannierise.py` | torch 批量 MV 最小化(Löwdin 初值、梯度、矩阵指数、SD 驱动;非方阵守卫) |
| `disentangle.py` | torch 批量 SMV 解纠缠(窗口/冻结、Z 矩阵不动点迭代、M0_reduced/U0 输出) |
| `check_against_w90.py` | 孤立带金标对照与逐 WF 偏差表 |
| `check_dis_against_w90.py` | 纠缠带金标对照(Ω_I + 解纠缠初态两道硬门) |
| `test_guards.py` | 纠缠守卫单测(`python test_guards.py` 或 pytest) |
| `bench_gpu.py` / `run_bench_v100.sbatch` | 固定迭代计时 + sai-pod V100 作业模板 |
| `fd_check.py` | 梯度符号/前因子的有限差分自检 |
| `DESIGN.md` | 架构、数值坑、disentangle 实现细节、金标数字、HBM 分块、LEOS 对接、待办 |
| `testdata/` | GaAs / Diamond(孤立带)+ vircef / li2mno3(纠缠带)金标算例 |
