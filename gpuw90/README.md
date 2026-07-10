# gpuw90 — 运行说明 (Phase 1: wannierise)

独立 GPU Wannier-化后处理工具原型。读 wannier90 标准格式(`.mmn`/`.amn`/`.win`),用 PyTorch
批量实现 Marzari–Vanderbilt 幺正最小化,输出 Wannier 中心与 spread。架构与设计见 `DESIGN.md`。

本阶段只做 **wannierise**(孤立带组,`num_bands == num_wann`)。CPU/CUDA 同码,本地用 CPU torch 即可复现。

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

## 3. 作为库调用

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

## 4. 测试数据来源

`testdata/` 内两组算例取自 sai-pod `~/apps/wannier90-3.1.0/test-suite/tests/`
(`testw90_example01` GaAs、`testw90_example05` diamond),`*_ref.wout` 为其 `benchmark.out`,
即完全收敛的金标终态(中心/spread/Ω 分解)。均为孤立带组、纯 wannierise 算例。

## 5. 有限差分梯度自检(可选)

梯度符号/前因子用有限差分方向导数钉死(见 `DESIGN.md` §4.1)。自检脚本 `fd_check.py`
比较解析 `Re Tr[g† X]` 与数值 `[Ω(U e^{εX})−Ω(U e^{−εX})]/2ε`,预期 ratio ≈ +1.0000。
移植/改梯度后建议先跑它:`python fd_check.py`。

## 6. 目录

| 文件 | 作用 |
|---|---|
| `io_w90.py` | `.mmn/.amn/.win/.nnkp` 解析、b 向量与 shell 权重重建、`.wout` 参考解析 |
| `wannierise.py` | torch 批量 MV 最小化(Löwdin 初值、梯度、矩阵指数、SD 驱动) |
| `check_against_w90.py` | 金标对照与逐 WF 偏差表 |
| `fd_check.py` | 梯度符号/前因子的有限差分自检 |
| `DESIGN.md` | 架构、数值坑、disentangle 接口预留、HBM 分块、LEOS 对接、下一阶段清单 |
| `testdata/` | GaAs / Diamond 金标算例 |
