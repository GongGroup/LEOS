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
