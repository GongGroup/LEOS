# CEH4 — I4/mmm (ThCr2Si2 型) CeH4,0 GPa 优化胞

- 来源:Bi & Zurek, *Chem. Eur. J.* **27**, 14858 (2021), DOI 10.1002/chem.202102679,
  Supporting Information 第 S2 节 "Structural Information for Metal Tetrahydrides at 0 GPa"
  (SI PDF 经 Wayback Machine 快照 2025-07-10 获取)。
- 空间群:I4/mmm (No. 139)。惯用体心四方胞,Z=2,共 10 原子。
- 晶格参数(0 GPa,PBE 优化):a = b = 3.545 Å, c = 6.149 Å。
- Wyckoff 坐标:Ce 2a (0,0,0);H 4d (0,1/2,1/4);H 4e (0,0,z), z = 0.6273。
- 压力点:**0 GPa**(该相热力学稳定区在高压,Bi & Zurek 在 0 GPa 下做静态优化并分析成键;
  Comas-Vilà JCTC 2025 对照集(案例 16-23)即采用这套 0 GPa 结构 —— 其 SI 中 H-H 距离
  1.566 Å 与本结构 |1-2z|·c = 1.566 Å 完全一致,故取 0 GPa 以保证对照可比)。
- 自检:化学计量 Ce2H8 = 2×CeH4 ✓;H2 分子单元 H-H = 1.566 Å(= 文献值)✓;
  最近 Ce-H = 1.566 Å;元素顺序沿用 data/v3/comas_gap/ceh4/POSCAR 的既有顺序
  (Ce H),便于 POTCAR 拼装逻辑复用。
- 替换对象:data/v3/comas_gap/ceh4/POSCAR(登记在案的 0 GPa 近似相占位结构)。
