# MGXE — CsCl 型 (B2) MgXe,200 GPa

- 相与压力:Miao et al., *J. Am. Chem. Soc.* **137**, 14122 (2015), DOI 10.1021/jacs.5b08162:
  MgXe 在 ≥125 GPa 进入凸包,>100 GPa 稳定构型为 CsCl 型(Pm-3m);
  Comas-Vilà JCTC 2025(案例 33)在 **200 GPa**、Pm-3m 下做 EOS 分析,本结构与之对齐。
- 空间群:Pm-3m (No. 221)。Mg 1a (0,0,0),Xe 1b (1/2,1/2,1/2),Z=1,共 2 原子。
- 晶格常数推定(文献未给数值表,ACS SI 不可获取):
  V(MgXe,200 GPa) = V_Mg + V_Xe + ΔV
  其中 V_Mg(200 GPa) = 9.07 Å^3(bcc Mg Vinet EOS:V0=23.1398 Å^3,K0=28.9 GPa,
  K'=4.80;Stinton et al., PRB **90**, 134105 (2014));
  V_Xe(200 GPa) = 15.66 Å^3(hcp Xe Vinet EOS:V0=57.41 Å^3/原子,K0=4.887 GPa,
  K'=6.2955,适用区 10-260 GPa;Dewaele et al. 2012 拟合);
  ΔV(200 GPa) ≈ -1.9±0.3 Å^3,读自 Miao 论文 Figure 3b(ΔV=V_MgXe-V_Mg-V_Xe 曲线)。
  => V = 22.83 Å^3,a = 2.837 Å(不确定度约 ±0.04 Å,即 ~1.5%)。
- 自洽校验:论文正文"125 GPa 时 MgXe 体积比 Mg+Xe 之和小约 10%"
  => ΔV(125) ≈ -2.88 Å^3,与 Figure 3b 读数(约 -2.6 至 -2.9 Å^3)一致 ✓。
- 自检:计量 Mg1Xe1 ✓;Mg-Xe = 2.457 Å(= a·√3/2,八配位 CsCl 几何);
  Xe-Xe = a = 2.837 Å,与 200 GPa 单质 Xe 的最近邻 (~2.8 Å) 相当,
  符合论文"Xe-Xe 间存在强相互作用"的描述。
- 注意:a 为推定值而非文献转录值,是本包中唯一非直接转录的晶格常数;
  若后续拿到 Miao SI(ja5b08162_si_001.pdf)应以其表值替换。
