基准氧化物 LEOS 结果库
56 个物质/结构标识，72 个计算版本。

组织方式：物质 / 参数版本 / 结果和输入文件。目录索引.tsv列出各版本。
参数文件夹依次表示泛函、U、ISPIN、轨道数、解纠缠上限、局域化上限、ISMEAR、运行方式。
详细电子态、磁矩、参数和来源记录在parameters_and_provenance.json及参数说明.txt。
按现有data/v3的轻量结果方式保存EOSCAR、POSCAR、INCAR、KPOINTS、CONTCAR、较小win/wout及中心、氧化态、beta区间结果。
POTCAR按用户要求保留。缺失文件在input_files_present中明确标记，不补造。
离线GPUw90的DFT来源输入单列于source_DFT_inputs，不能误当作离线优化设置。
不上传OUTCAR、超过1 MiB的Wannier日志、逐步优化历史、完整beta逐点评估和大型波函数/密度/矩阵；原件保留在SAI。
失败尝试也按原状归档，不等于经过科学验收。原始文件哈希与来源可在文件清单中追溯。
当前仓库为私有；POTCAR的使用与分发遵循原赝势许可。
