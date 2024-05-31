2024-1-15比2023-9-25修改部分：
SegFinder.sh脚本中
（1） [--min_multi]... 修改为 [--min_rdrp_multi]... [--min_nordrp_multi]... （我们对rdrp序列的载量和同一簇中非rdrp的载量做了限定）
（2）由于上述参数改变导致min_multi=20 变为min_rdrp_multi=50  min_nordrp_multi=20
（3）由于上述参数改变导致 Rscript Cor_contigs_extract.R  ${cor} ${library_ID} ${min_multi} ${min_TPM}变为 	Rscript Cor_contigs_extract.R  ${cor} ${library_ID}  ${min_TPM} ${min_rdrp_multi} ${min_nordrp_multi}

3.coefficient-matrix.R
脚本3.coefficient-matrix.R有修改，优化了excel表格的读取输出问题，相关性计算速度问题，脚本直接替换即可


Cor_contigs_extract.R
（1）修改了一些小问题：NR_idetity变为NR_identity（单词拼错了）
（2）添加参数后的修改rdrp <- which(confidence$`RdRp(yes/no)` == 'yes')
                                       non_rdrp <- setdiff(1:length(contigs_name),rdrp)
（3）confidence1 <- confidence1[confidence1$Frequency < 2,]中的变为3（阈值设定）


