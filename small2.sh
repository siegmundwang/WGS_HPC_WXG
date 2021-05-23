#!/usr/bin/env bash
## software: FastQC (version: 0.11.8), multiqc (version: 1.9)
## software: trim_galore (version: 0.6.5), cutadapt (version: 1.18)
## software: STAR (version: 2.7.2b)
## software: RSEM (version: 1.3.1)
# conda install -c bioconda fastqc=0.11.8 trim-galore=0.6.5 rsem=1.3.1 cutadapt=1.18 star=2.7.2b # 注意conda没有下划线
# 
# set -x
# shopt -s extglob
touch /data/Manager/RNA-Seq_Manage/Finished_STAR
ALLOC=/data/Manager/RNA-Seq_Manage/RNA$(printf "%02d" $(hostname|cut -f2 -d"-")) # 自动获取对应的文件安排
FINISHED=/data/Manager/RNA-Seq_Manage/Finished_STAR # 记录分析完成的样本名
BAD=/data/Manager/RNA-Seq_Manage/BAD
STAR_INDEX=/data/common_data/reference/STAR_INDEX_GRCh38
RSEM_INDEX=/data/common_data/reference/RSEM_INDEX_GRCh38
RESULT=/data/Results/RNA-Seq_STAR
#INDEX=/data/common_data/kallisto.idx


for i in $(cat ${ALLOC} | head -n 10)
do
        if [ -f ${RESULT}/QC/${i}/${i}.stat/${i}.cnt ] && [ -f ${RESULT}/QUANT/${i}/${i}.gene.results ]; then
                echo "This sample has been processed! And have proper quant file & QC file!"
                echo "${i}" >> ${FINISHED}
                echo "${i} 已被写入完成列表"  
                continue
done