#!/usr/bin/env bash

#set -x # show each running command
# @File    :   DNA-Seq.sh
# @Time    :   2021/05/13 16:02:09
# @Author  :   Wang Xiangeng 
# @Version :   1.0
# @Contact :   meinetbosekampf@gmail.com
# @License :   (C)Copyright 2017-2018, Liugroup-NLPR-CASIA
# should be run inside /local_data/WGS


########################################
#  
# GLOBAL CONTROL
# 
########################################

shopt -s extglob
touch /data/Manager/DNA-Seq_Manage/Finished_WGS_PREP
touch /data/Manager/DNA-Seq_Manage/FQ
touch /data/Manager/DNA-Seq_Manage/BQSR
#ALLOC=/data/Manager/DNA-Seq-Manage/WGS$(printf "%02d" $(hostname|cut -f2 -d"-")) # 自动获取对应的文件安排
FINISHED=/data/Manager/DNA-Seq_Manage/Finished_WGS_PREP # 记录分析完成的样本名
RESQ=/data/Manager/DNA-Seq_Manage/RESQ
FQ_UPLOADED=/data/Manager/DNA-Seq_Manage/FQ
BQSRED=/data/Manager/DNA-Seq_Manage/BQSR



########################################
#
# WGS CONTROL
#
########################################


# the path for references 需要根据项目修改
REF_DIR=/data/common_data/gatk_bundle
REF=${REF_DIR}/Homo_sapiens_assembly38.fasta
#dbsnp=${REF_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf
DBSNP=${REF_DIR}/dbsnp_151.hg38.vcf.gz
INDEL_1000=${REF_DIR}/1000G_phase1.snps.high_confidence.hg38.vcf.gz
INDEL_MILLS=${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

# sentieon server and installation dir
export SENTIEON_LICENSE=192.168.10.112:8990
SENTIEON_INSTALL_DIR=/data/biosoft/sentieon-genomics-202010.01
SENTIEON=${SENTIEON_INSTALL_DIR}/bin/sentieon

# cores
# nt=$( cat /proc/cpuinfo |grep processor|wc -l ) by hands, maybe better by command
N=$( nproc )

# for better compression
bam_option="--bam_compression 1"
  
# 14TS0H0110TR1

${SENTIEON} driver -r ${REF} -t ${N} -i ${1}_sorted.bam \
    --algo MeanQualityByCycle ALL_QC/mq_metrics.txt --algo QualDistribution ALL_QC/qd_metrics.txt \
    --algo GCBias --summary ALL_QC/gc_summary.txt ALL_QC/gc_metrics.txt 
    --algo AlignmentStat --adapter_seq '' ALL_QC/aln_metrics.txt --algo InsertSizeMetricAlgo ALL_QC/is_metrics.txt

${SENTIEON} plot GCBias -o ALL_QC/gc-report.pdf ALL_QC/gc_metrics.txt
${SENTIEON} plot QualDistribution -o ALL_QC/qd-report.pdf ALL_QC/qd_metrics.txt
${SENTIEON} plot MeanQualityByCycle -o ALL_QC/mq-report.pdf ALL_QC/mq_metrics.txt
${SENTIEON} plot InsertSizeMetricAlgo -o ALL_QC/is-report.pdf ALL_QC/is_metrics.txt

${SENTIEON} driver -t ${N} -i ${1}_sorted.bam --algo LocusCollector --fun score_info ALL_QC/score.txt

${SENTIEON} driver -t ${N} -i ${1}_sorted.bam --algo Dedup --rmdup --score_info ALL_QC/score.txt \
    --metrics ALL_QC/dedup_metrics.txt ${1}_deduped.bam 

${SENTIEON} driver -r $REF -t ${N} -i ${1}_deduped.bam --algo QualCal \
    -k ${DBSNP} -k ${INDEL_1000} -k ${INDEL_MILLS} ALL_QC/recal_data.table

${SENTIEON} driver -r $REF -t ${N} -i ${1}_deduped.bam \
    -q ALL_QC/recal_data.table --algo QualCal -k ${DBSNP} -k ${INDEL_1000} -k ${INDEL_MILLS} ALL_QC/recal_data.table.post

${SENTIEON} driver -t ${N} --algo QualCal \
    --plot --before ALL_QC/recal_data.table --after ALL_QC/recal_data.table.post ALL_QC/recal.csv   

${SENTIEON} plot QualCal -o ALL_QC/recal_plots.pdf ALL_QC/recal.csv
${SENTIEON} driver -r ${REF} -t ${N} -i ${1}_deduped.bam -q ALL_QC/${1}_recal_data.table --algo ReadWriter ${1}_recaled.bam

aws s3 cp . s3://jiguang2021/DNA_RESULT/${1}/ --quiet --recursive --exclude "*" \
        --include "*_recaled.bam" --include "*_deduped.bam" --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com && (echo "${i}" ${BQSRED})


aws s3 cp ALL_QC s3://jiguang2021/DNA_RESULT/14TS0H0110TR1/QC/ --quiet --recursive --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com