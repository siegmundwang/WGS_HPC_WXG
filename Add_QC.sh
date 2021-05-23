#!/usr/bin/env bash

#set -x # show each running command
# @File    :   DNA-Seq.sh
# @Time    :   2021/05/13 16:02:09
# @Author  :   Wang Xiangeng 
# @Version :   1.0
# @Contact :   meinetbosekampf@gmail.com
# @License :   (C)Copyright 2017-2018, Liugroup-NLPR-CASIA
# should be run inside /local_data/WGS
# 一些心得体会：
#   删除一定谨慎明确

########################################
#  
# GLOBAL CONTROL
# 
########################################

# shopt -s extglob
# touch /data/Manager/DNA-Seq_Manage/Finished_WGS_PREP
# touch /data/Manager/DNA-Seq_Manage/Finished_GERMLINE
#touch /data/Manager/DNA-Seq_Manage/FQ
#touch /data/Manager/DNA-Seq_Manage/BQSR
ALLOC=/data/Manager/DNA-Seq_Manage/WGS #$(printf "%02d" $(hostname|cut -f2 -d"-")) # 自动获取对应的文件安排
FINISHED=/data/Manager/DNA-Seq_Manage/Finished_WGS_PREP # 记录分析完成的样本名
GERMLINE=/data/Manager/DNA-Seq_Manage/Finished_GERMLINE
RESQ=/data/Manager/DNA-Seq_Manage/RESQ
RESULT=/data/Results/DNA-Seq
#FQ_UPLOADED=/data/Manager/DNA-Seq_Manage/FQ
#BQSRED=/data/Manager/DNA-Seq_Manage/BQSR

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

# nt=$( cat /proc/cpuinfo |grep processor|wc -l ) by hands, maybe better by command
N=$( nproc )

# for better compression
bam_option="--bam_compression 1"

Add_QC_fn(){

    ${SENTIEON} driver -r ${REF} -t ${N} -i ${1}_sorted.bam \
        --algo MeanQualityByCycle ${RESULT}/QC/${1}/mq_metrics.txt --algo QualDistribution ${RESULT}/QC/${1}/qd_metrics.txt \
        --algo GCBias --summary ${RESULT}/QC/${1}/gc_summary.txt ${RESULT}/QC/${1}/gc_metrics.txt \
        --algo AlignmentStat --adapter_seq '' ${RESULT}/QC/${1}/aln_metrics.txt --algo InsertSizeMetricAlgo ${RESULT}/QC/${1}/is_metrics.txt

    ${SENTIEON} plot GCBias -o ${RESULT}/QC/${1}/gc-report.pdf ${RESULT}/QC/${1}/gc_metrics.txt
    ${SENTIEON} plot QualDistribution -o ${RESULT}/QC/${1}/qd-report.pdf ${RESULT}/QC/${1}/qd_metrics.txt # 
    ${SENTIEON} plot MeanQualityByCycle -o ${RESULT}/QC/${1}/mq-report.pdf ${RESULT}/QC/${1}/mq_metrics.txt  # 
    ${SENTIEON} plot InsertSizeMetricAlgo -o ${RESULT}/QC/${1}/is-report.pdf ${RESULT}/QC/${1}/is_metrics.txt 
}


for i in $(cat ${ALLOC});do
    aws s3api head-object --bucket jiguang2021 --key DNA_RESULT/${i}/${i}_sorted.bam --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com
    if [ $? -eq 0 ]; then
        echo "${i} have good ${i}_sorted.bam"
    fi
done
