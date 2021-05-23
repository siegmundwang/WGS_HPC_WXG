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
touch /data/Manager/DNA-Seq_Manage/Finished_WGS_PREP
touch /data/Manager/DNA-Seq_Manage/Finished_GERMLINE
#touch /data/Manager/DNA-Seq_Manage/FQ
#touch /data/Manager/DNA-Seq_Manage/BQSR
ALLOC=/data/Manager/DNA-Seq_Manage/WGS$(printf "%02d" $(hostname|cut -f2 -d"-")) # 自动获取对应的文件安排
FINISHED=/data/Manager/DNA-Seq_Manage/Finished_WGS_PREP # 记录分析完成的样本名
GERMLINE=/data/Manager/DNA-Seq_Manage/Finished_GERMLINE
RESQ=/data/Manager/DNA-Seq_Manage/RESQ
RESULT=/data/Results/DNA-Seq
#FQ_UPLOADED=/data/Manager/DNA-Seq_Manage/FQ
#BQSRED=/data/Manager/DNA-Seq_Manage/BQSR

i=14TS0H0125NR1



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

GERMLINE_fn(){
    if [ ! -d ${RESULT}/VCF/${1} ]; then
        mkdir -p ${RESULT}/VCF/${1}
    fi
    ${SENTIEON} -r ${REF} -t ${REF} -i ${1}_deduped.bam -q ${RESULT}/QC/${1}/recal_data.table --algo Haplotyper \
        -d ${DBSNP} --emit_conf=30 --call_conf=30 ${RESULT}/VCF/${1}/${1}_N.vcf.gz
}


    SECONDS=0
    if [[ ${i} =~ "NR" ]]; then
            GERMLINE_fn ${i}
    else
        cd /local_data/WGS # 明确！
        rm -rf /local_data/WGS/${i} & # 可以放置到后台
        #continue
    fi 

    if [ -f ${RESULT}/VCF/${1}/${1}_N.vcf.gz ]; then
        echo "${i}" >> ${GERMLINE}
        duration_min=$(($SECONDS/60))
        echo "Using ${duration_min} seconds to GERMLINE for ${i}"
        cd /local_data/WGS 
        rm -rf /local_data/WGS/${i} & # 可以放置到后台
    fi


    if [ -f ${RESULT}/VCF/${1}/${1}_N.vcf.gz ]; then
        echo "${i}" >> ${GERMLINE}
        duration_min=$(($SECONDS/60))
        echo "Using ${duration_min} seconds to GERMLINE for ${i}"
        cd /local_data/WGS 
        rm -rf /local_data/WGS/${i} & # 可以放置到后台
    fi