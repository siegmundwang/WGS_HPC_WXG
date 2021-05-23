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
touch /data/Manager/DNA-Seq_Manage/Finished_WGS_VCF
ALLOC=/data/Manager/DNA-Seq-Manage/TN$(printf "%02d" $(hostname|cut -f2 -d"-")) # 自动获取对应的文件安排
FINISHED=/data/Manager/DNA-Seq_Manage/Finished_WGS_VCF # 记录分析完成的样本名
RESULT=/data/Results/DNA-Seq
#RESQ=/data/Manager/DNA-Seq-Manage/RESQ
#TN=/data/Manager/DNA-Seq_Manage/GOOD_TN

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

########################################
#
# MODULES
#
########################################
DIR_fun(){
    # 暂时用不到
        if [ ! -d ${RESULT}/VCF/${1} ]; then
                mkdir -p ${RESULT}/VCF/${1}
        fi
        # if [ ! -d ${RESULT}/VCF/${1} ]; then
        #         mkdir -p ${RESULT}/QUANT/${1}
        # fi
}


DOWNLOAD_fn() {
    # ${1} tumor; ${2} normal
    # # It will take so long time!
    # if [ `grep -c "${1}" "${RESQ}"` -ne '0' ]; then # whether resq?

    #     echo "sample ${i} is resequenced!"

    #     aws s3 cp  --quiet --recursive s3://jiguang2021/ReSequence/${1} /local_data/WGS/${1} \
    #         --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com 
    # fi
    # aws s3 cp  --quiet --recursive s3://jiguang2021/WGS/${1} /local_data/WGS/${1} \
    #     --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com # ivolces for 内网；volces for 外网

    aws s3 cp -quiet --recursive s3://jiguang2021/DNA_RESULT/${1}/ --quiet --exclude "*" \
        --include "*_deduped.bam" --include "*_recal_data.table" --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com
    aws s3 cp -quiet --recursive s3://jiguang2021/DNA_RESULT/${2}/ --quiet --exclude "*" \
        --include "*_deduped.bam" --include "*_recal_data.table" --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com
}

TN_fn(){
    # ${1} tumor; ${2} normal; ${3} patient

    ${SENTIEON} -t ${N} -r ${REF} \
   -i ${1}_deduped.bam -q ${1}_recal_data.table \
   -i ${2}_deduped.bam -q ${2}_recal_data.table \
   --algo TNhaplotyper2 --tumor_sample ${1} \
      --normal_sample ${2} \
    #   [--germline_vcf GERMLINE_RESOURCE] \
    #   [--pon PANEL_OF_NORMAL] \
      ${RESULT}/VCF/${3}/${3}_TN.vcf \
#    [ --algo OrientationBias --tumor_sample TUMOR_SAMPLE_NAME \
#       ORIENTATION_DATA ] \
#    [ --algo ContaminationModel --tumor_sample TUMOR_SAMPLE_NAME \
#       --normal_sample NORMAL_SAMPLE_NAME \
#       --vcf GERMLINE_RESOURCE \
#       --tumor_segments CONTAMINATION_DATA.segments \
#       CONTAMINATION_DATA ]

${SENTIEON} -r ${REF} \
   --algo TNfilter --tumor_sample ${1} \
   --normal_sample ${2} \
   -v ${RESULT}/VCF/${3}/${3}_TN_RAW.vcf \
#    [--contamination CONTAMINATION_DATA] \
#    [--tumor_segments CONTAMINATION_DATA.segments] \
#    [--orientation_priors ORIENTATION_DATA] \
   ${RESULT}/VCF/${3}/${3}_TN_FILTERED.vcf
}






OLDIFS=$IFS
IFS=','
while read Patient Normal Tumor; do
    #
    #  for i in $(cat ${ALLOC}); do
    if [ `grep -c "${Patient}" "${FINISHED}"` -ne '0' ]; then
        echo "This TN pair has been processed!"
    else
        DIR_fun ${Patient}
        DOWNLOAD_fn  ${Tumor} ${Normal}
        TN_fn  ${Tumor} ${Normal} ${Patient}
done < $ALLOC
IFS=$OLDIFS

