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
        if [ ! -d ${RESULT}/QC/${1} ]; then
                mkdir -p ${RESULT}/QC/${1}
        fi
        # if [ ! -d ${RESULT}/VCF/${1} ]; then
        #         mkdir -p ${RESULT}/QUANT/${1}
        # fi
}

DOWNLOAD_fn() {
    # It will take so long time!
    # if [ `grep -c "${1}" "${FQ_UPLOADED}"` -gt '0' ]; then
    #     aws s3 cp  --quiet --recursive s3://jiguang2021/DNA_RESULT/${1}/ /local_data/WGS/${1} \
    #         --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com 
    # fi
    if [ `grep -c "${1}" "${RESQ}"` -ne '0' ]; then # whether resq?

        echo "sample ${i} is resequenced!"

        aws s3 cp  --quiet --recursive s3://jiguang2021/ReSequence/${1} /local_data/WGS/${1} \
            --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com 
    fi
    aws s3 cp  --quiet --recursive s3://jiguang2021/WGS/${1} /local_data/WGS/${1} \
        --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com # ivolces for 内网；volces for 外网
    
}

QC_fn(){
    # if [ ! -d ALL_QC/ ]; then # 不存在则
    #     mkdir ALL_QC
    # fi
    # if [ `grep -c "${1}" "${FQ_UPLOADED}"` -eq '0' ]; then
    # 直接输入到目标文件
    fastp -i ${1}_1.fq.gz -I ${1}_2.fq.gz \
        -o out_${1}_1.fq.gz -O out_${1}_2.fq.gz \
        -q 25 -u 10 -l 50 -y -x -w $N -j ${RESULT}/QC/${1}/${1}_fastp.json -h ${RESULT}/QC/${1}/${1}_fastp.html -R "${1} Fastp Report"

    # 文件太大，应该删除原始文件
    rm ${1}_1.fq.gz ${1}_2.fq.gz & 
    # 上传处理后的，可能要在文件比对之前
    aws s3 cp . s3://jiguang2021/DNA_RESULT/${1}/ --recursive --quiet --exclude "*" \
        --include "out_${1}_[12].fq.gz"  --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com && (echo "SUCCESS upload out_${i}_[12]fq.gz ")
}

PRE_fn(){
    # if [ ! -d ALL_QC/ ]; then # 不存在则
    #     mkdir ALL_QC
    # fi
    #  除了超大的score.txt都在本地
    (${SENTIEON} bwa mem -M -R "@RG\tID:${1}\tSM:${1}\tPL:BGI\tCN:BGI" \
    -t ${N} -K 10000000 $REF out_${1}_1.fq.gz out_${1}_1.fq.gz || echo -n 'error' ) |\
     ${SENTIEON} util sort -r ${REF} -o ${1}_sorted.bam -t ${N} --sam2bam -i -
    # 
    ${SENTIEON} driver -r ${REF} -t ${N} -i ${1}_sorted.bam \
        --algo MeanQualityByCycle ${RESULT}/QC/${1}/mq_metrics.txt \
        --algo QualDistribution ${RESULT}/QC/${1}/qd_metrics.txt \
        --algo GCBias --summary ${RESULT}/QC/${1}/gc_summary.txt ${RESULT}/QC/${1}/gc_metrics.txt \
        --algo AlignmentStat --adapter_seq '' ${RESULT}/QC/${1}/aln_metrics.txt \
        --algo InsertSizeMetricAlgo ${RESULT}/QC/${1}/is_metrics.txt

    ${SENTIEON} plot GCBias -o ${RESULT}/QC/${1}/gc-report.pdf ${RESULT}/QC/${1}/gc_metrics.txt
    ${SENTIEON} plot QualDistribution -o ${RESULT}/QC/${1}/qd-report.pdf ${RESULT}/QC/${1}/qd_metrics.txt # 
    ${SENTIEON} plot MeanQualityByCycle -o ${RESULT}/QC/${1}/mq-report.pdf ${RESULT}/QC/${1}/mq_metrics.txt  # 
    ${SENTIEON} plot InsertSizeMetricAlgo -o ${RESULT}/QC/${1}/is-report.pdf ${RESULT}/QC/${1}/is_metrics.txt 

    ${SENTIEON} driver -t ${N} -i ${1}_sorted.bam --algo LocusCollector --fun score_info score.txt

    ${SENTIEON} driver -t ${N} -i ${1}_sorted.bam --algo Dedup --rmdup --score_info score.txt \
        --metrics ${RESULT}/QC/${1}/dedup_metrics.txt ${1}_deduped.bam 

    ${SENTIEON} driver -r $REF -t ${N} -i ${1}_deduped.bam --algo QualCal \
        -k ${DBSNP} -k ${INDEL_1000} -k ${INDEL_MILLS} ${RESULT}/QC/${1}/recal_data.table

    ${SENTIEON} driver -r $REF -t ${N} -i ${1}_deduped.bam \
        -q ${RESULT}/QC/${1}/recal_data.table --algo QualCal -k ${DBSNP} -k ${INDEL_1000} -k ${INDEL_MILLS} ${RESULT}/QC/${1}/recal_data.table.post

    ${SENTIEON} driver -t ${N} --algo QualCal \
        --plot --before ${RESULT}/QC/${1}/recal_data.table --after ${RESULT}/QC/${1}/recal_data.table.post ${RESULT}/QC/${1}/recal.csv   

    ${SENTIEON} plot QualCal -o ${RESULT}/QC/${1}/recal_plots.pdf ${RESULT}/QC/${1}/recal.csv
    ${SENTIEON} driver -r ${REF} -t ${N} -i ${1}_deduped.bam -q ${RESULT}/QC/${1}/recal_data.table --algo ReadWriter ${1}_recaled.bam

    # do not forget to add index file!
    aws s3 cp . s3://jiguang2021/DNA_RESULT/${1}/ --recursive --quiet --exclude "*" \
            --include "*_recaled.bam*" --include "*_deduped.bam*" --include "score.txt*" \
            --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com  && (echo "SUCCESS upload BAMs")
}
GERMLINE_fn(){
    if [ ! -d ${RESULT}/VCF/${1} ]; then
        mkdir -p ${RESULT}/VCF/${1}
    fi
    ${SENTIEON} -r ${REF} -t ${REF} -i ${1}_deduped.bam -q ${RESULT}/QC/${1}/recal_data.table --algo Haplotyper \
        -d ${DBSNP} --emit_conf=30 --call_conf=30 ${RESULT}/VCF/${1}/${1}_N.vcf.gz
}


########################################
#
# MAIN LOOP
#
########################################

for i in $(cat ${ALLOC});do


	if [ `grep -c "${i}" "${FINISHED}"` -ne '0' ]; then # 总体判断
		echo "This sample ${i} has been processed!"
        continue
    else
        # 必须有recursive
        start_time=`date +%F_%R`
        # start_time_s=`date +%s`
        SECONDS=0
        echo "start processing ${i} at ${start_time}"
        DIR_fun ${i}

        if [ ! -e ${i}/${i}_3.fq.gz ]; then
            DOWNLOAD_fn ${i} && (echo "SUCCESS download FQs")
        fi

        # if [ `grep -c "${i}" "${RESQ}"` -ne '0' ]; then # whether resq?
        #     aws s3 cp  --quiet --recursive s3://jiguang2021/ReSequence/${i} /local_data/WGS/${i} \
        #         --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com 
        # fi
        # aws s3 cp  --quiet --recursive s3://jiguang2021/WGS/${i} /local_data/WGS/${i} \
        #     --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com # ivolces for 内网；volces for 外网

        if [ $? -eq 0 ]; then # 下载完成
            cd ${i}
            touch ${i}_3.fq.gz # 此表名下载完成！
            #ls *fq.gz # 是否在根目录
            if [ $(ls -1d */|wc -l) -gt 0 ]; then # 有子目录
                mv */*.fq.gz .
            fi
            if [ $(ls -1 *_1.fq.gz|wc -l) -gt 1 ];then
                cat *_1.fq.gz > ${i}_1.fq.gz &
                cat *_2.fq.gz > ${i}_2.fq.gz & 
            fi

        fi
        waite 
        find . -type f -and -not -name "${i}_?.fq.gz" -exec rm -f {} \;
        #rm -rf !(${i}_1.fq.gz|${i}_2.fq.gz|DOWNLOAD_OK_${i}) & # 删除lane文件夹



        QC_fn ${i}
        PRE_fn ${i}
        # DIR_fun ${i}
        # aws s3 cp . s3://jiguang2021/DNA_RESULT/${i}/ --quiet --recursive --exclude "*" \
        #     --include "*_recaled.bam" --include "out_*_[12].fq.gz" --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com

        # from WT workflow
        # fastqc -o . -t 88 -q ${i}_1.fq.gz ${i}_2.fq.gz &
        # trim_galore -j 88 -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --paired  ${i}_1.fq.gz ${i}_2.fq.gz -o .
        # fastqc -o . -t 88 -q ${i}_1_val_1.fq.gz ${i}_2_val_2.fq.gz &
        #fastp -i ${i}_1.fq.gz -I ${i}_2.fq.gz -o ${i}_1_out.fq.gz -O ${i}_2_out.fq.gz -w 16 -h ${RESULT}/QC/${i}.html
        # STAR_fun ${i}
        # #kallisto quant -i $INDEX -o ${RESULT}/QUANT/${i} -b 100 ${i}_1_out.fq.gz ${i}_2_out.fq.gz -t 88 

        # if [ $? -eq 0 ]; then
        #         if [ ! -d RSEM ]; then  # RSEM 失败
        #                 mkdir RSEM
        #         fi
        #         rsem-calculate-expression --paired-end --no-bam-output --alignments -p 88 -q \
        #         ${i}_Aligned.toTranscriptome.out.bam ${RSEM_INDEX}/GRCh38 RSEM/${i}   # 观察输出文件
        # fi
        # DIR_fun ${i}
        # # if [ ! -d ${RESULT}/QC/${i} ]; then
        # #         mkdir -p ${RESULT}/QC/${i}
        # # fi
        # # if [ ! -d ${RESULT}/QUANT/${i} ]; then
        # #         mkdir -p ${RESULT}/QUANT/${i}
        # # fi

        # aws s3 cp ALL_QC s3://jiguang2021/DNA_RESULT/${i}/QC/ --quiet --recursive --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com
        # mv ALL_QC/* -t ${RESULT}/QC/${i}/
        # mv *.txt *.html *.zip *_Log* RSEM/${i}.stat -t ${RESULT}/QC/${i}/
        # mv *.tab  RSEM/*results -t ${RESULT}/QUANT/${i}/

        # aws s3 cp . s3://jiguang2021/RNA_RESULT/${i}/ --quiet --recursive --exclude "*" \
        # --include "*.bam" --include "*val_[12].fq.gz" --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com


        # 判断是否成功上传
        aws s3api head-object --bucket jiguang2021 --key DNA_RESULT/${i}/score.txt --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com

        if [ $? -eq 0 ]; then
                echo "${i}" >> ${FINISHED}
                echo "${i} 已被写入完成列表"  
        else 
            echo "SOMETHIGN FISHY!" 
            exit
        fi

  
        end_time=`date +%F_%H:%M:%S`
        echo "end processing ${i} at ${end_time}"
        duration_min=$(($SECONDS/60))
        echo "Using ${duration_min} seconds to QC and Prepare for ${i}"


        SECONDS=0
        #echo "start GERMLINE ${i} at ${start_time}"
        if [[ ${i} =~ "NR" ]]; then
            GERMLINE_fn ${i}
        else
            cd /local_data/WGS # 明确！
            rm -rf /local_data/WGS/${i} & # 可以放置到后台
            continue
        fi 

        if [ -f ${RESULT}/VCF/${1}/${1}_N.vcf.gz ]; then
            echo "${i}" >> ${GERMLINE}
            duration_min=$(($SECONDS/60))
            echo "Using ${duration_min} seconds to GERMLINE for ${i}"
            cd /local_data/WGS 
            rm -rf /local_data/WGS/${i} & # 可以放置到后台
        fi
	fi
done



########################################
#
# DEAD CODES
#
########################################
        #fastqc -t $N ${RNA_FQ}/${$1}_1.fq.gz ${RNA_FQ}/${$1}_2.fq.gz
        #fastqc -t $N ${RNA_FQ}/out_${$1}_1.fq.gz ${RNA_FQ}/out_${$1}_2.fq.gz
        



            # if [ $? -eq 0 ];then # 在根目录

            #         if [ $(ls -1 *_1.fq.gz|wc -l) -gt 1 ];then
            #                 cat *_1.fq.gz > ${i}_1.fq.gz
            #                 cat *_2.fq.gz > ${i}_2.fq.gz
            #         else
            #                 mv *_1.fq.gz ${i}_1.fq.gz
            #                 mv *_2.fq.gz ${i}_2.fq.gz
            #         fi

            #         # cd ..
            # else # 有可能多个文件夹！
            #         mv */*.fq.gz .
            #         #cd *
            #         if [ $(ls -1 *_1.fq.gz|wc -l) -gt 1 ];then
            #                 cat *_1.fq.gz > ${i}_1.fq.gz
            #                 cat *_2.fq.gz > ${i}_2.fq.gz
            #         else
            #                 mv *_1.fq.gz ${i}_1.fq.gz
            #                 mv *_2.fq.gz ${i}_2.fq.gz
            #         fi
            #         #cd ..
            # fi