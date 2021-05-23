#!/usr/bin/env bash
## software: FastQC (version: 0.11.8), multiqc (version: 1.9)
## software: trim_galore (version: 0.6.5), cutadapt (version: 1.18)
## software: STAR (version: 2.7.2b)
## software: RSEM (version: 1.3.1)
# conda install -c bioconda fastqc=0.11.8 trim-galore=0.6.5 rsem=1.3.1 cutadapt=1.18 star=2.7.2b # 注意conda没有下划线
# 
# set -x
shopt -s extglob
touch /data/Manager/RNA-Seq_Manage/Finished_STAR
ALLOC=/data/Manager/RNA-Seq_Manage/RNA$(printf "%02d" $(hostname|cut -f2 -d"-")) # 自动获取对应的文件安排
FINISHED=/data/Manager/RNA-Seq_Manage/Finished_STAR # 记录分析完成的样本名
BAD=/data/Manager/RNA-Seq_Manage/BAD
RERUN=
#INDEX=/data/common_data/kallisto.idx
STAR_INDEX=/data/common_data/reference/STAR_INDEX_GRCh38
RSEM_INDEX=/data/common_data/reference/RSEM_INDEX_GRCh38
RESULT=/data/Results/RNA-Seq_STAR
STAR --runThreadN 88 --runMode genomeGenerate \
        --genomeFastaFiles ./GRCh38.p13.genome.fa \
        --sjdbGTFfile ./gencode.v33.annotation.gtf \
        --sjdbOverhang 99 --genomeDir ./STAR_INDEX_GRCh38/
rsem-prepare-reference --gtf gencode.v33.annotation.gtf \
       GRCh38.p13.genome.fa RSEM_INDEX_GRCh38/GRCh38


STAR_fun(){
        # 注意输出文件
        STAR --genomeDir ${STAR_INDEX} --runThreadN 88 \
        --readFilesCommand zcat --readFilesIn ${1}_1_val_1.fq.gz ${1}_2_val_2.fq.gz \
        --outSAMtype BAM Unsorted --outFileNamePrefix ${1}_ \
        --outFilterType BySJout --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 \
        --quantMode TranscriptomeSAM GeneCounts
}

DIR_fun(){
        if [ ! -d ${RESULT}/QC/${1} ]; then
                mkdir -p ${RESULT}/QC/${1}
        fi
        if [ ! -d ${RESULT}/QUANT/${1} ]; then
                mkdir -p ${RESULT}/QUANT/${1}
        fi
}

for i in $(cat ${ALLOC} | head -n 10)
do
        if [ -f ${RESULT}/QC/${i}/${i}.stat/${i}.cnt ] && [ -f ${RESULT}/QUANT/${i}/${i}.gene.results ]; then
                echo "This sample has been processed! And have proper quant file & QC file!"
                echo "${i}" >> ${FINISHED}
                echo "${i} 已被写入完成列表"  
                continue

        elif [ `grep -c "${i}" "${BAD}"` -ne '0' ]; then
                echo "This sample is a bad one!"
                continue
        else
                # 必须有recursive
                start_time=`date +%F_%R`
                # start_time_s=`date +%s`
                SECONDS=0
                echo "start processing ${i} at ${start_time}"

                aws s3 cp  --quiet --recursive s3://jiguang2021/RNA/${i} /local_data/${i} \
                --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com # ivolces for 内网；volces for 外网

                if [ $? -eq 0 ]; then
                cd ${i}
                ls *fq.gz
                if [ $? -eq 0 ];then
                        if [ $(ls -1 *_1.fq.gz|wc -l) -gt 1 ];then
                                cat *_1.fq.gz > ${i}_1.fq.gz
                                cat *_2.fq.gz > ${i}_2.fq.gz
                        else
                                mv *_1.fq.gz ${i}_1.fq.gz
                                mv *_2.fq.gz ${i}_2.fq.gz
                        fi
                        # cd ..
                else # 有可能多个文件夹！
                        mv */*.fq.gz .
                        #cd *
                        if [ $(ls -1 *_1.fq.gz|wc -l) -gt 1 ];then
                                cat *_1.fq.gz > ${i}_1.fq.gz
                                cat *_2.fq.gz > ${i}_2.fq.gz
                        else
                                mv *_1.fq.gz ${i}_1.fq.gz
                                mv *_2.fq.gz ${i}_2.fq.gz
                        fi
                        #cd ..
                fi
                fi
                rm -rf !(*.fq.gz) & # 删除lane文件夹

                # from WT workflow
                fastqc -o . -t 88 -q ${i}_1.fq.gz ${i}_2.fq.gz &
                trim_galore -j 88 -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --paired  ${i}_1.fq.gz ${i}_2.fq.gz -o .
                fastqc -o . -t 88 -q ${i}_1_val_1.fq.gz ${i}_2_val_2.fq.gz &
                #fastp -i ${i}_1.fq.gz -I ${i}_2.fq.gz -o ${i}_1_out.fq.gz -O ${i}_2_out.fq.gz -w 16 -h ${RESULT}/QC/${i}.html
                STAR_fun ${i}
                #kallisto quant -i $INDEX -o ${RESULT}/QUANT/${i} -b 100 ${i}_1_out.fq.gz ${i}_2_out.fq.gz -t 88 

                if [ $? -eq 0 ]; then
                        if [ ! -d RSEM ]; then  # RSEM 失败
                                mkdir RSEM
                        fi
                        rsem-calculate-expression --paired-end --no-bam-output --alignments -p 88 -q \
                        ${i}_Aligned.toTranscriptome.out.bam ${RSEM_INDEX}/GRCh38 RSEM/${i}   # 观察输出文件
                fi
                DIR_fun ${i}
                # if [ ! -d ${RESULT}/QC/${i} ]; then
                #         mkdir -p ${RESULT}/QC/${i}
                # fi
                # if [ ! -d ${RESULT}/QUANT/${i} ]; then
                #         mkdir -p ${RESULT}/QUANT/${i}
                # fi

                mv *.txt *.html *.zip *_Log* RSEM/${i}.stat -t ${RESULT}/QC/${i}/
                mv *.tab  RSEM/*results -t ${RESULT}/QUANT/${i}/

                aws s3 cp . s3://jiguang2021/RNA_RESULT/${i}/ --quiet --recursive --exclude "*" \
                --include "*.bam" --include "*val_[12].fq.gz" --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com


                if [ -f ${RESULT}/QC/${i}/${i}.stat/${i}.cnt ]; then
                        echo "${i}" >> ${FINISHED}
                        echo "${i} 已被写入完成列表"    
                fi

                cd .. 
                end_time=`date +%F_%H:%M:%S`
                echo "end processing ${i} at ${end_time}"
                duration_min=$(($SECONDS/60))
                echo "Using ${duration_min} seconds to QC and quantification ${i}"
                rm -rf ${i} & # 可以放置到后台

	fi
done