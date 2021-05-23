https://support.sentieon.com/quick_start/

https://support.sentieon.com/manual/

https://support.sentieon.com/appnotes/

BQSR的文件：/mnt/minio/node75/wt/5.GATK/BQSR

/mnt/minio/node75/wt/Readme



# https://support.sentieon.com/manual/DNAseq_usage/dnaseq/#dnaseq-step-usage
#! /bin/bash

NUMBER_THREADS=48
NUMBER_THREADS_MEM=48
NUMBER_THREADS_SORT=48
REFERENCE="/mnt/minio/node75/wt/genome/Human/Gencode/v33.assembly/GRCh38.primary_assembly.genome.fa"
PLATFORM="Illumina"
KNOWN_SITES="/mnt/minio/node75/wt/genome/Human/GATK/snp/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
DBSNP="/mnt/minio/node75/wt/genome/Human/GATK/snp/dbsnp_146.hg38.vcf.gz"
BWA_IND="/mnt/minio/node75/wt/genome/Human/Gencode/v33.assembly/BWA2IndexForWGS/GRCh38.primary_assembly.genome"



SAMPLE_NAME_N=$1
SAMPLE_NAME_T=$2
mkdir fq2VCF/$SAMPLE_NAME_T
SORTED_BAM_N="fq2VCF/$SAMPLE_NAME_T/$SAMPLE_NAME_N.sort.bam"
SORTED_BAM_T="fq2VCF/$SAMPLE_NAME_T/$SAMPLE_NAME_T.sort.bam"
DEDUPED_BAM_N="fq2VCF/$SAMPLE_NAME_T/$SAMPLE_NAME_N.sort.deduped.bam"
DEDUPED_BAM_T="fq2VCF/$SAMPLE_NAME_T/$SAMPLE_NAME_T.sort.deduped.bam"
REALIGNED_BAM_N="fq2VCF/$SAMPLE_NAME_T/$SAMPLE_NAME_N.sort.deduped.realgin.bam"
REALIGNED_BAM_T="fq2VCF/$SAMPLE_NAME_T/$SAMPLE_NAME_T.sort.deduped.realgin.bam"
OUT_TN_VCF="fq2VCF/$SAMPLE_NAME_T/$SAMPLE_NAME_T.raw.vcf.gz"

# Normal
#####################################################
# Map reads to reference
time (sentieon bwa mem -M -R "@RG\\tID:$SAMPLE_NAME_N\\tSM:$SAMPLE_NAME_N\\tPL:$PLATFORM" \
  -t $NUMBER_THREADS_MEM $BWA_IND rawfq/${SAMPLE_NAME_N}.R1.fastq.gz rawfq/${SAMPLE_NAME_N}.R2.fastq.gz || echo -n 'error' ) \
  | sentieon util sort -r $REFERENCE -o $SORTED_BAM_N -t $NUMBER_THREADS_SORT --sam2bam -i -
echo "Time ----- Map and sort" 

# Calculate data metric
time sentieon driver -t $NUMBER_THREADS -r $REFERENCE -i $SORTED_BAM_N \
  --algo GCBias --summary ${SORTED_BAM_N}.GC_SUMMARY_TXT ${SORTED_BAM_N}.GC_METRIC_TXT \
  --algo MeanQualityByCycle ${SORTED_BAM_N}.MQ_METRIC_TXT \
  --algo QualDistribution ${SORTED_BAM_N}.QD_METRIC_TXT \
  --algo InsertSizeMetricAlgo ${SORTED_BAM_N}.IS_METRIC_TXT  \
  --algo AlignmentStat ${SORTED_BAM_N}.ALN_METRIC_TXT
echo "Time ----- Calculate data metric"

time sentieon plot GCBias -o ${SORTED_BAM_N}.GC_METRIC_PDF ${SORTED_BAM_N}.GC_METRIC_TXT
time sentieon plot MeanQualityByCycle -o ${SORTED_BAM_N}.MQ_METRIC_PDF ${SORTED_BAM_N}.MQ_METRIC_TXT
time sentieon plot QualDistribution -o ${SORTED_BAM_N}.QD_METRIC_PDF ${SORTED_BAM_N}.QD_METRIC_TXT
time sentieon plot InsertSizeMetricAlgo -o ${SORTED_BAM_N}.IS_METRIC_PDF ${SORTED_BAM_N}.IS_METRIC_TXT

# Remove or mark duplicates
time sentieon driver -t $NUMBER_THREADS -i ${SORTED_BAM_N} \
  --algo LocusCollector --fun score_info ${SORTED_BAM_N}.SCORE.gz
echo "Time ----- Calculated duplicates matrix" 
time sentieon driver -t $NUMBER_THREADS -i ${SORTED_BAM_N} \
  --algo Dedup --score_info ${SORTED_BAM_N}.SCORE.gz  \
  --metrics ${SORTED_BAM_N}.DEDUP_METRIC_TXT $DEDUPED_BAM_N
echo "Time ----- remove duplicate"

# Indel realignment (optional)
time sentieon driver -t $NUMBER_THREADS -r $REFERENCE \
  -i $DEDUPED_BAM_N --algo Realigner -k $KNOWN_SITES $REALIGNED_BAM_N
echo "Time ----- Indel realgin"

# Base quality score recalibration (BQSR)
time sentieon driver -t $NUMBER_THREADS -r $REFERENCE \
  -i $REALIGNED_BAM_N --algo QualCal -k $KNOWN_SITES ${REALIGNED_BAM_N}.RECAL_DATA.TABLE
echo "Time ----- BQSR"
##################################################
  
# Tumor
##################################################
# Map reads to reference
time (sentieon bwa mem -M -R "@RG\\tID:$SAMPLE_NAME_T\\tSM:$SAMPLE_NAME_T\\tPL:$PLATFORM" \
  -t $NUMBER_THREADS_MEM $BWA_IND rawfq/${SAMPLE_NAME_T}.R1.fastq.gz rawfq/${SAMPLE_NAME_T}.R2.fastq.gz || echo -n 'error' ) \
  | sentieon util sort -r $REFERENCE -o $SORTED_BAM_T -t $NUMBER_THREADS_SORT --sam2bam -i -
echo "Time ----- Map and sort" 


# Calculate data metric
time sentieon driver -t $NUMBER_THREADS -r $REFERENCE -i $SORTED_BAM_T \
  --algo GCBias --summary ${SORTED_BAM_T}.GC_SUMMARY_TXT ${SORTED_BAM_T}.GC_METRIC_TXT \
  --algo MeanQualityByCycle ${SORTED_BAM_T}.MQ_METRIC_TXT \
  --algo QualDistribution ${SORTED_BAM_T}.QD_METRIC_TXT \
  --algo InsertSizeMetricAlgo ${SORTED_BAM_T}.IS_METRIC_TXT  \
  --algo AlignmentStat ${SORTED_BAM_T}.ALN_METRIC_TXT
echo "Time ----- Calculate data metric"

time sentieon plot GCBias -o ${SORTED_BAM_T}.GC_METRIC_PDF ${SORTED_BAM_T}.GC_METRIC_TXT
time sentieon plot MeanQualityByCycle -o ${SORTED_BAM_T}.MQ_METRIC_PDF ${SORTED_BAM_T}.MQ_METRIC_TXT
time sentieon plot QualDistribution -o ${SORTED_BAM_T}.QD_METRIC_PDF ${SORTED_BAM_T}.QD_METRIC_TXT
time sentieon plot InsertSizeMetricAlgo -o ${SORTED_BAM_T}.IS_METRIC_PDF ${SORTED_BAM_T}.IS_METRIC_TXT

# Remove or mark duplicates
time sentieon driver -t $NUMBER_THREADS -i ${SORTED_BAM_T} \
  --algo LocusCollector --fun score_info ${SORTED_BAM_T}.SCORE.gz
echo "Time ----- Calculated duplicates matrix"
time sentieon driver -t $NUMBER_THREADS -i ${SORTED_BAM_T} \
  --algo Dedup --score_info ${SORTED_BAM_T}.SCORE.gz  \
  --metrics ${SORTED_BAM_T}.DEDUP_METRIC_TXT $DEDUPED_BAM_T
echo "Time ----- remove duplicate"


# Indel realignment (optional)
time sentieon driver -t $NUMBER_THREADS -r $REFERENCE \
  -i $DEDUPED_BAM_T --algo Realigner -k $KNOWN_SITES $REALIGNED_BAM_T
echo "Time ----- Indel realgin"

# Base quality score recalibration (BQSR)
time sentieon driver -t $NUMBER_THREADS -r $REFERENCE \
  -i $REALIGNED_BAM_T --algo QualCal -k $KNOWN_SITES ${REALIGNED_BAM_T}.RECAL_DATA.TABLE
echo "Time ----- BQSR"
##################################################  
  
# co-realgin and call
##################################################
# 参考下面的或者one command to perform both
time sentieon driver -t $NUMBER_THREADS -r $REFERENCE -i $REALIGNED_BAM_T \
  -q ${REALIGNED_BAM_T}.RECAL_DATA.TABLE -i $REALIGNED_BAM_N \
  -q ${REALIGNED_BAM_N}.RECAL_DATA.TABLE --algo TNhaplotyper \
  --tumor_sample $SAMPLE_NAME_T --normal_sample $SAMPLE_NAME_N \
  --dbsnp $DBSNP $OUT_TN_VCF
echo "Time ----- co-realgin and call"



