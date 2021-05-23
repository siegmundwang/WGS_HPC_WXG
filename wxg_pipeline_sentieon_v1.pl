#!/usr/bin/perl -w 
use strict;
use Getop::Long;
my $usage=<<'USAGE';

Cancer pipeline of hg38 for BGISEQ-500 data and Sentieon Platform.
Version 1.0: 2019-9-30
Author : lifuqiang && zhengzetian, CANCER INSTITUTE
    
Options:
    -list            STR  reads info list.  <require>
    -config          STR  configuration of parameters to run the pipeline. <require>
    -outdir          STR  the output directory. <require>
    -rmCleanReads    STR  whether remove clean reads after finish Sentieon mapping. default is true. [optional: true|false]
    -mergeSI         STR  whether merge mutli-tools somatic snv/indel results. default is true. [optional: true|false]
    -oneSentieon     STR  whether to run some sentieon steps (mkdup, indel-realignment, BQSR, QC and SNP/INDEL calling) in one shell. default is false. [optional: true|false]
    -bam2cram        STR  whether convert bam to cram (could reduce file size at least 35%), and then remove bam. default is true. [optional: true|false]
    -bamList         STR  existed bam/cram list. cram will be converted to bam. [optional:bam/cram file list]
    -h                    Help Information
		 
!!! This pipeline is wrote base on chromosome name with "chr", be careful when replace default files!!!!!
USAGE

my($List,$config,$outdir,$rmCleanReads,$mergeSI,$oneSentieon,$bam2cram,$bamList);