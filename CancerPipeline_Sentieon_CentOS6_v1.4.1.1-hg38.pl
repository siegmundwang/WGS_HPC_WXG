#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my $usage=<<'USAGE';
    
Cancer pipeline of hg38 for BGISEQ-500 data and Sentieon Platform.
Version 1.4: 2019-9-30
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

GetOptions(
    'list=s' => \$List,
    'config=s' => \$config,
    'outdir=s' => \$outdir,
    'rmCleanReads=s' => \$rmCleanReads,
    'mergeSI=s' => \$mergeSI,
    'oneSentieon=s' => \$oneSentieon,
    'bam2cram=s' => \$bam2cram,
    'bamList=s' => \$bamList,

    'h|?' => sub{die "$usage\n";},
);
die "$usage\n" unless (defined($List) && defined($config) && defined($outdir));
$rmCleanReads ||= "true";
die "String input for '-rmCleanReads' is typo: must be true or false\n" unless ($rmCleanReads=~/true/i || $rmCleanReads=~/false/i);
$mergeSI ||= "true";
die "String input for '-mergeSI' is typo: must be true or false\n" unless ($mergeSI=~/true/i || $mergeSI=~/false/i);
$oneSentieon ||= "false";
die "String input for '-oneSentieon' is typo: must be true or false\n" unless ($oneSentieon=~/true/i || $oneSentieon=~/false/i);
$bam2cram ||= "true";
die "String input for '-bam2cram' is typo: must be true or false\n" unless ($bam2cram=~/true/i || $bam2cram=~/false/i);

my(%cleanList, %samples, %sampleBam, %hash_sampair,%ALL_SAMP);

my $platform="COMPLETE_GENOMICS"; ## sequencing platform, fixed[COMPLETE_GENOMICS]
my($seqType); # library constructions type [target|wgs]
my($BIN,$DATAPATH,$reference,$Sentieonindexpath,$targetRegion,$queue,$priority,$subProjectName,$trname,$Sentieon,$SentieonAuto,$SentieonCPU); ## basic
my($SOAPnuke,$SOAPnukePara,$fqcheck);# fq clean 
my($samtools,$bamtools,$bcftools,$bedtools,$bgzip,$tabix,$bam_readcount,$bedtools_2_26_0); # bam & vcf tools
my($covdepPara); ## covdep
my($Dbsnp,$Dbsnpcommon,$kgIndel,$millIndel,$gnomad,$Cosmic,$ExACcommon); # extended database
my($somatic,$MuSE,$MuTect,$Platypus,$GATK2,$GATK3,$GATK4,$novoBreak,$factera,$blast,$bam_somaticsniper,$Manta,$Strelka2,$Lancet,$Svaba); # somatic calling tools 
my($Use_TNcaller,$Use_MuSE,$Use_MuTect,$Use_MuTect2,$Use_Platypus,$Use_FACETS,$Use_novoBreak,$Use_SomaticSniper,$Use_Lancet,$Use_Svaba,$Use_Manta,$Use_Strelka2,$Use_MantaAndStrelka2);# wihch somatic calling tool be used
my($MuSE_TRsplitNum,$MuTect_TRsplitNum,$MuTect2_TRsplitNum,$SomaticSniper_TRsplitNum,$Platypus_TRsplitNum,$Strelka2_TRsplitNum);
my($PlatypusPara,$PlatypussmtStrictPara,$PlatypusSomaticIndelFilter,$vcfutils,$Vt,$AnnotSV);# data and post processing for Platypus
my($pythonBin,$ldlib,$pythonPath); ## python support to Platypus,must in conformity with Platypus
my($GATKsomaticIndelPara,$GATKsmtFlitPara);# post processing for somaticindeldetector of gatk2
my($bammathPATH, $hapmapAF, $PanCanQC); # qc : bam-matcher & ContEst & PanCanQC
my($MuTectFilter); # post processing for MuTect
my($annovardb,$annovar); ## annotated by annovar
my($java7,$java8,$Picard); ## java7: GATK2, GATK3(<=3.5), MuTect(1.1.7), picard(<=1.141), VarScan2; java8: GATK4, GATK3(>=3.6), picard(>=2.0.1).
my($minSnvTools,$minIndelTools);

################################ Read configuration file ###################################

open CFG,"$config" or die $!;
while(<CFG>){
    chomp;
    next if (/^#/);
    next unless (/\w+/);

    if($_=~/^queue=(.*?);/){$queue=$1;}
    elsif($_ =~ /^priority=(.*?);/){$priority=$1;}
    elsif($_ =~ /^subProjectName=(.*?);/){$subProjectName=$1;}
    elsif($_ =~ /^platform=(.*?);/){$platform=$1;}
    elsif($_ =~ /^seqType=(.*?);/){$seqType=$1;}
    elsif($_ =~ /^targetRegion=(.*?);/){$targetRegion=$1;}
    elsif($_ =~ /^somatic=(.*?);/){$somatic=$1;}
    elsif($_ =~ /^reference=(.*?);/){$reference=$1;}
    elsif($_ =~ /^Sentieonindexpath=(.*?);/){$Sentieonindexpath=$1;}
    elsif($_ =~ /^DATAPATH=(.*?);/){$DATAPATH=$1;}
    elsif($_ =~ /^BIN=(.*?);/){$BIN=$1;}
    elsif($_ =~ /^Sentieon=(.*?);/){$Sentieon=$1;}
    elsif($_ =~ /^SentieonAuto=(.*?);/){$SentieonAuto=$1;}
    elsif($_ =~ /^SentieonCPU=(.*?);/){$SentieonCPU=$1;}
    elsif($_ =~ /^annovardb=(.*?);/){$annovardb=$1;}
    elsif($_ =~ /^annovar=(.*?);/){$annovar=$1;}
    elsif($_ =~ /^SOAPnuke=(.*?);/){$SOAPnuke=$1;}
    elsif($_ =~ /^SOAPnukePara=(.*?);/){$SOAPnukePara=$1;}
    elsif($_ =~ /^MuSE=(.*?);/){$MuSE=$1;}
    elsif($_ =~ /^MuTect=(.*?);/){$MuTect=$1;}
    elsif($_ =~ /^bam_somaticsniper=(.*?);/){$bam_somaticsniper=$1;}
    elsif($_ =~ /^Manta=(.*?);/){$Manta=$1;}
    elsif($_ =~ /^Strelka2=(.*?);/){$Strelka2=$1;}
    elsif($_ =~ /^Lancet=(.*?);/){$Lancet=$1;}
    elsif($_ =~ /^Svaba=(.*?);/){$Svaba=$1;}
    elsif($_ =~ /^GATK2=(.*?);/){$GATK2=$1;}
    elsif($_ =~ /^GATK3=(.*?);/){$GATK3=$1;}
    elsif($_ =~ /^Platypus=(.*?);/){$Platypus=$1;}
    elsif($_ =~ /^novoBreak=(.*?);/){$novoBreak=$1;}
    elsif($_ =~ /^Use_TNcaller=(.*?);/){$Use_TNcaller=$1;}
    elsif($_ =~ /^Use_MuSE=(.*?);/){$Use_MuSE=$1;}
    elsif($_ =~ /^Use_MuTect=(.*?);/){$Use_MuTect=$1;}
    elsif($_ =~ /^Use_MuTect2=(.*?);/){$Use_MuTect2=$1;}
    elsif($_ =~ /^Use_Manta=(.*?);/){$Use_Manta=$1;}
    elsif($_ =~ /^Use_Strelka2=(.*?);/){$Use_Strelka2=$1;}
    elsif($_ =~ /^Use_MantaAndStrelka2=(.*?);/){$Use_MantaAndStrelka2=$1;}
    elsif($_ =~ /^Use_SomaticSniper=(.*?);/){$Use_SomaticSniper=$1;}
    elsif($_ =~ /^Use_Lancet=(.*?);/){$Use_Lancet=$1;}
    elsif($_ =~ /^Use_Platypus=(.*?);/){$Use_Platypus=$1;}
    elsif($_ =~ /^Use_FACETS=(.*?);/){$Use_FACETS=$1;}
    elsif($_ =~ /^Use_novoBreak=(.*?);/){$Use_novoBreak=$1;}
    elsif($_ =~ /^Use_Svaba=(.*?);/){$Use_Svaba=$1;}
    elsif($_ =~ /^MuSE_TRsplitNum=(.*?);/){$MuSE_TRsplitNum=$1;}
    elsif($_ =~ /^MuTect_TRsplitNum=(.*?);/){$MuTect_TRsplitNum=$1;}
    elsif($_ =~ /^MuTect2_TRsplitNum=(.*?);/){$MuTect2_TRsplitNum=$1;}
    elsif($_ =~ /^SomaticSniper_TRsplitNum=(.*?);/){$SomaticSniper_TRsplitNum=$1;}
    elsif($_ =~ /^Platypus_TRsplitNum=(.*?);/){$Platypus_TRsplitNum=$1;}
    elsif($_ =~ /^Strelka2_TRsplitNum=(.*?);/){$Strelka2_TRsplitNum=$1;}
    elsif($_ =~ /^GATKsomaticIndelPara=(.*?);/){$GATKsomaticIndelPara=$1;}
    elsif($_ =~ /^GATKsmtFlitPara=(.*?);/){$GATKsmtFlitPara=$1;}
    elsif($_ =~ /^PlatypusPara=(.*?);/){$PlatypusPara=$1;}
    elsif($_ =~ /^PlatypussmtStrictPara=(.*?);/){$PlatypussmtStrictPara=$1;}
    elsif($_ =~ /^java7=(.*?);/){$java7=$1;}
    elsif($_ =~ /^java8=(.*?);/){$java8=$1;}
    elsif($_ =~ /^Dbsnp=(.*?);/){$Dbsnp=$1;}
    elsif($_ =~ /^kgIndel=(.*?);/){$kgIndel=$1;}
    elsif($_ =~ /^millIndel=(.*?);/){$millIndel=$1;}
    elsif($_ =~ /^gnomad=(.*?);/){$gnomad=$1;}
    elsif($_ =~ /^Dbsnpcommon=(.*?);/){$Dbsnpcommon=$1;}
    elsif($_ =~ /^ExACcommon=(.*?);/){$ExACcommon=$1;}
    elsif($_ =~ /^Cosmic=(.*?);/){$Cosmic=$1;}
    elsif($_ =~ /^minSnvTools=(.*?);/){$minSnvTools=$1;}
    elsif($_ =~ /^minIndelTools=(.*?);/){$minIndelTools=$1;}
}

################################### main parameters #####################################
$BIN ||= "/hwfssz1/ST_CANCER/POL/SHARE/CancerPipeline/Sentieon_V1/Bin_CentOS6";
$DATAPATH ||= "/hwfssz1/ST_CANCER/POL/SHARE/CancerPipeline/Sentieon_V1/DataBase_hg38";
$seqType ||="wgs"; 
die "The bin directory $BIN does not exist!\n" if(!-e $BIN);
die "The database directory $DATAPATH does not exist!\n" if(!-e $DATAPATH);

$Sentieonindexpath ||= "$DATAPATH/Sentieon_hg38_v201808.06";
die "The hash index $Sentieonindexpath does not exist!\n" if (!-e $Sentieonindexpath);
$reference ||= "$Sentieonindexpath/Homo_sapiens_assembly38.fasta";
die "The reference $reference does not exist!\n" if(!-e $reference);
my $reffai = "$reference.fai";
die "The reference index $reffai does not exist!\n" if(!-e $reffai);
my $refdict = "$reference.dict";
if(!-e $refdict){
    $refdict =~ s/\.[^\.]+\.dict$/.dict/;
    if(!-e $refdict){
        die "The reference dictionary $refdict or $reference.dict does not exist!\n";
    }
}
my $ref2bit = "$reference.2bit";
if(!-e $ref2bit){
    $ref2bit =~ s/\.[^\.]+\.2bit$/.2bit/;
    if(!-e $ref2bit){
        die "The 2bit reference $ref2bit or $reference.2bit does not exist!\n";
    }
}

$Dbsnp ||= "$DATAPATH/dbsnp_151.hg38.vcf.gz";
$kgIndel ||= "$DATAPATH/1000G_phase1.snps.high_confidence.hg38.vcf.gz";
$millIndel ||= "$DATAPATH/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz";
$gnomad ||= "$DATAPATH/af-only-gnomad.hg38.vcf.gz";
$Cosmic ||= "$DATAPATH/cosmic_v89.hg38.vcf.gz";
$Dbsnpcommon ||= "$DATAPATH/dbsnp_151.common.hg38.vcf.gz";
$ExACcommon ||= "$DATAPATH/small_exac_common_3.hg38.vcf.gz";
$targetRegion ||= "$DATAPATH/chrmosomes.bed";
$minSnvTools ||="2";
$minIndelTools ||="2";

##annovar
$annovar ||="$BIN/ANNOVAR/v20180416";
$annovardb ||="$BIN/ANNOVAR/humandb20180504";
die "$annovar does not exist!\n" if (!-e "$annovar");
die "$annovardb does not exist!\n" if (!-e "$annovardb");

my %softwares;
$Use_TNcaller ||= "false";
$Use_MuSE ||= "false"; $softwares{'MuSE'}=$Use_MuSE;
$Use_MuTect ||= "false";
$Use_MuTect2 ||= "false"; $softwares{'MuTect2'}=$Use_MuTect2;
$Use_Platypus ||= "false";
$Use_FACETS ||= "false";
$Use_novoBreak ||= "false";
$Use_Strelka2 ||= "false"; $softwares{'Strelka2'}=$Use_Strelka2;
$Use_SomaticSniper ||= "false"; $softwares{'SomaticSniper'}=$Use_SomaticSniper;
$Use_Lancet ||= "false"; $softwares{'Lancet'}=$Use_Lancet;
$Use_Manta ||= "false";
$Use_Svaba ||= "false"; $softwares{'Svaba'}=$Use_Svaba;
$Use_MantaAndStrelka2 ||="true";
my %SoftTRs;
$MuSE_TRsplitNum ||= 1; if($Use_MuSE=~/true/i){$SoftTRs{'MuSE'}="$MuSE_TRsplitNum";}
$MuTect_TRsplitNum ||= 1; if($Use_MuTect=~/true/i){$SoftTRs{'MuTect'}="$MuTect_TRsplitNum";}
$MuTect2_TRsplitNum ||= 1; if($Use_MuTect2=~/true/i){$SoftTRs{'MuTect2'}="$MuTect2_TRsplitNum";}
$SomaticSniper_TRsplitNum ||= 1; if($Use_SomaticSniper=~/true/i){$SoftTRs{'SomaticSniper'}="$SomaticSniper_TRsplitNum";}
$Platypus_TRsplitNum ||= 1; if($Use_Platypus=~/true/i){$SoftTRs{'Platypus'}="$Platypus_TRsplitNum";}
$Strelka2_TRsplitNum ||= 1; if($Use_Strelka2=~/true/i){$SoftTRs{'Strelka2'}="$Strelka2_TRsplitNum";}

## tools
$Sentieon ||= "$BIN/Sentieon";
$SentieonAuto ||= "$BIN/SentieonAuto";
$SentieonCPU ||= 48;
$SOAPnuke ||= "$BIN/SOAPnuke2.0.7";
$fqcheck ||= "$BIN/fqcheck33";
$samtools ||= "$BIN/samtools-1.5/bin/samtools";
$vcfutils ||= "$BIN/samtools-1.5/bin/vcfutils.pl";
$Vt ||= "$BIN/vt-0.57721";
$bamtools ||="$BIN/bamtools-2.5.1/bin/bamtools";
$bcftools ||="$BIN/bcftools-1.9";
$bedtools ||="$BIN/bedtools-2.27.1";
$bedtools_2_26_0 ||="$BIN/bedtools-2.26.0";
$bgzip ||= "$BIN/samtools-1.5/bin/bgzip";
$tabix ||= "$BIN/samtools-1.5/bin/tabix";
$java7 ||="/usr/java/latest/bin/java";
$java8 ||="/hwfssz1/ST_CANCER/POL/SHARE/tools/JAVA/jre1.8.0_92/bin/java";
$GATK2 ||="$BIN/GATK/GenomeAnalysisTK-2.3-9-ge5ebf34/GenomeAnalysisTK.jar";
$GATK3 ||="$BIN/GATK/GenomeAnalysisTK-3.8-1.jar";
$GATK4 ||="$BIN/GATK/gatk-4.1.0.0/gatk";
$Picard ||="$BIN/Picard/picard-2.18.27.jar";
$MuTect ||="$BIN/GATK/mutect-1.1.7.jar";
$MuSE ||="$BIN/MuSE";
$bam_readcount ||= "$BIN/bam-readcount-0.8.0";
$bam_somaticsniper ||= "$BIN/bam-somaticsniper";
$Manta ||= "$BIN/manta-1.5.0/bin/configManta.py";
$Strelka2 ||= "$BIN/strelka-2.9.9/bin/configureStrelkaSomaticWorkflow.py";
$Lancet ||= "$BIN/lancet";
$Svaba ||= "$BIN/svaba-0.2.1";
$AnnotSV ||="$BIN/AnnotSV/v2.0/bin/AnnotSV";

$Platypus ||="$BIN/Platypus-0.8.1/bin";
$novoBreak ||="$BIN/novoBreak-1.1.3";
$factera ||="$BIN/FACTERA-v1.4.4";
$blast ||="$BIN/blast+-v2.7.1";

$MuTectFilter="$BIN/mutectFilter.pl";
$PlatypusSomaticIndelFilter="$BIN/PlatypusSomaticIndelFilter.V2.pl";
$PanCanQC ||="$BIN/PanCanQC-1.2.2";

## python path
$pythonBin="/opt/python/bin";
$ldlib="$BIN/Platypus-0.8.1/htslib-1.3.2:$BIN/Platypus-0.8.1/xz-5.2.3/lib:/opt/python/lib";
$pythonPath="$BIN/Platypus-0.8.1/pythonLIB";

## clean, coverage
$SOAPnukePara ||= "-l 5 -q 0.5 -n 0.1 -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG";
$covdepPara ||= "--maxdepth 1000 -q 1";

## Platypus parameters
$PlatypusPara ||= "--filterDuplicates=1 --nCPU=4 --logFileName=/dev/null";
$PlatypussmtStrictPara ||= "--minPosterior 5";

## somatic QC
$bammathPATH="$BIN/BAM-matcher-v20170604";
$hapmapAF="$DATAPATH/hapmap_3.3_grch38_pop_stratified_af.vcf.gz";

my @r = ('a'..'z','A'..'Z');
$subProjectName ||= join '', map { $r[int rand @r] } 0..6;

################################ prepare Target Region #####################################
if(defined $targetRegion && -e $targetRegion){
    $trname = &parseTR($targetRegion);
  
}else{
    die "the target region file for the WES/Panel/WGS sequence $targetRegion must be specified or does not exist!\n";
}

############################### unify chromosome names ####################################

my @chrs=();
if (`grep chr $reffai`){
    for (my $i=1; $i<=22; $i++) {push @chrs,"chr$i";}
    push @chrs,"chrX";
    push @chrs,"chrY";
    push @chrs,"chrM";
}
else{
    for (my $i=1; $i<=22; $i++) {push @chrs,"$i";}
    push @chrs,"X";
    push @chrs,"Y";
    push @chrs,"MT";
}
my $chrlistcat = join " ", @chrs;
my $chrlistrm = $chrlistcat;
$chrlistcat =~ s/^1\s+//;
$chrlistcat =~ s/^chr1\s+//;

################################ Process sequence ##########################################

open LL,$List or die $!;
while(<LL>){
    chomp;
    my($q,$id);
    my ($samp,$qual,$readsId,$lib,$lane,$lane_barcode,$fqprefix)=split /\s+/, $_;
    if($qual==33){ $q=2; }else{ $q=1; }
    if($readsId =~ /old/i){ $id=0 ; }elsif($readsId =~ /new/i){ $id=1 ; }#else{ die "please offer reads id is old or new !";}
    `mkdir -p $outdir/$samp`;
	`mkdir -p $outdir/$samp/Sentieon`;
	`mkdir -p $outdir/$samp/shell`;
    `mkdir -p $outdir/$samp/$lane_barcode`;
    `mkdir -p $outdir/$samp/$lane_barcode/clean`;
    `mkdir -p $outdir/$samp/$lane_barcode/Sentieon`;
    `mkdir -p $outdir/$samp/$lane_barcode/shell`;
	
	$samples{$samp}{$lib}{$lane_barcode} = 1;
	
	my $rgid = $lane_barcode;
	$rgid=~s/_read//;
    my $fq1="$outdir/$samp/$lane_barcode/clean/$lane_barcode\_1.clean.fq.gz";
    my $fq2="$outdir/$samp/$lane_barcode/clean/$lane_barcode\_2.clean.fq.gz";
	
	&peLaneClean($samp,$fqprefix,$lane_barcode,$q,$id);
    &peLaneAln($samp,$rgid,$lib,$lane,$lane_barcode,$platform,$fq1,$fq2);
    my $lane_barcode_bam = "$outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.bam";

    push @{$cleanList{$samp}}, "$rgid,$samp,$lib,$platform,$fq1,$fq2,$lane_barcode_bam";
}
close LL;

for my $samp (sort keys %cleanList){
    my @bamarray;
	my $inbamcsv = "$outdir/$samp/Sentieon/01bamList4Germline.csv";
	open BAMCSV, ">$inbamcsv" || die $!;
    for my $csvline (@{$cleanList{$samp}}){
        print BAMCSV "$csvline\n";
        my @F = split /,/, $csvline;
        push @bamarray, $F[-1];
    }
	close BAMCSV;

	&RunGermline($samp, \@bamarray);
}

my %cram2bam;
if($bamList){
    open BAMLIST,"$bamList" or die $!;
    while(<BAMLIST>){
        chomp;
        my @arr=split/\s+/,$_;
        my($samp,$bam)=($arr[0],$arr[1]);
        my $bamRG=`$samtools view -H $bam | grep "^\@RG"|head -1|awk -F "SM:" '{print \$2}'|cut -f1`;
        chomp $bamRG;
        if($samp ne $bamRG){
            die "the sample name $samp you given for existed bam is different from it's RG $bamRG\n";
        }
        if(exists $samples{$samp}){
            die "sample $samp already given in $List, please don't given repeatedly\n";
        }
        $sampleBam{$samp}=$bam if($bam=~/\.bam$/);
        if($bam=~/\.bam$/ && ! -s "$bam.bai"){
        	die "the bam index (*.bam.bai) should be in your given bam path: $bam\n";
        }
        elsif($bam=~/\.cram$/ && ! -s "$bam.crai"){
            die "the cram index (*.cram.crai) should be in your given cram path: $bam\n";
        }
        else{
        	`mkdir -p $outdir/$samp`;
            `mkdir -p $outdir/$samp/shell`;
            if($bam=~/\.cram$/){
                my ($prefix) = $bam =~ /([^\/]+)\.cram/;
                $sampleBam{$samp}="$outdir/$samp/$prefix.bam" if($bam=~/\.cram$/);
                $cram2bam{$samp}{"S"} = "$outdir/$samp/shell/Sentieon.cram2bam.$samp.start.sh";
                open CVSH, ">$outdir/$samp/shell/Sentieon.cram2bam.$samp.start.sh" || die $!;
                print CVSH "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
                print CVSH "$SentieonAuto driver -t $SentieonCPU -r $reference -i $bam --algo ReadWriter $outdir/$samp/$prefix.bam --algo AlignmentStat $bam.ALN_METRIC.txt && \\\n";
                print CVSH "$SentieonAuto driver -t $SentieonCPU -r $reference -i $outdir/$samp/$prefix.bam --algo AlignmentStat $outdir/$samp/$prefix.bam.ALN_METRIC.txt && \\\n";
                print CVSH "LineNum=`grep -wFf $bam.ALN_METRIC.txt $outdir/$samp/$prefix.bam.ALN_METRIC.txt | wc -l` && \\\n";
                print CVSH "if [ \$LineNum -eq 6 ];then echo cram was converted to bam successfully;else echo cram was converted to bam unsuccessfully 1>&2;exit;fi && \\\n";
                print CVSH echostring("$outdir/$samp/shell/Sentieon.cram2bam.$samp.start.sh");
                close CVSH;

                $cram2bam{$samp}{"E"} = "$outdir/$samp/shell/Sentieon.cram2bam.$samp.end.sh";
                open CLSH, ">$outdir/$samp/shell/Sentieon.cram2bam.$samp.end.sh" || die $!;
                print CLSH "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
                print CLSH "rm -rf $outdir/$samp/$prefix.bam $outdir/$samp/$prefix.bam.bai $bam.ALN_METRIC.txt $outdir/$samp/$prefix.bam.ALN_METRIC.txt && \\\n";
                print CLSH echostring("$outdir/$samp/shell/Sentieon.cram2bam.$samp.end.sh");
                close CLSH;
            }
        }
    }
}

######################## somatic variant calling: Ssnv, Sindel #########################
my %NT;
if($somatic){
    open SMTF,"$somatic" or die $!;
    while(<SMTF>){
        chomp;
	      next unless (/\w+/);
        my @line=split /\s+/,$_;
        my $normal=$line[0];
        my $tumor=$line[1];

        $ALL_SAMP{$normal}=0;
        $ALL_SAMP{$tumor}=0;
        ### used already available bam result ###
        if(!$sampleBam{$normal}){
            my @MF = glob("$outdir/$normal/Sentieon/$normal*.bam");
            for my $tmp_bf (@MF){
                if ($tmp_bf =~ /\.bqsr\./){$sampleBam{$normal} = $tmp_bf}
                elsif ($tmp_bf =~ /\.realn\./){$sampleBam{$normal} = $tmp_bf}
                elsif ($tmp_bf =~ /\.mkdup\./){$sampleBam{$normal} = $tmp_bf}
                elsif ($tmp_bf =~ /$normal/){$sampleBam{$normal} = $tmp_bf}
            }
        }
        if(!$sampleBam{$tumor}){
            my @MF = glob("$outdir/$tumor/Sentieon/$tumor*.bam");
            for my $tmp_bf (@MF){
                if ($tmp_bf =~ /\.bqsr\./){$sampleBam{$tumor} = $tmp_bf}
                elsif ($tmp_bf =~ /\.realn\./){$sampleBam{$tumor} = $tmp_bf}
                elsif ($tmp_bf =~ /\.mkdup\./){$sampleBam{$tumor} = $tmp_bf}
                elsif ($tmp_bf =~ /$tumor/){$sampleBam{$tumor} = $tmp_bf}
            }
        }        
        ### used already available bam result ###       
        die "Bam file of $normal do not exists !\n" if(!$sampleBam{$normal});
        die "Bam file of $tumor do not exists !\n" if(!$sampleBam{$tumor});
        my $nbam=$sampleBam{$normal};
        my $tbam=$sampleBam{$tumor};
        my $samplePair="$normal-VS-$tumor";
        $NT{$samplePair}=0;
        $hash_sampair{$normal}{$tumor} = $samplePair;
        $hash_sampair{$tumor}{$normal} = $samplePair;
   
        `mkdir -p $outdir/somatic/$samplePair`;
        &SplitTarget;
        &somaticQC($normal,$tumor,$nbam,$tbam);
		if($Use_TNcaller=~/true/i){&TNcallerCalling($normal,$tumor,$nbam,$tbam);}
        if($Use_MuTect=~/true/i){&MuTectCalling($normal,$tumor,$nbam,$tbam);}    
        if($Use_MuTect2=~/true/i){&MuTect2Calling($normal,$tumor,$nbam,$tbam); &mergeBed($normal,$tumor,$nbam,$tbam,"MuTect2","snv_shortindel");}               
        if($Use_MuSE=~/true/i){&MuSECalling($normal,$tumor,$nbam,$tbam); &mergeBed($normal,$tumor,$nbam,$tbam,"MuSE","snv");}       
        if($Use_Strelka2=~/true/i){&Strelka2Calling($normal,$tumor,$nbam,$tbam); &mergeBed($normal,$tumor,$nbam,$tbam,"Strelka2","snv,indel");}
        if($Use_SomaticSniper=~/true/i){&SomaticSniperCalling($normal,$tumor,$nbam,$tbam); &mergeBed($normal,$tumor,$nbam,$tbam,"SomaticSniper","snv");}
        if($Use_Lancet=~/true/i){&LancetCalling($normal,$tumor,$nbam,$tbam); &mergeBed($normal,$tumor,$nbam,$tbam,"Lancet","snv_shortindel");}
        if($Use_Manta=~/true/i){&MantaCalling($normal,$tumor,$nbam,$tbam);}
        if($Use_Svaba=~/true/i){&runSvaba($normal,$tumor,$nbam,$tbam);}    
        if($Use_Platypus=~/true/i){&PlatypusCalling($normal,$tumor,$nbam,$tbam);}
        if($Use_FACETS=~/true/i){&RunFACETS($normal,$tumor,$nbam,$tbam);}
        if($Use_novoBreak=~/true/i){&runnovoBreak($normal,$tumor,$nbam,$tbam);}
        if($mergeSI=~/true/i){&mergeMulTools($normal,$tumor,$nbam,$tbam);} #merge mutli-tools somatic snv/indel results
    }
    close SMTF;
}

############################### monitor to control the pipeline #######################

&edgeList();  ##creat scrips to monitor

#########################################################################################

## date
sub date{
    my $dat='';
    my $tmp=`date`;
    my @tmp=split /\s+/,$tmp;
    $dat.=$tmp[2];
    my @time=split /:/,$tmp[3];
    foreach (@time){$dat.=$_;}
    return $dat;
}

## echo sign for monitor  
sub echostring{
    my $sh=shift;
    my $ostr="echo ==========end at : `date` ========== && \\\n";
    $ostr.="echo Still_waters_run_deep 1>&2 && \\\n";
    $ostr.="echo Still_waters_run_deep > $sh.sign\n";
    return $ostr;
}

## parse targetRegion for gatk
sub parseTR{
    my $tr=shift;
    my $name = (split /\//, $tr)[-1];
    `mkdir -p $outdir/TR`;
    `cp $tr $outdir/TR/$name `;  #`cp $tr $outdir/TR/allchr.bed `;
    `export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH`;
    `cat $outdir/TR/$name |$bgzip -f > $outdir/TR/$name.gz`;
    `$tabix -f -p bed $outdir/TR/$name.gz`;
    open TR,"$tr" or die $!;
    open OUT,">$outdir/TR/$name.intervals" or die $!;  
    while(<TR>){
        chomp;
        my @c=split /\s+/,$_;
        $c[1] += 1;
        print OUT "$c[0]:$c[1]-$c[2]\n";
    }
    close TR;
    close OUT;
    return ($name);
}

## fastq Clean by SOAPnuke
my %fqnum; #check fq input redundancy

sub peLaneClean{
    my ($samp,$fqprefix,$lane_barcode,$q,$id)=@_;
    my ($fq1,$fq2);
    my @fqfiles=`ls $fqprefix\_*.fq*`;
    foreach(@fqfiles){
        chomp;
        if($_=~/\_1.fq$/ || $_=~/\_1.fq.gz$/){$fq1=$_;}
        if($_=~/\_2.fq$/ || $_=~/\_2.fq.gz$/){$fq2=$_;}
    }
    $fqnum{$fq1}++;
    $fqnum{$fq2}++;
    die "Do not exists $fqprefix*.fq.gz\n" unless ( -e $fq1 && -e $fq2);
    die "found redundancy FQ prefix $fqprefix, please check your info file\n" unless ($fqnum{$fq1} eq 1 && $fqnum{$fq2} eq 1);
 
    my $cleandir="$outdir/$samp/$lane_barcode/clean";
    my $fqout1="$lane_barcode\_1.clean.fq.gz";
    my $fqout2="$lane_barcode\_2.clean.fq.gz";
    my $cleanCmd="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    $cleanCmd.="if [ -e \"$cleandir/$fqout1\" ];then rm -rf $cleandir/$fqout1;fi && \\\n";
    $cleanCmd.="if [ -e \"$cleandir/$fqout2\" ];then rm -rf $cleandir/$fqout2;fi && \\\n";
	$cleanCmd.="export LD_LIBRARY_PATH=/share/app/libz/zlib-1.2.11:\$LD_LIBRARY_PATH && \\\n";
    $cleanCmd.="$SOAPnuke filter -1 $fq1 -2 $fq2 $SOAPnukePara -Q $q -G 2 --seqType $id -o $cleandir -C $fqout1 -D $fqout2\n"; 
    $cleanCmd.="gzip -t $cleandir/$fqout1 && \\\n";
	$cleanCmd.="gzip -t $cleandir/$fqout2 && \\\n";
	$cleanCmd.="/share/app/perl-5.22.0/bin/perl $BIN/soapnuke2_stat.pl $cleandir/Basic_Statistics_of_Sequencing_Quality.txt $cleandir/Statistics_of_Filtered_Reads.txt > $cleandir/$lane_barcode.clean.stat && \\\n";
    open CSH,">$outdir/$samp/$lane_barcode/shell/clean.sh";
    $cleanCmd.=echostring("$outdir/$samp/$lane_barcode/shell/clean.sh");
    print CSH $cleanCmd;
    close CSH;
	
	my $checkCmd="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
	$checkCmd.="$fqcheck -r $cleandir/$fqout1 -c $cleandir/$lane_barcode\_1.clean.fqcheck && \\\n";
    $checkCmd.="$fqcheck -r $cleandir/$fqout2 -c $cleandir/$lane_barcode\_2.clean.fqcheck && \\\n";
    $checkCmd.="export GNUPLOT_PS_DIR=/share/app/gnuplot-4.6.7/share/gnuplot/4.6/PostScript && \\\n";
    $checkCmd.="export PERL5LIB=$BIN/../lib:\$PERL5LIB && \\\n";
    $checkCmd.="/share/app/perl-5.22.0/bin/perl $BIN/fqcheck_distribute.pl $cleandir/$lane_barcode\_1.clean.fqcheck $cleandir/$lane_barcode\_2.clean.fqcheck -o $cleandir/$lane_barcode.clean. && \\\n";
    $checkCmd.="rm -rf $cleandir/*.txt $cleandir/log && \\\n";
	open CSH,">$outdir/$samp/$lane_barcode/shell/fqcheck.sh";
    $checkCmd.=echostring("$outdir/$samp/$lane_barcode/shell/fqcheck.sh");
    print CSH $checkCmd;
    close CSH;
}

## clean reads align by Sentieon (BWA-mem)
sub peLaneAln{
    my ($samp,$rgid,$lib,$lane,$lane_barcode,$platform,$fq1,$fq2) = @_;
    open ALNSH, ">$outdir/$samp/$lane_barcode/shell/Sentieon.bwamem.sh" || die $!;
    print ALNSH "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    print ALNSH "$Sentieon bwa mem -M -t $SentieonCPU -R \"\@RG\\tID:$rgid\\tLB:$lib\\tPU:$lane\\tSM:$samp\\tPL:$platform\\tCN:BGI\" $reference $fq1 $fq2 -o $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.bam && \\\n";
    print ALNSH "$SentieonAuto driver -t $SentieonCPU -r $reference -i $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.bam --algo GCBias --summary $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.total.GC_SUMMARY.txt $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.total.GC_METRIC.txt --algo MeanQualityByCycle $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.total.MQ_METRIC.txt --algo QualDistribution $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.total.QD_METRIC.txt --algo InsertSizeMetricAlgo $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.total.IS_METRIC.txt --algo AlignmentStat $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.total.ALN_METRIC.txt --algo WgsMetricsAlgo $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.total.WGS_METRIC.txt --algo SequenceArtifactMetricsAlgo --dbsnp $Dbsnp $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.total.ARTIFACT_METRICS_BASE.txt && \\\n";
    print ALNSH "$SentieonAuto plot GCBias -o $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.total.GC_METRIC.pdf $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.total.GC_METRIC.txt && \\\n";
    print ALNSH "$SentieonAuto plot MeanQualityByCycle -o $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.total.MQ_METRIC.pdf $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.total.MQ_METRIC.txt && \\\n";
    print ALNSH "$SentieonAuto plot QualDistribution -o $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.total.QD_METRIC.pdf $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.total.QD_METRIC.txt && \\\n";
    print ALNSH "$SentieonAuto plot InsertSizeMetricAlgo -o $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.total.IS_METRIC.pdf $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.total.IS_METRIC.txt && \\\n";
    if($seqType =~ /target/i){ # WES & Panel
        print ALNSH "$SentieonAuto driver -t $SentieonCPU --interval $outdir/TR/$trname -r $reference -i $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.bam --algo GCBias --summary $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.target.GC_SUMMARY.txt $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.target.GC_METRIC.txt --algo MeanQualityByCycle $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.target.MQ_METRIC.txt --algo QualDistribution $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.target.QD_METRIC.txt --algo InsertSizeMetricAlgo $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.target.IS_METRIC.txt --algo AlignmentStat $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.target.ALN_METRIC.txt && \\\n";
        print ALNSH "$SentieonAuto plot GCBias -o $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.target.GC_METRIC.pdf $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.target.GC_METRIC.txt && \\\n";
        print ALNSH "$SentieonAuto plot MeanQualityByCycle -o $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.target.MQ_METRIC.pdf $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.target.MQ_METRIC.txt && \\\n";
        print ALNSH "$SentieonAuto plot QualDistribution -o $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.target.QD_METRIC.pdf $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.target.QD_METRIC.txt && \\\n";
        print ALNSH "$SentieonAuto plot InsertSizeMetricAlgo -o $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.target.IS_METRIC.pdf $outdir/$samp/$lane_barcode/Sentieon/$lane_barcode.sorted.target.IS_METRIC.txt && \\\n";
    }
    print ALNSH echostring("$outdir/$samp/$lane_barcode/shell/Sentieon.bwamem.sh");
    close ALNSH;
}

##################
if($bamList){    
    open(EB,">$outdir/TR/shell/pipeline_with_existed_bam.sh") or die $!;
    my $EBSH="echo this is only used for somatic mutation calling with existed sample pairs && \\\n";
    $EBSH.=echostring("$outdir/TR/shell/pipeline_with_existed_bam.sh");
    print EB $EBSH;
    close EB;
}
#################

my %TRs;

#target number for each samplepair and corresponding software
my %Pair_TRs;
sub SplitTarget{
    `mkdir -p $outdir/TR/shell`;
    if(! -e "$outdir/TR/soft_TRs.record.txt"){
        open(RD,">$outdir/TR/soft_TRs.record.txt") or die "not target records\n";
        my $time=`date`;
        chomp $time;
        $time="== $time ==";
        foreach my $pair(sort keys %NT){    
            foreach my $soft(keys %SoftTRs){
                next unless (exists $SoftTRs{$soft});
                $TRs{$SoftTRs{$soft}}=0;
                print RD "$pair\t$soft\t$SoftTRs{$soft}\t$time\n";                
            }
        }            
    close RD;    
    }
    else{
        my %REC;
        open(RD,"<$outdir/TR/soft_TRs.record.txt") or die "not record file\n";
        while(my $record = <RD>){
            chomp $record;
            my @ar=split/\t/,$record;
            my($pair,$soft,$num,$time)=($ar[0],$ar[1],$ar[2],$ar[3]);
            $REC{$soft}{$pair}=$num;
        }
        close RD;   
        my $check_num=0; 
        open(CONFLICT,">$outdir/TR/soft_TRs.conflict.txt") or die $!;
        foreach my $nt(sort keys %NT){
            foreach my $soft(sort keys %SoftTRs){
                next unless (exists $SoftTRs{$soft});
                if(exists $REC{$soft}{$nt}){
                    if ($REC{$soft}{$nt} ne $SoftTRs{$soft}){
                        my @files=glob("$outdir/somatic/$nt/$soft/*/*");
                        if($files[0]){
                            $check_num++;
                            print CONFLICT "the $soft target number for $nt was $REC{$soft}{$nt}, but now is $SoftTRs{$soft}\n";                       
                        }
                    }
                }
            }
        }
        close CONFLICT;
        if($check_num > 0){
            print "target regoin number for some sample pairs are conflict with previous one. if still need to change that, please enter \<YES\>, and delete all temporary files, including bash, raw and result files, if \<NO\>,please change split number as previous one in your configure file\n";
            my $answer = <STDIN>;
            chomp $answer;
            if($answer =~ /YES/i){
                die "please seen conflict samples and corresponding software detail files at $outdir/TR/soft_TRs.conflict.txt\n";
            }
            else{
                die "change your software target split number, see previous records details at $outdir/TR/soft_TRs.record.txt\n";
            }
        }        
        if($check_num eq 0){
            `rm -fr $outdir/TR/soft_TRs.conflict.txt`;
            open(RD,">>$outdir/TR/soft_TRs.record.txt") or die "not target records\n";
            foreach my $nt(sort keys %NT){
                foreach my $soft(sort keys %SoftTRs){
                    my $time=`date`;
                    chomp $time;
                    $time="== $time ==";
                    foreach my $pair(sort keys %NT){
                        foreach my $soft(keys %SoftTRs){
                            $TRs{$SoftTRs{$soft}}=0;
                            if(exists $REC{$soft}{$pair} && $REC{$soft}{$pair} eq $SoftTRs{$soft}){
                            }
                            else{
                                print RD "$pair\t$soft\t$SoftTRs{$soft}\t$time\n";
                            }   
                        }
                    }
                }
            }
        }    
        close RD; 
    }                                       
    foreach my $TRsplitNum(sort keys %TRs){
        open(SP,">$outdir/TR/shell/split_intervals.$TRsplitNum.sh") or die "$!\n";
        my $splitSH="";         
       
        if ($TRsplitNum eq 1){
            $splitSH.="echo hello \\\n";
        }
        elsif($TRsplitNum > 1){ 
            `mkdir -p $outdir/TR/split_intervals$TRsplitNum`;       
            $splitSH.="rm -fr $outdir/TR/split_intervals$TRsplitNum/* && \\\n";
            $splitSH.="export PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/JAVA/jre1.8.0_92/bin:\$PATH && \\\n";
            $splitSH.="$GATK4 SplitIntervals -R $reference -L $targetRegion --scatter-count $TRsplitNum -O $outdir/TR/split_intervals$TRsplitNum/ && \\\n";
            $splitSH.="ls $outdir/TR/split_intervals$TRsplitNum/|/share/app/perl-5.22.0/bin/perl -we 'while(<>){chomp;\$_=~/^0+([0-9]+-scattered.interval_list)/;my \$tr=\$1;`mv $outdir/TR/split_intervals$TRsplitNum/\$_ $outdir/TR/split_intervals$TRsplitNum/\$tr.bed`;}' && \\\n";
            my $Num=$TRsplitNum-1;
            $splitSH.="export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH && \\\n";
            $splitSH.="for i in {0..$Num};do cat $outdir/TR/split_intervals$TRsplitNum/\$i-scattered.interval_list.bed|grep -v \"\@\"|awk '{print \$1\"\\t\"\$2-1\"\\t\"\$3}'|$bgzip > $outdir/TR/split_intervals$TRsplitNum/\$i-scattered.intervals.bed.gz;$tabix -f $outdir/TR/split_intervals$TRsplitNum/\$i-scattered.intervals.bed.gz;done  && \\\n";
            $splitSH.="for i in {0..$Num};do zcat $outdir/TR/split_intervals$TRsplitNum/\$i-scattered.intervals.bed.gz > $outdir/TR/split_intervals$TRsplitNum/\$i-scattered.intervals.bed;done && \\\n"; 
            $splitSH.="rm $outdir/TR/split_intervals$TRsplitNum/*interval_list.bed && \\\n";
        }
        $splitSH.=echostring("$outdir/TR/shell/split_intervals.$TRsplitNum.sh");
        print SP $splitSH;
        close SP;
    }       
} 

## Merge, indel-realignment, BQSR, QC and SNP/INDEL Calling by Sentieon, then readcounter
sub RunGermline{
    my ($samp, $bamarray) = @_;

    # merge multiple bams 
    my $instring = join " -i ", @$bamarray;
    my $rmstring = "";
    for my $tmpfile (@$bamarray){
	    $rmstring .= "$tmpfile $tmpfile.bai ";
    }
	if($oneSentieon=~/false/i){
    	open MSH, ">$outdir/$samp/shell/Sentieon.mkdup.$samp.sh" || die $!;
    	print MSH "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    	print MSH "$Sentieon driver -t $SentieonCPU -i $instring --algo LocusCollector --fun score_info $outdir/$samp/Sentieon/$samp.mkdup.score.txt.gz && \\\n";
    	print MSH "$Sentieon driver -t $SentieonCPU -i $instring --algo Dedup --score_info $outdir/$samp/Sentieon/$samp.mkdup.score.txt.gz --metrics $outdir/$samp/Sentieon/$samp.mkdup.metrics.txt $outdir/$samp/Sentieon/$samp.mkdup.bam && \\\n";
    	print MSH "rm -rf $rmstring && \\\n";
   		print MSH echostring("$outdir/$samp/shell/Sentieon.mkdup.$samp.sh");
    	close MSH;

    	# indel-realignment
    	open IRSH, ">$outdir/$samp/shell/Sentieon.realn.$samp.sh" || die $!;
    	print IRSH "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    	print IRSH "$SentieonAuto driver -t $SentieonCPU -r $reference -i $outdir/$samp/Sentieon/$samp.mkdup.bam --algo Realigner -k $kgIndel -k $millIndel";
    	if($seqType =~ /target/i){ # WES & Panel
        	print IRSH " --interval_list $outdir/TR/$trname";
    	}
    	print IRSH " $outdir/$samp/Sentieon/$samp.realn.bam && \\\n";
    	print IRSH "rm -rf $outdir/$samp/Sentieon/$samp.mkdup.bam* $outdir/$samp/Sentieon/$samp.mkdup.score.txt.gz* $outdir/$samp/Sentieon/$samp.mkdup.metrics.txt && \\\n";
    	print IRSH echostring("$outdir/$samp/shell/Sentieon.realn.$samp.sh");
    	close IRSH;

    	# bqsr & snp, indel calling
    	open BRSH, ">$outdir/$samp/shell/Sentieon.bqsr.$samp.sh" || die $!;
    	print BRSH "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    	print BRSH "$SentieonAuto driver -t $SentieonCPU -r $reference -i $outdir/$samp/Sentieon/$samp.realn.bam --algo QualCal -k $Dbsnp -k $kgIndel -k $millIndel $outdir/$samp/Sentieon/$samp.bqsr.RECAL_PRE.txt && \\\n";
    	print BRSH "$SentieonAuto driver -t $SentieonCPU -r $reference -i $outdir/$samp/Sentieon/$samp.realn.bam -q $outdir/$samp/Sentieon/$samp.bqsr.RECAL_PRE.txt --algo QualCal -k $Dbsnp -k $kgIndel -k $millIndel $outdir/$samp/Sentieon/$samp.bqsr.RECAL_POST.txt --algo ReadWriter $outdir/$samp/Sentieon/$samp.bqsr.bam --algo Haplotyper -d $Dbsnp  $outdir/$samp/Sentieon/$samp.Haplotyper.vcf.gz && \\\n";
    	print BRSH "$SentieonAuto driver -t $SentieonCPU --algo QualCal --plot --before $outdir/$samp/Sentieon/$samp.bqsr.RECAL_PRE.txt --after $outdir/$samp/Sentieon/$samp.bqsr.RECAL_POST.txt $outdir/$samp/Sentieon/$samp.bqsr.RECAL_RESULT.csv && \\\n";
    	print BRSH "$SentieonAuto plot QualCal -o $outdir/$samp/Sentieon/$samp.bqsr.RECAL_RESULT.pdf $outdir/$samp/Sentieon/$samp.bqsr.RECAL_RESULT.csv && \\\n";
    	print BRSH "rm -rf $outdir/$samp/Sentieon/$samp.realn.bam* && \\\n";
    	print BRSH echostring("$outdir/$samp/shell/Sentieon.bqsr.$samp.sh");
    	close BRSH;

    	$sampleBam{$samp} = "$outdir/$samp/Sentieon/$samp.bqsr.bam";

    	# qc by Sentieon
    	open QCSH, ">$outdir/$samp/shell/Sentieon.qc.$samp.sh" || die $!;
    	print QCSH "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    	print QCSH "$SentieonAuto driver -t $SentieonCPU -r $reference -i $outdir/$samp/Sentieon/$samp.bqsr.bam --algo GCBias --summary $outdir/$samp/Sentieon/$samp.bqsr.total.GC_SUMMARY.txt $outdir/$samp/Sentieon/$samp.bqsr.total.GC_METRIC.txt --algo MeanQualityByCycle $outdir/$samp/Sentieon/$samp.bqsr.total.MQ_METRIC.txt --algo QualDistribution $outdir/$samp/Sentieon/$samp.bqsr.total.QD_METRIC.txt --algo InsertSizeMetricAlgo $outdir/$samp/Sentieon/$samp.bqsr.total.IS_METRIC.txt --algo AlignmentStat $outdir/$samp/Sentieon/$samp.bqsr.total.ALN_METRIC.txt --algo WgsMetricsAlgo $outdir/$samp/Sentieon/$samp.bqsr.total.WGS_METRIC.txt --algo SequenceArtifactMetricsAlgo --dbsnp $Dbsnp $outdir/$samp/Sentieon/$samp.bqsr.total.ARTIFACT_METRICS && \\\n";
    	print QCSH "$SentieonAuto plot GCBias -o $outdir/$samp/Sentieon/$samp.bqsr.total.GC_METRIC.pdf $outdir/$samp/Sentieon/$samp.bqsr.total.GC_METRIC.txt && \\\n";
    	print QCSH "$SentieonAuto plot MeanQualityByCycle -o $outdir/$samp/Sentieon/$samp.bqsr.total.MQ_METRIC.pdf $outdir/$samp/Sentieon/$samp.bqsr.total.MQ_METRIC.txt && \\\n";
    	print QCSH "$SentieonAuto plot QualDistribution -o $outdir/$samp/Sentieon/$samp.bqsr.total.QD_METRIC.pdf $outdir/$samp/Sentieon/$samp.bqsr.total.QD_METRIC.txt && \\\n";
    	print QCSH "$SentieonAuto plot InsertSizeMetricAlgo -o $outdir/$samp/Sentieon/$samp.bqsr.total.IS_METRIC.pdf $outdir/$samp/Sentieon/$samp.bqsr.total.IS_METRIC.txt && \\\n";
    	if($seqType =~ /target/i){ # WES & Panel
       		print QCSH "$SentieonAuto driver -t $SentieonCPU --interval $outdir/TR/$trname -r $reference -i $outdir/$samp/Sentieon/$samp.bqsr.bam --algo GCBias --summary $outdir/$samp/Sentieon/$samp.bqsr.target.GC_SUMMARY.txt $outdir/$samp/Sentieon/$samp.bqsr.target.GC_METRIC.txt --algo MeanQualityByCycle $outdir/$samp/Sentieon/$samp.bqsr.target.MQ_METRIC.txt --algo QualDistribution $outdir/$samp/Sentieon/$samp.bqsr.target.QD_METRIC.txt --algo InsertSizeMetricAlgo $outdir/$samp/Sentieon/$samp.bqsr.target.IS_METRIC.txt --algo AlignmentStat $outdir/$samp/Sentieon/$samp.bqsr.target.ALN_METRIC.txt && \\\n";
        	print QCSH "$SentieonAuto plot GCBias -o $outdir/$samp/Sentieon/$samp.bqsr.target.GC_METRIC.pdf $outdir/$samp/Sentieon/$samp.bqsr.target.GC_METRIC.txt && \\\n";
        	print QCSH "$SentieonAuto plot MeanQualityByCycle -o $outdir/$samp/Sentieon/$samp.bqsr.target.MQ_METRIC.pdf $outdir/$samp/Sentieon/$samp.bqsr.target.MQ_METRIC.txt && \\\n";
        	print QCSH "$SentieonAuto plot QualDistribution -o $outdir/$samp/Sentieon/$samp.bqsr.target.QD_METRIC.pdf $outdir/$samp/Sentieon/$samp.bqsr.target.QD_METRIC.txt && \\\n";
        	print QCSH "$SentieonAuto plot InsertSizeMetricAlgo -o $outdir/$samp/Sentieon/$samp.bqsr.target.IS_METRIC.pdf $outdir/$samp/Sentieon/$samp.bqsr.target.IS_METRIC.txt && \\\n";
    	}
    	print QCSH echostring("$outdir/$samp/shell/Sentieon.qc.$samp.sh");
    	close QCSH;
	}
    elsif($oneSentieon=~/true/i){
        open MSH, ">$outdir/$samp/shell/Sentieon.mkdup.realign.bqsr.qc.$samp.sh" || die $!;
        print MSH "#!/bin/bash\nset -euo pipefail\necho ==========mkdup start at : `date` ==========\n";
        print MSH "$Sentieon driver -t $SentieonCPU -i $instring --algo LocusCollector --fun score_info $outdir/$samp/Sentieon/$samp.mkdup.score.txt.gz && \\\n";
        print MSH "$Sentieon driver -t $SentieonCPU -i $instring --algo Dedup --score_info $outdir/$samp/Sentieon/$samp.mkdup.score.txt.gz --metrics $outdir/$samp/Sentieon/$samp.mkdup.metrics.txt $outdir/$samp/Sentieon/$samp.mkdup.bam && \\\n";
        print MSH "echo ==========mkdup end at : `date` ==========\n";
        print MSH "echo ==========indel realignment start at : `date` ==========\n"; 
        print MSH "$SentieonAuto driver -t $SentieonCPU -r $reference -i $outdir/$samp/Sentieon/$samp.mkdup.bam --algo Realigner -k $kgIndel -k $millIndel";
        if($seqType =~ /target/i){ # WES & Panel
            print MSH " --interval_list $outdir/TR/$trname";
        }
        print MSH " $outdir/$samp/Sentieon/$samp.realn.bam && \\\n";
        print MSH "echo ==========indel realignment end at : `date` ==========\n";
        print MSH "echo ==========bqsr, snp and indel calling start at : `date` ==========\n";
        print MSH "$SentieonAuto driver -t $SentieonCPU -r $reference -i $outdir/$samp/Sentieon/$samp.realn.bam --algo QualCal -k $Dbsnp -k $kgIndel -k $millIndel $outdir/$samp/Sentieon/$samp.bqsr.RECAL_PRE.txt && \\\n";
        print MSH "$SentieonAuto driver -t $SentieonCPU -r $reference -i $outdir/$samp/Sentieon/$samp.realn.bam -q $outdir/$samp/Sentieon/$samp.bqsr.RECAL_PRE.txt --algo QualCal -k $Dbsnp -k $kgIndel -k $millIndel $outdir/$samp/Sentieon/$samp.bqsr.RECAL_POST.txt --algo ReadWriter $outdir/$samp/Sentieon/$samp.bqsr.bam --algo Haplotyper -d $Dbsnp  $outdir/$samp/Sentieon/$samp.Haplotyper.vcf.gz && \\\n";
        print MSH "$SentieonAuto driver -t $SentieonCPU --algo QualCal --plot --before $outdir/$samp/Sentieon/$samp.bqsr.RECAL_PRE.txt --after $outdir/$samp/Sentieon/$samp.bqsr.RECAL_POST.txt $outdir/$samp/Sentieon/$samp.bqsr.RECAL_RESULT.csv && \\\n";
        print MSH "$SentieonAuto plot QualCal -o $outdir/$samp/Sentieon/$samp.bqsr.RECAL_RESULT.pdf $outdir/$samp/Sentieon/$samp.bqsr.RECAL_RESULT.csv && \\\n";
        $sampleBam{$samp} = "$outdir/$samp/Sentieon/$samp.bqsr.bam";
        print MSH "echo ==========bqsr, snp and indel calling end at : `date` ==========\n";
        print MSH "echo ==========QC start at : `date` ==========\n";
        print MSH "$SentieonAuto driver -t $SentieonCPU -r $reference -i $outdir/$samp/Sentieon/$samp.bqsr.bam --algo GCBias --summary $outdir/$samp/Sentieon/$samp.bqsr.total.GC_SUMMARY.txt $outdir/$samp/Sentieon/$samp.bqsr.total.GC_METRIC.txt --algo MeanQualityByCycle $outdir/$samp/Sentieon/$samp.bqsr.total.MQ_METRIC.txt --algo QualDistribution $outdir/$samp/Sentieon/$samp.bqsr.total.QD_METRIC.txt --algo InsertSizeMetricAlgo $outdir/$samp/Sentieon/$samp.bqsr.total.IS_METRIC.txt --algo AlignmentStat $outdir/$samp/Sentieon/$samp.bqsr.total.ALN_METRIC.txt --algo WgsMetricsAlgo $outdir/$samp/Sentieon/$samp.bqsr.total.WGS_METRIC.txt --algo SequenceArtifactMetricsAlgo --dbsnp $Dbsnp $outdir/$samp/Sentieon/$samp.bqsr.total.ARTIFACT_METRICS && \\\n";
        print MSH "$SentieonAuto plot GCBias -o $outdir/$samp/Sentieon/$samp.bqsr.total.GC_METRIC.pdf $outdir/$samp/Sentieon/$samp.bqsr.total.GC_METRIC.txt && \\\n";
        print MSH "$SentieonAuto plot MeanQualityByCycle -o $outdir/$samp/Sentieon/$samp.bqsr.total.MQ_METRIC.pdf $outdir/$samp/Sentieon/$samp.bqsr.total.MQ_METRIC.txt && \\\n";
        print MSH "$SentieonAuto plot QualDistribution -o $outdir/$samp/Sentieon/$samp.bqsr.total.QD_METRIC.pdf $outdir/$samp/Sentieon/$samp.bqsr.total.QD_METRIC.txt && \\\n";
        print MSH "$SentieonAuto plot InsertSizeMetricAlgo -o $outdir/$samp/Sentieon/$samp.bqsr.total.IS_METRIC.pdf $outdir/$samp/Sentieon/$samp.bqsr.total.IS_METRIC.txt && \\\n";
        if($seqType =~ /target/i){ # WES & Panel
            print MSH "$SentieonAuto driver -t $SentieonCPU --interval $outdir/TR/$trname -r $reference -i $outdir/$samp/Sentieon/$samp.bqsr.bam --algo GCBias --summary $outdir/$samp/Sentieon/$samp.bqsr.target.GC_SUMMARY.txt $outdir/$samp/Sentieon/$samp.bqsr.target.GC_METRIC.txt --algo MeanQualityByCycle $outdir/$samp/Sentieon/$samp.bqsr.target.MQ_METRIC.txt --algo QualDistribution $outdir/$samp/Sentieon/$samp.bqsr.target.QD_METRIC.txt --algo InsertSizeMetricAlgo $outdir/$samp/Sentieon/$samp.bqsr.target.IS_METRIC.txt --algo AlignmentStat $outdir/$samp/Sentieon/$samp.bqsr.target.ALN_METRIC.txt && \\\n";
            print MSH "$SentieonAuto plot GCBias -o $outdir/$samp/Sentieon/$samp.bqsr.target.GC_METRIC.pdf $outdir/$samp/Sentieon/$samp.bqsr.target.GC_METRIC.txt && \\\n";
            print MSH "$SentieonAuto plot MeanQualityByCycle -o $outdir/$samp/Sentieon/$samp.bqsr.target.MQ_METRIC.pdf $outdir/$samp/Sentieon/$samp.bqsr.target.MQ_METRIC.txt && \\\n";
            print MSH "$SentieonAuto plot QualDistribution -o $outdir/$samp/Sentieon/$samp.bqsr.target.QD_METRIC.pdf $outdir/$samp/Sentieon/$samp.bqsr.target.QD_METRIC.txt && \\\n";
            print MSH "$SentieonAuto plot InsertSizeMetricAlgo -o $outdir/$samp/Sentieon/$samp.bqsr.target.IS_METRIC.pdf $outdir/$samp/Sentieon/$samp.bqsr.target.IS_METRIC.txt && \\\n";
        }
        print MSH "rm -rf $rmstring && \\\n";
        print MSH "rm -rf $outdir/$samp/Sentieon/$samp.mkdup.bam* $outdir/$samp/Sentieon/$samp.mkdup.score.txt.gz* $outdir/$samp/Sentieon/$samp.mkdup.metrics.txt && \\\n";
        print MSH "rm -rf $outdir/$samp/Sentieon/$samp.realn.bam* && \\\n";
        print MSH "echo ==========QC end at : `date` ==========\n";
        print MSH echostring("$outdir/$samp/shell/Sentieon.mkdup.realign.bqsr.qc.$samp.sh");
        close MSH;
    }

    open ANNSH, ">$outdir/$samp/shell/variant_Ann.$samp.sh" || die $!;
    print ANNSH "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    print ANNSH "export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH && \\\n";
    print ANNSH "/share/app/perl-5.22.0/bin/perl $annovar/table_annovar.pl $outdir/$samp/Sentieon/$samp.Haplotyper.vcf.gz $annovardb -buildver hg38 -out $outdir/$samp/Sentieon/$samp.Haplotyper -remove -protocol refGene,ensGene,knownGene,cytoBand,genomicSuperDups,gwasCatalog,phastConsElements100way,wgRna,avsnp150,clinvar_20190305,cosmic89,dbnsfp35c,exac03nonpsych,exac03nontcga,gnomad_exome,gnomad_genome,intervar_20180118 -operation g,g,g,r,r,r,r,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput && \\\n";
    print ANNSH "gzip -f $outdir/$samp/Sentieon/$samp.Haplotyper.hg38_multianno.txt && \\\n";
    print ANNSH "rm -rf $outdir/$samp/Sentieon/$samp.Haplotyper.hg38_multianno.vcf $outdir/$samp/Sentieon/$samp.Haplotyper.avinput && \\\n";
    print ANNSH echostring("$outdir/$samp/shell/variant_Ann.$samp.sh");
    close ANNSH;
	
    #read count per 1000bp
    my $outwig = $sampleBam{$samp};
    $outwig =~ s/\.bam$/.win1k.wig.gz/;
    open RCSH, ">$outdir/$samp/shell/bam_readcount.$samp.sh" || die $!;
    print RCSH "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    if (`grep chr $reffai`){
        print RCSH "$BIN/readCounter -w 1000 -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY $sampleBam{$samp} | gzip > $outwig && \\\n";
    }
    else{
        print RCSH "$BIN/readCounter -w 1000 -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y $sampleBam{$samp} | gzip > $outwig && \\\n";
    }
    print RCSH echostring("$outdir/$samp/shell/bam_readcount.$samp.sh");
    close RCSH;

=header # PanCanQC not support hg38
    # calculate the quality control values used in the ICGC PanCan project. Generally for WGS data.
    # modify form https://github.com/eilslabs/PanCanQC
    my $qcdir = "$outdir/$samp/QC";
    my $ACEseq = "$outdir/$samp/QC/ACEseq";
    my $qcshell = "$outdir/$samp/QC/shell";
    unless(-e $qcdir){system("mkdir -p $qcdir")}
    unless(-e $ACEseq){system("mkdir -p $ACEseq")}
    unless(-e $qcshell){system("mkdir -p $qcshell")}
    open STEP1, ">$qcshell/flagstat.$samp.sh" || die $!;
    print STEP1 "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    print STEP1 "export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH && \\\n";
    print STEP1 "$samtools flagstat $sampleBam{$samp} > $qcdir/$samp.flagstat && \\\n";
    print STEP1 echostring("$qcshell/flagstat.$samp.sh");
    close STEP1;
	
    open STEP2, ">$qcshell/coverageQc.$samp.sh" || die $!;
    print STEP2 "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    print STEP2 "$PanCanQC/coverageQc --alignmentFile=$sampleBam{$samp} --outputFile=$qcdir/$samp.genome_coverage --processors=1 --basequalCutoff=0 --ungappedSizes=$PanCanQC/hg19.chrLenOnlyACGT_realChromosomes.tab && \\\n";
    print STEP2 echostring("$qcshell/coverageQc.$samp.sh");
    close STEP2;
	
    open STEP3, ">$qcshell/genomeCoverage.$samp.sh" || die $!;
    print STEP3 "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    print STEP3 "export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:/opt/python/lib:\$LD_LIBRARY_PATH && \\\n";
    print STEP3 "$PanCanQC/genomeCoverage --alignmentFile=$sampleBam{$samp} --outputFile=/dev/stdout --processors=4 --mode=countReads --windowSize=1 | /share/app/perl-5.22.0/bin/perl $PanCanQC/filter_readbins.pl - $PanCanQC/hg19.chrLenOnlyACGT_realChromosomes.tab | gzip > $qcdir/$samp.readbin_coverage.gz && \\\n";  
    print STEP3 "gzip -dc $qcdir/$samp.readbin_coverage.gz | awk '{print \$1,\$2,\$2+999,\$3}' | sed 's/ /\\t/g' |  sed '1i\\#chr\\tpos\\tend\\tcoverage' | /share/app/perl-5.22.0/bin/perl $PanCanQC/annotate_vcf.pl -a - --aFileType=custom --aChromColumn chr --aPosColumn pos --aEndColumn end -b $PanCanQC/wgEncodeCrgMapabilityAlign100mer_chr.bedGraph.gz --tabix_bin $tabix --bFileType=bed --reportBFeatCoord --columnName map | $pythonBin/python $PanCanQC/addMappability.py -o $ACEseq/$samp.readbin_coverage.Mappability.gz && \\\n";
    print STEP3 "$pythonBin/python $PanCanQC/merge_and_filter_cnv.py --inputfile $ACEseq/$samp.readbin_coverage.Mappability.gz --output $ACEseq/$samp.readbin_coverage.Mappability.filtered.gz --coverage 0 --mappability 1000 --NoOfWindows 5 && \\\n";
    print STEP3 "/share/app/R-3.3.2/bin/Rscript $PanCanQC/correctGCBias.R --windowFile $ACEseq/$samp.readbin_coverage.Mappability.filtered.gz --timefile $PanCanQC/ReplicationTime_10cellines_mean_10KB.Rda --chrLengthFile $PanCanQC/chrlengths.txt --pid $samp --outfile $samp.corrected.txt --corPlot $samp.gc_corrected.png --corTab $samp.qc_gc_corrected.tsv --qcTab $samp.qc_gc_corrected.slim.txt --gcFile $PanCanQC/hg19_GRch37_100genomes_gc_content_10kb.txt --outDir $ACEseq --lowess_f 0.1 --scaleFactor 0.9 --coverageYlims 4 && \\\n";
    print STEP3 echostring("$qcshell/genomeCoverage.$samp.sh");
    close STEP3;
	
    open STEP4, ">$qcshell/bam_stats.$samp.sh" || die $!;
    print STEP4 "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    print STEP4 "/share/app/glibc-2.17/lib/ld-2.17.so --library-path /share/app/glibc-2.17/lib:/share/app/libz/zlib-1.2.11 $PanCanQC/bam_stats -i $sampleBam{$samp} -o $qcdir/$samp.read_edits && \\\n";
    print STEP4 echostring("$qcshell/bam_stats.$samp.sh");
    close STEP4;
	
    open STEP5, ">$qcshell/summary.$samp.sh" || die $!;
    print STEP5 "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    print STEP5 "export PERL5LIB=$PanCanQC/perl_modules/lib/site_perl/5.22.0:\$PERL5LIB && \\\n";
    print STEP5 "/share/app/perl-5.22.0/bin/perl $PanCanQC/summaryQC.pl $qcdir/$samp.flagstat $qcdir/$samp.genome_coverage $ACEseq/$samp.qc_gc_corrected.slim.txt $qcdir/$samp.read_edits $qcdir/$samp.readbin_coverage.gz $PanCanQC/hg19_reducedGenome.n300l5M.sorted.bed $qcdir/$samp.PanCanQC.summary.txt && \\\n";
    print STEP5 "tar -czf $qcdir/ACEseq.tar.gz --directory=$qcdir ACEseq && rm -rf $ACEseq && \\\n";
    print STEP5 echostring("$qcshell/summary.$samp.sh");
    close STEP5;
=cut

    # rapidly determine whether two BAM files represent samples from the same biological source by comparing their genotypes
    # the result maybe different with ContEst by GATK3 in soamticQC (maybe cause by different population allele frequencies source)!!!!
    my $qcdir = "$outdir/$samp/QC";
    my $qcshell = "$outdir/$samp/QC/shell";
    unless(-e $qcdir){system("mkdir -p $qcdir")}
    unless(-e $qcshell){system("mkdir -p $qcshell")}
    open STEP6, ">$qcshell/ContEst.$samp.sh" || die $!;
    print STEP6 "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    print STEP6 "mkdir -p $outdir/$samp/QC/tmp && \\\n";
    print STEP6 "export PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/JAVA/jre1.8.0_92/bin:\$PATH && \\\n";
    print STEP6 "$GATK4 --java-options '-Xmx4g -XX:ParallelGCThreads=1 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$outdir/$samp/QC/tmp' GetPileupSummaries -I $sampleBam{$samp} -V $ExACcommon -L $ExACcommon -O $outdir/$samp/QC/$samp.getpileupsummaries.txt && \\\n";
    print STEP6 "$GATK4 --java-options '-Xmx4g -XX:ParallelGCThreads=1 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$outdir/$samp/QC/tmp' CalculateContamination -I $outdir/$samp/QC/$samp.getpileupsummaries.txt -O $outdir/$samp/QC/$samp.calculatecontamination.txt && \\\n";
    print STEP6 "gzip -f $outdir/$samp/QC/$samp.getpileupsummaries.txt && \\\n";
    print STEP6 "rm -rf $outdir/$samp/QC/tmp && \\\n";
    print STEP6 echostring("$qcshell/ContEst.$samp.sh");
    close STEP6;

    # WES & Panel, calculate coverage and depth information
	if ($seqType =~ /target/i){ 
		`mkdir -p $outdir/$samp/Coverage`;
		open COVDEPSH,">$outdir/$samp/shell/covdep.$samp.sh" || die $!;
		print COVDEPSH "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
		print COVDEPSH "export LD_LIBRARY_PATH=$BIN/bamtools-2.5.1/lib64:/share/app/gcc-4.9.3/lib64:/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH && \\\n";
		print COVDEPSH "$BIN/bamdst $covdepPara -p $outdir/TR/$trname -o $outdir/$samp/Coverage $sampleBam{$samp} && \\\n";
        print COVDEPSH "perl $BIN/bamdstPlot.pl -i $outdir/$samp/Coverage/depth_distribution.plot -c $outdir/$samp/Coverage/coverage.report -o $outdir/$samp/Coverage && \\\n";
        print COVDEPSH "$samtools view -bu -F 4 $sampleBam{$samp} -L $outdir/TR/$trname | $bamtools stats -in /dev/stdin -insert > $outdir/$samp/Coverage/Target.report && \\\n";
		print COVDEPSH "gzip -f $outdir/$samp/Coverage/depth_distribution.plot $outdir/$samp/Coverage/insertsize.plot $outdir/$samp/Coverage/uncover.bed && \\\n";
		print COVDEPSH "rm -rf $outdir/$samp/Coverage/depth.tsv.gz && \\\n";
		print COVDEPSH echostring("$outdir/$samp/shell/covdep.$samp.sh");
        close COVDEPSH;
	}

    # remove clean reads
    if ($rmCleanReads=~/true/i){
        open RMSH, ">$outdir/$samp/shell/rmCleanReads.$samp.sh" || die $!;
        print RMSH "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
        print RMSH "/share/app/perl-5.22.0/bin/perl $BIN/rm_cleanReads_Sentieon_pipeline.pl $sampleBam{$samp} $outdir/$samp/Sentieon/01bamList4Germline.csv && \\\n";
        print RMSH echostring("$outdir/$samp/shell/rmCleanReads.$samp.sh");
        close RMSH;
    }

    # convert bam to cram
    if($bam2cram=~/true/i){
        open CVSH, ">$outdir/$samp/shell/Sentieon.bam2cram.$samp.sh" || die $!;
        print CVSH "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
        print CVSH "$SentieonAuto driver -t $SentieonCPU -r $reference -i $outdir/$samp/Sentieon/$samp.bqsr.bam --algo ReadWriter $outdir/$samp/Sentieon/$samp.bqsr.cram && \\\n";
        print CVSH "$SentieonAuto driver -t $SentieonCPU -r $reference -i $outdir/$samp/Sentieon/$samp.bqsr.cram --algo AlignmentStat $outdir/$samp/Sentieon/$samp.bqsr.cram.ALN_METRIC.txt && \\\n";
        print CVSH "LineNum=`grep -wFf $outdir/$samp/Sentieon/$samp.bqsr.total.ALN_METRIC.txt $outdir/$samp/Sentieon/$samp.bqsr.cram.ALN_METRIC.txt | wc -l` && \\\n";
        print CVSH "if [ \$LineNum -eq 6 ];then echo bam was converted to cram successfully and deleted;rm -rf $outdir/$samp/Sentieon/$samp.bqsr.bam $outdir/$samp/Sentieon/$samp.bqsr.bam.bai $outdir/$samp/Sentieon/$samp.bqsr.cram.ALN_METRIC.txt;else echo bam was converted to cram unsuccessfully 1>&2;exit;fi && \\\n";
        print CVSH echostring("$outdir/$samp/shell/Sentieon.bam2cram.$samp.sh");
        close CVSH;
    }
}

## Determining whether two BAM contain reads sequenced from the same sample or patient by counting genotype matches at common SNPs.
sub somaticQC{
    my ($normal,$tumor,$nbam,$tbam)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $QCDir="$outdir/somatic/$samplePair/QC";
	
    `mkdir -p $QCDir`;
    `mkdir -p $QCDir/tmp`;

    open SH, ">$QCDir/somaticQC.sh" or die $!;
    print SH "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    print SH "export PATH=$pythonBin:\$PATH && \\\n";
    print SH "export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:/opt/python/lib:\$LD_LIBRARY_PATH && \\\n";
    print SH "export PYTHONPATH=$bammathPATH/pythonLIB/lib/python2.7/site-packages && \\\n";
    print SH "$pythonBin/python $bammathPATH/bam-matcher.py --bam1 $nbam --bam2 $tbam --output $QCDir/bam_matcher.report.txt --vcf $bammathPATH/1kg.exome.highAF.merge.hg38.vcf.gz --caller gatk --dp-threshold 15 --gatk-mem-gb 4 --gatk-nt 1 --reference $reference --do-not-cache --experimental --scratch-dir $QCDir/tmp && \\\n";
    if($seqType =~ /target/i){
        print SH "$java8 -Xms4g -Xmx4g -Djava.io.tmpdir=$QCDir/tmp -XX:MaxMetaspaceSize=128m -XX:-UseGCOverheadLimit -XX:ParallelGCThreads=1 -jar $GATK3 -T ContEst -R $reference -I:eval $tbam -I:genotype $nbam --popfile $hapmapAF -L $hapmapAF -L $outdir/TR/$trname -isr INTERSECTION -o $QCDir/ContEst.report.txt && \\\n";
    }
    elsif($seqType =~ /wgs/i){
        print SH "$java8 -Xms4g -Xmx4g -Djava.io.tmpdir=$QCDir/tmp -XX:MaxMetaspaceSize=128m -XX:-UseGCOverheadLimit -XX:ParallelGCThreads=1 -jar $GATK3 -T ContEst -R $reference -I:eval $tbam -I:genotype $nbam --popfile $hapmapAF -L $hapmapAF -isr INTERSECTION -o $QCDir/ContEst.report.txt && \\\n";
    }
    print SH "rm -rf $QCDir/tmp && \\\n";
    print SH echostring("$QCDir/somaticQC.sh");
    close SH;
}

sub SnvIndelAnn{    
    #annovar 
    my ($soft,$samplePair,$type)=@_;
    open ANNO, ">$outdir/somatic/$samplePair/$soft/shell/$soft.$samplePair.ann.sh" || die $!;
    my $AnnotationCmd="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    die "$annovar does not exist!\n" if (!-e "$annovar");
    die "$annovardb does not exist!\n" if (!-e "$annovardb");
    my $normal=(split/-VS-/,$samplePair)[0];
    my $tumor=(split/-VS-/,$samplePair)[1];
    my @TP=split/,/,$type;
    my $table_annovar="$annovar/table_annovar.pl";
    if ($soft=~/Strelka2/){
        $table_annovar="$BIN/table_annovar_Strelka2.pl";
    }    
    foreach my $T(@TP){ 
        $AnnotationCmd.="/share/app/perl-5.22.0/bin/perl $table_annovar $outdir/somatic/$samplePair/$soft/result/$samplePair.$soft.$T.vcf.gz $annovardb -buildver hg38 -out $outdir/somatic/$samplePair/$soft/result/$samplePair.$soft.$T -remove -protocol refGene,ensGene,knownGene,cytoBand,genomicSuperDups,gwasCatalog,phastConsElements100way,wgRna,avsnp150,clinvar_20190305,cosmic89,dbnsfp35c,exac03nonpsych,exac03nontcga,gnomad_exome,gnomad_genome,intervar_20180118 -operation g,g,g,r,r,r,r,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput && \\\n";
        $AnnotationCmd.="gzip -f $outdir/somatic/$samplePair/$soft/result/$samplePair.$soft.$T.hg38_multianno.txt && \\\n";
        $AnnotationCmd.="rm -rf $outdir/somatic/$samplePair/$soft/result/$samplePair.$soft.$T.avinput $outdir/somatic/$samplePair/$soft/result/$samplePair.$soft.$T.hg38_multianno.vcf && \\\n";
    }
    if($soft=~/MuTect$/){
        $AnnotationCmd.="for i in {0..$MuTect_TRsplitNum};do rm -fr $outdir/somatic/$samplePair/MuTect/raw/$samplePair.\$i.MuTect.vcf.gz $outdir/somatic/$samplePair/MuTect/raw/$samplePair.\$i.MuTect.vcf.gz.idx $outdir/somatic/$samplePair/MuTect/raw/$samplePair.\$i.MuTect.txt;done && \\\n";
    }
    elsif($soft=~/SomaticSniper/){
        $AnnotationCmd.="for i in {0..$SomaticSniper_TRsplitNum};do rm -fr $outdir/somatic/$samplePair/SomaticSniper/raw/$normal.indel.\$i.pileup $outdir/somatic/$samplePair/SomaticSniper/raw/$tumor.indel.\$i.pileup;done && \\\n";
    }
    elsif($soft=~/MuTect2/){
        $AnnotationCmd.="for i in {0..$MuTect2_TRsplitNum};do rm -fr $outdir/somatic/$samplePair/MuTect2/raw/$samplePair.\$i.MuTect2.vcf.gz $outdir/somatic/$samplePair/MuTect2/raw/$samplePair.\$i.MuTect2.vcf.gz.tbi;done && \\\n";
    }
    elsif($soft=~/Strelka2/){
        $AnnotationCmd.="for i in {0..$Strelka2_TRsplitNum};do rm -fr $outdir/somatic/$samplePair/Strelka2/raw/RAW\$i;done && \\\n";
    }
    elsif($soft=~/MuSE/){
        $AnnotationCmd.="for i in {0..$MuSE_TRsplitNum};do rm -fr $outdir/somatic/$samplePair/MuSE/raw/$samplePair.\$i.MuSE.txt;done && \\\n";
    }
    elsif($soft=~/Platypus/){
        $AnnotationCmd.="for i in {0..$Platypus_TRsplitNum};do rm -fr $outdir/somatic/$samplePair/Platypus/raw/$samplePair.\$i.Platypus.vcf.gz;done && \\\n";
    }          
    $AnnotationCmd.=echostring("$outdir/somatic/$samplePair/$soft/shell/$soft.$samplePair.ann.sh");
    print ANNO $AnnotationCmd;
    close ANNO;
}

## calling somatic SNV & INDEL by Sentieon TNsnv (MuTect 1.1.5), TNhaplotyper (MuTect2, GATK 3.7/GATK3.8), TNhaplotyper2 (GATK 4.0.2.1 Mutect2)
sub TNcallerCalling{
    my ($normal,$tumor,$nbam,$tbam)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $TNcallerDir="$outdir/somatic/$samplePair/Sentieon";

    `mkdir -p $TNcallerDir`;
    `mkdir -p $TNcallerDir/shell`;
    `mkdir -p $TNcallerDir/raw`;
    `mkdir -p $TNcallerDir/result`;

    open TNSH, ">$TNcallerDir/shell/Sentieon.TNcaller.$samplePair.sh" || die $!;
    print TNSH "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    print TNSH "$SentieonAuto driver -t $SentieonCPU -r $reference";
    if($seqType =~ /target/i){ # WES & Panel
        print TNSH " --interval $outdir/TR/$trname";
    }
    print TNSH " -i $nbam -i $tbam --algo TNsnv --normal_sample $normal --tumor_sample $tumor --dbsnp $Dbsnp --cosmic $Cosmic --call_stats_out $TNcallerDir/raw/$samplePair.TNsnv.txt.gz $TNcallerDir/raw/$samplePair.TNsnv.vcf.gz --algo TNhaplotyper --normal_sample $normal --tumor_sample $tumor --dbsnp $Dbsnp --cosmic $Cosmic $TNcallerDir/raw/$samplePair.TNhaplotyper.vcf.gz --algo TNhaplotyper2 --normal_sample $normal --tumor_sample $tumor --germline_vcf $gnomad --default_af 3.125e-05 $TNcallerDir/raw/$samplePair.TNhaplotyper2.vcf.gz && \\\n";
    print TNSH "$SentieonAuto tnhapfilter --normal_sample $normal --tumor_sample $tumor -v $TNcallerDir/raw/$samplePair.TNhaplotyper2.vcf.gz $TNcallerDir/raw/$samplePair.TNhaplotyper2.Filter.vcf.gz && \\\n";
    print TNSH echostring("$TNcallerDir/shell/Sentieon.TNcaller.$samplePair.sh");
    close TNSH;
}

## calling somatic SNV by MuTect
sub MuTectCalling{
    my ($normal,$tumor,$nbam,$tbam)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $MuTectDir="$outdir/somatic/$samplePair/MuTect";
    my $Num=$MuTect_TRsplitNum-1;
    
    `mkdir -p $MuTectDir`;
    `mkdir -p $MuTectDir/shell`;
    `mkdir -p $MuTectDir/raw`;
    `mkdir -p $MuTectDir/result`;
    my $all_vcfs="";
    if ($seqType =~ /target/i && $MuTect_TRsplitNum eq 1){
        open MT,">$MuTectDir/shell/MuTect.$samplePair.sh" or die $!;
        my $MuTectCmd="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
        $MuTectCmd.="$java7 -Xms8g -Xmx8g -Djava.io.tmpdir=$MuTectDir/tmp -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -XX:ParallelGCThreads=1 -jar $MuTect -T MuTect -R $reference --input_file:normal $nbam --input_file:tumor $tbam --normal_sample_name $normal --tumor_sample_name $tumor --vcf $MuTectDir/raw/$samplePair.MuTect.vcf.gz --out $MuTectDir/raw/$samplePair.MuTect.txt -L $outdir/TR/$trname.intervals";
        $MuTectCmd.=" --dbsnp $Dbsnp" if(-f $Dbsnp);
        $MuTectCmd.=" --cosmic $Cosmic" if(-f $Cosmic);
        $MuTectCmd.=" --enable_extended_output --downsampling_type NONE && \\\n";
        $MuTectCmd.="rm -rf $MuTectDir/tmp && \\\n";
        $MuTectCmd.=echostring("$MuTectDir/shell/MuTect.$samplePair.sh");
        print MT $MuTectCmd;
        close MT;    
    }

    elsif(($seqType =~ /wgs/i) || ($seqType =~ /target/i && $MuTect_TRsplitNum > 1)){
        foreach my $n(0 .. $Num){
            my $bed="$outdir/TR/split_intervals$MuTect_TRsplitNum/$n-scattered.intervals.bed";     
            open MT,">$MuTectDir/shell/MuTect.$samplePair.$n.sh" or die $!;
            my $MuTectCmd="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
            $MuTectCmd.="$java7 -Xms8g -Xmx8g -Djava.io.tmpdir=$MuTectDir/tmp -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -XX:ParallelGCThreads=1 -jar $MuTect -T MuTect -R $reference --input_file:normal $nbam --input_file:tumor $tbam --normal_sample_name $normal --tumor_sample_name $tumor --vcf $MuTectDir/raw/$samplePair.$n.MuTect.vcf.gz -L $bed --out $MuTectDir/raw/$samplePair.$n.MuTect.txt";
            $MuTectCmd.=" --dbsnp $Dbsnp" if(-f $Dbsnp);
            $MuTectCmd.=" --cosmic $Cosmic" if(-f $Cosmic);
            $MuTectCmd.=" --enable_extended_output --downsampling_type NONE && \\\n";
            $MuTectCmd.=echostring("$MuTectDir/shell/MuTect.$samplePair.$n.sh");
            print MT $MuTectCmd;
            close MT;
            $all_vcfs.="I=$MuTectDir/raw/$samplePair.$n.MuTect.vcf.gz ";
        }   
    } 

    open MTPRO,">$MuTectDir/shell/MuTect.$samplePair.process.sh" or die $!;
    my $ProcessCmd="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n"; 
    $ProcessCmd.="export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH && \\\n";
    if(($seqType =~ /wgs/i) || ($seqType =~ /target/i && $MuTect_TRsplitNum > 1)){    
        $ProcessCmd.="$java8 -jar $Picard SortVcf SD=$refdict $all_vcfs O=$MuTectDir/raw/$samplePair.MuTect.vcf.gz && \\\n";
        $ProcessCmd.="$tabix -f -p vcf $MuTectDir/raw/$samplePair.MuTect.vcf.gz && \\\n";
        $ProcessCmd.="cat $MuTectDir/raw/$samplePair.0.MuTect.txt > $MuTectDir/raw/$samplePair.MuTect.txt && \\\n";
        $ProcessCmd.="for i in {1..$Num};do cat $MuTectDir/raw/$samplePair.\$i.MuTect.txt;done | grep -v  \"^#\" >> $MuTectDir/raw/$samplePair.MuTect.txt && \\\n";
    }
    $ProcessCmd.="gzip -f $MuTectDir/raw/$samplePair.MuTect.txt && \\\n";
    $ProcessCmd.="rm -rf $MuTectDir/tmp && \\\n";
 	
    #filter raw snv
    $ProcessCmd.="/share/app/perl-5.22.0/bin/perl $MuTectFilter $MuTectDir/raw/$samplePair.MuTect.txt.gz $MuTectDir/raw/$samplePair.MuTect.vcf.gz 0.02 0.05 $MuTectDir/result/$samplePair.MuTect.snv.txt.gz $MuTectDir/result/$samplePair.MuTect.snv.vcf.gz && \\\n";
    $ProcessCmd.=echostring("$MuTectDir/shell/MuTect.$samplePair.process.sh");
    print MTPRO $ProcessCmd;
    close MTPRO;
    #annovar
    &SnvIndelAnn("MuTect",$samplePair,"snv");
}

## calling somatic SNV by MuTect2
sub MuTect2Calling{
    my ($normal,$tumor,$nbam,$tbam)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $MuTect2Dir="$outdir/somatic/$samplePair/MuTect2";
    my $Num=$MuTect2_TRsplitNum-1;
    
    `mkdir -p $MuTect2Dir`;
    `mkdir -p $MuTect2Dir/shell`;
    `mkdir -p $MuTect2Dir/raw`;
    `mkdir -p $MuTect2Dir/result`;
    `mkdir -p $MuTect2Dir/tmp`;
   
    if ($seqType =~ /target/i && $MuTect2_TRsplitNum eq 1){
        open MT2,">$MuTect2Dir/shell/MuTect2.$samplePair.sh" or die $!;
        my $MuTect2Cmd="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
        $MuTect2Cmd.="export PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/JAVA/jre1.8.0_92/bin:\$PATH && \\\n";
        $MuTect2Cmd.="$GATK4 --java-options '-Xmx6g -XX:ParallelGCThreads=1 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$MuTect2Dir/tmp' Mutect2 -R $reference -I $tbam -tumor $tumor -I $nbam -normal $normal --germline-resource $gnomad --af-of-alleles-not-in-resource 3.125e-05 -L $outdir/TR/$trname.intervals -O $MuTect2Dir/raw/$samplePair.MuTect2.vcf.gz && \\\n";
        $MuTect2Cmd.=echostring("$MuTect2Dir/shell/MuTect2.$samplePair.sh");
        print MT2 $MuTect2Cmd;
        close MT2;
    }
    
    my $all_vcfs="";
    if(($seqType =~ /wgs/i) || ($seqType =~ /target/i && $MuTect2_TRsplitNum > 1)){        
        foreach my $n(0 .. $Num){
            my $bed="$outdir/TR/split_intervals$MuTect2_TRsplitNum/$n-scattered.intervals.bed"; 
            open MT2,">$MuTect2Dir/shell/MuTect2.$samplePair.$n.sh" or die $!;
            my $MuTect2Cmd="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
            $MuTect2Cmd.="export PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/JAVA/jre1.8.0_92/bin:\$PATH && \\\n"; 
            $MuTect2Cmd.="$GATK4 --java-options '-Xmx6g -XX:ParallelGCThreads=1 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$MuTect2Dir/tmp' Mutect2 -R $reference -I $tbam -tumor $tumor -I $nbam -normal $normal --germline-resource $gnomad --af-of-alleles-not-in-resource 3.125e-05 -L $bed -O $MuTect2Dir/raw/$samplePair.$n.MuTect2.vcf.gz && \\\n";
            $MuTect2Cmd.=echostring("$MuTect2Dir/shell/MuTect2.$samplePair.$n.sh");      
            print MT2 $MuTect2Cmd;
            close MT2;    
            $all_vcfs.="I=$MuTect2Dir/raw/$samplePair.$n.MuTect2.vcf.gz ";
        }
    }
   
    open MT2PRO,">$MuTect2Dir/shell/MuTect2.$samplePair.process.sh" or die $!;
    my $ProcessCmd="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    $ProcessCmd.="export PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/JAVA/jre1.8.0_92/bin:\$PATH && \\\n";
    $ProcessCmd.="export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH && \\\n";
    $ProcessCmd.="$java8 -jar $Picard SortVcf SD=$refdict $all_vcfs O=$MuTect2Dir/raw/$samplePair.MuTect2.vcf.gz && \\\n" if(($seqType =~ /wgs/i) || ($seqType =~ /target/i && $MuTect2_TRsplitNum > 1));
    $ProcessCmd.="$GATK4 --java-options '-Xmx6g -XX:ParallelGCThreads=1 -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$MuTect2Dir/tmp' FilterMutectCalls -V $MuTect2Dir/raw/$samplePair.MuTect2.vcf.gz -O $MuTect2Dir/raw/$samplePair.MuTect2.Filter.vcf.gz && \\\n";
    $ProcessCmd.="zcat $MuTect2Dir/raw/$samplePair.MuTect2.Filter.vcf.gz|awk '\$1~/^#/ || \$7~/PASS/ {print}' | $bgzip > $MuTect2Dir/result/$samplePair.MuTect2.snv_shortindel.vcf.gz && \\\n";
    $ProcessCmd.="$tabix -f -p vcf $MuTect2Dir/raw/$samplePair.MuTect2.vcf.gz && \\\n";
    $ProcessCmd.="$tabix -f -p vcf $MuTect2Dir/result/$samplePair.MuTect2.snv_shortindel.vcf.gz && \\\n";
    $ProcessCmd.=echostring("$MuTect2Dir/shell/MuTect2.$samplePair.process.sh");
    print MT2PRO $ProcessCmd;
    close MT2PRO;
             
    &SnvIndelAnn("MuTect2",$samplePair,"snv_shortindel");     
}

sub SomaticSniperCalling{
    my ($normal,$tumor,$nbam,$tbam)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $SomaticSniperDir="$outdir/somatic/$samplePair/SomaticSniper"; 
    my $Num=$SomaticSniper_TRsplitNum-1;
    `mkdir -p $SomaticSniperDir`;
    `mkdir -p $SomaticSniperDir/shell`;
    `mkdir -p $SomaticSniperDir/raw`;
    `mkdir -p $SomaticSniperDir/result`;
    
        open SomaticSniper, ">$SomaticSniperDir/shell/SomaticSniper.$samplePair.sh" || die $!;
        my $SspCmd="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
        $SspCmd.="$bam_somaticsniper/build/bin/bam-somaticsniper -q 1 -Q 40 -L -G -F vcf -t $tumor -n $normal -f $reference $tbam $nbam $SomaticSniperDir/raw/$samplePair.SomaticSniper.vcf.gz && \\\n";
        $SspCmd.=echostring("$SomaticSniperDir/shell/SomaticSniper.$samplePair.sh");
        print SomaticSniper $SspCmd;
        close SomaticSniper;
        
    if ($seqType =~ /target/i && $SomaticSniper_TRsplitNum eq 1){
        open Mpileup, ">$SomaticSniperDir/shell/SomaticSniper.$samplePair.mpileup.sh" || die $!;
        my $mpileupsh="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
        $mpileupsh.="export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH  && \\\n";
        $mpileupsh.="$bcftools mpileup -A -B -f $reference $nbam|$bcftools call -c|grep -v \"#\"|$vcfutils varFilter -Q 20|/share/app/perl-5.22.0/bin/perl $bam_somaticsniper/src/scripts/snpfilter.pl --snp-file $SomaticSniperDir/raw/$samplePair.SomaticSniper.vcf.gz --indel-file - --out-file $SomaticSniperDir/raw/$normal.SNPfilter && \\\n";
        $mpileupsh.="$bcftools mpileup -A -B -f $reference $tbam|$bcftools call -c|grep -v \"#\"|$vcfutils varFilter -Q 20|/share/app/perl-5.22.0/bin/perl $bam_somaticsniper/src/scripts/snpfilter.pl --snp-file $SomaticSniperDir/raw/$normal.SNPfilter --indel-file - --out-file $SomaticSniperDir/raw/$samplePair.SNPfilter && \\\n";
        $mpileupsh.=echostring("$SomaticSniperDir/shell/SomaticSniper.$samplePair.sh");
        print Mpileup $mpileupsh;
        close Mpileup;
    }

    if(($seqType =~ /wgs/i) || ($seqType =~ /target/i && $SomaticSniper_TRsplitNum > 1)){         
        foreach my $n(0 .. $Num){
            open Mpileup, ">$SomaticSniperDir/shell/SomaticSniper.$samplePair.$n.mpileup.sh" || die $!;
            my $bed="$outdir/TR/split_intervals$SomaticSniper_TRsplitNum/$n-scattered.intervals.bed"; 
            my $mpileupsh="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
            $mpileupsh.="export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH  && \\\n";
            $mpileupsh.="$bcftools mpileup -A -B -f $reference $nbam -R $bed|$bcftools call -c|grep -v \"#\"|$vcfutils varFilter -Q 20 > $SomaticSniperDir/raw/$normal.indel.$n.pileup && \\\n";
            $mpileupsh.="$bcftools mpileup -A -B -f $reference $tbam -R $bed|$bcftools call -c|grep -v \"#\"|$vcfutils varFilter -Q 20 > $SomaticSniperDir/raw/$tumor.indel.$n.pileup && \\\n";
            $mpileupsh.=echostring("$SomaticSniperDir/shell/SomaticSniper.$samplePair.$n.mpileup.sh");
            print Mpileup $mpileupsh;
            close Mpileup;
        }
    }
    open Filter,">$SomaticSniperDir/shell/SomaticSniper.$samplePair.filter.sh" || die $!; 
    my $filtersh="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    if(($seqType =~ /wgs/i) || ($seqType =~ /target/i && $SomaticSniper_TRsplitNum > 1)){ 
        $filtersh.="for i in {0..$Num};do cat $SomaticSniperDir/raw/$normal.indel.\$i.pileup >> $SomaticSniperDir/raw/$normal.indel.pileup;done && \\\n";
        $filtersh.="for i in {0..$Num};do cat $SomaticSniperDir/raw/$tumor.indel.\$i.pileup >> $SomaticSniperDir/raw/$tumor.indel.pileup;done && \\\n";
    
        $filtersh.="/share/app/perl-5.22.0/bin/perl $bam_somaticsniper/src/scripts/snpfilter.pl --snp-file $SomaticSniperDir/raw/$samplePair.SomaticSniper.vcf.gz --indel-file $SomaticSniperDir/raw/$normal.indel.pileup --out-file $SomaticSniperDir/raw/$normal.SNPfilter && \\\n";
        $filtersh.="/share/app/perl-5.22.0/bin/perl $bam_somaticsniper/src/scripts/snpfilter.pl --snp-file $SomaticSniperDir/raw/$normal.SNPfilter --indel-file $SomaticSniperDir/raw/$tumor.indel.pileup --out-file $SomaticSniperDir/raw/$samplePair.SNPfilter && \\\n";
    }

    $filtersh.="/share/app/perl-5.22.0/bin/perl $bam_somaticsniper/src/scripts/prepare_for_readcount.pl --snp-file $SomaticSniperDir/raw/$samplePair.SNPfilter --out-file $SomaticSniperDir/raw/$samplePair.SNPfilter.pos && \\\n"; 
    $filtersh.="export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/bam-readcount/v0.8.0/lib:/share/app/gcc-4.9.3/lib64:\$LD_LIBRARY_PATH && \\\n";
    $filtersh.="$bam_readcount -w 1 -b 15 -q 1 -f $reference -l $SomaticSniperDir/raw/$samplePair.SNPfilter.pos $tbam > $SomaticSniperDir/raw/$samplePair.readcounts && \\\n";
    $filtersh.="/share/app/perl-5.22.0/bin/perl $bam_somaticsniper/src/scripts/fpfilter.pl --snp-file $SomaticSniperDir/raw/$samplePair.SNPfilter --readcount-file $SomaticSniperDir/raw/$samplePair.readcounts && \\\n";
    ##--out-file, snp output file after filter
    if($seqType =~ /wgs/i){
        $filtersh.="/share/app/perl-5.22.0/bin/perl $bam_somaticsniper/src/scripts/highconfidence.pl --min-mapping-quality 40 --min-somatic-score 40 --snp-file $SomaticSniperDir/raw/$samplePair.SNPfilter.fp_pass --lq-output $SomaticSniperDir/result/$samplePair.low_qual.vcf.gz --out-file $SomaticSniperDir/result/$samplePair.SomaticSniper.snv.vcf && \\\n";
    }
    elsif($seqType=~/target/i){
        $filtersh.="export LD_LIBRARY_PATH=/share/app/gcc-4.9.3/lib64:/share/app/libz/zlib-1.2.11:\$LD_LIBRARY_PATH && \\\n";
        $filtersh.="/share/app/perl-5.22.0/bin/perl $bam_somaticsniper/src/scripts/highconfidence.pl --min-mapping-quality 40 --min-somatic-score 40 --snp-file $SomaticSniperDir/raw/$samplePair.SNPfilter.fp_pass --lq-output $SomaticSniperDir/result/$samplePair.low_qual.vcf.gz --out-file $SomaticSniperDir/raw/$samplePair.SomaticSniper.highconf.snv.vcf && \\\n";
        $filtersh.="$bedtools intersect -header -a $SomaticSniperDir/raw/$samplePair.SomaticSniper.highconf.snv.vcf -b $targetRegion > $SomaticSniperDir/result/$samplePair.SomaticSniper.snv.vcf && \\\n";
    }
    $filtersh.="export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH && \\\n";
    $filtersh.="$bgzip -f $SomaticSniperDir/result/$samplePair.SomaticSniper.snv.vcf  && \\\n";
    $filtersh.="$tabix -f -p vcf $SomaticSniperDir/result/$samplePair.SomaticSniper.snv.vcf.gz  && \\\n"; 
    $filtersh.=echostring("$SomaticSniperDir/shell/SomaticSniper.$samplePair.filter.sh");
    print Filter $filtersh; 
    close Filter; 
    #annovar 
    &SnvIndelAnn("SomaticSniper",$samplePair,"snv");   
}

## calling somatic snv by MuSE
sub MuSECalling{
    my ($normal,$tumor,$nbam,$tbam)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $MuSEDir="$outdir/somatic/$samplePair/MuSE";
    my $Num=$MuSE_TRsplitNum-1;
    
    `mkdir -p $MuSEDir`;
    `mkdir -p $MuSEDir/shell`;
    `mkdir -p $MuSEDir/raw`;
    `mkdir -p $MuSEDir/result`;
	  
    if ($seqType =~ /target/i && $MuSE_TRsplitNum eq 1){
        open SHELL, ">$MuSEDir/shell/MuSE.$samplePair.sh" || die $!;
        print SHELL "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
        print SHELL "export LD_LIBRARY_PATH=/share/app/gcc-4.9.3/lib64:/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH && \\\n";
        print SHELL "$MuSE call -f $reference -l $outdir/TR/$trname -O $MuSEDir/raw/$samplePair $tbam $nbam && \\\n";
        print SHELL "colnum=`cat $MuSEDir/raw/$samplePair.MuSE.txt|tail -1|awk -F \"\\t\" '{print NF}'`\n";
        print SHELL "$MuSE sump -I $MuSEDir/raw/$samplePair.MuSE.txt -E -O $MuSEDir/raw/$samplePair.MuSE.vcf -D $Dbsnp && \\\n";
        print SHELL echostring("$MuSEDir/shell/MuSE.$samplePair.sh");
        print SHELL "else echo MuSE output is not completed 1>&2;fi\n";
        close SHELL;
    }

    if(($seqType =~ /wgs/i) || ($seqType =~ /target/i && $MuSE_TRsplitNum > 1)){
        foreach my $n(0 .. $Num){
            my $bed="$outdir/TR/split_intervals$MuSE_TRsplitNum/$n-scattered.intervals.bed"; 
            open SHELL,">$MuSEDir/shell/MuSE.$samplePair.$n.sh" or die $!;
            print SHELL "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
            print SHELL "export LD_LIBRARY_PATH=/share/app/gcc-4.9.3/lib64:/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH && \\\n";
            print SHELL "$MuSE call -f $reference -l $bed -O $MuSEDir/raw/$samplePair.$n $tbam $nbam && \\\n";
            print SHELL "colnum=`cat $MuSEDir/raw/$samplePair.$n.MuSE.txt|tail -1|awk -F \"\\t\" '{print NF}'`\n";
            print SHELL "if [ \$colnum == 30 ];then\n";
            print SHELL echostring("$MuSEDir/shell/MuSE.$samplePair.$n.sh");
            print SHELL "else echo MuSE output is not completed 1>&2;fi\n";
            close SHELL;
        }
    }	

    #filter raw snv.  From experience, we suggest using calls up to Tier 4 for WES data, and calls up to Tier 5 for WGS data. 
    open MSPRO,">$MuSEDir/shell/MuSE.$samplePair.process.sh" or die $!;
    my $ProcessCmd="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    $ProcessCmd.="export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH && \\\n";
    if(($seqType =~ /wgs/i) || ($seqType =~ /target/i && $MuSE_TRsplitNum > 1)){
        $ProcessCmd.="cat $MuSEDir/raw/$samplePair.0.MuSE.txt > $MuSEDir/raw/$samplePair.MuSE.txt && \\\n"; 
        $ProcessCmd.="for i in {1..$Num};do cat $MuSEDir/raw/$samplePair.\$i.MuSE.txt;done | grep -v  \"^#\" >> $MuSEDir/raw/$samplePair.MuSE.txt && \\\n";
        $ProcessCmd.="$MuSE sump -I $MuSEDir/raw/$samplePair.MuSE.txt -O $MuSEDir/raw/$samplePair.MuSE.vcf -D $Dbsnp";
        $ProcessCmd.=" -G && \\\n" if($seqType =~ /wgs/i);
        $ProcessCmd.=" -E && \\\n" if($seqType =~ /target/i);
    }
    $ProcessCmd.="$bgzip -f $MuSEDir/raw/$samplePair.MuSE.vcf && \\\n";
    $ProcessCmd.="$bgzip -f $MuSEDir/raw/$samplePair.MuSE.txt && \\\n";
    $ProcessCmd.="$tabix -f -p vcf $MuSEDir/raw/$samplePair.MuSE.vcf.gz && \\\n";
    $ProcessCmd.="zcat $MuSEDir/raw/$samplePair.MuSE.vcf.gz | awk '\$1~/^#/ || \$7~/PASS|Tier1|Tier2|Tier3|Tier4/ {print}' | $bgzip > $MuSEDir/result/$samplePair.MuSE.snv.vcf.gz && \\\n" if($seqType =~ /target/i); 
    $ProcessCmd.="zcat $MuSEDir/raw/$samplePair.MuSE.vcf.gz | awk '\$1~/^#/ || \$7~/PASS|Tier1|Tier2|Tier3|Tier4|Tier5/ {print}' | $bgzip > $MuSEDir/result/$samplePair.MuSE.snv.vcf.gz && \\\n" if($seqType =~ /wgs/i); 
    $ProcessCmd.="$tabix -f -p vcf $MuSEDir/result/$samplePair.MuSE.snv.vcf.gz && \\\n";
    $ProcessCmd.=echostring("$MuSEDir/shell/MuSE.$samplePair.process.sh");
    print MSPRO $ProcessCmd;
    close MSPRO;
    ##annovar
    &SnvIndelAnn("MuSE",$samplePair,"snv");
}

sub MantaCalling{
    my($normal,$tumor,$nbam,$tbam)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $MantaDir="$outdir/somatic/$samplePair/Manta";
    `mkdir -p $MantaDir/shell`;
    `mkdir -p $MantaDir/result`;
    #calling structural variants (SVs) and indel
    open MANTA,">$MantaDir/shell/Manta.$samplePair.sh";
    my $MantaCmd="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    $MantaCmd.="export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:/opt/python/lib:\$LD_LIBRARY_PATH && \\\n";
    $MantaCmd.="rm -rf $MantaDir/result/* && \\\n";
    if($seqType =~ /target/i){
        $MantaCmd.="$pythonBin/python $Manta --normalBam $nbam --tumorBam $tbam --referenceFasta $reference --runDir $MantaDir/result --exome --callRegions $outdir/TR/$trname.gz --generateEvidenceBam --outputContig && \\\n";
    }
    elsif($seqType =~ /wgs/i){
        $MantaCmd.="$pythonBin/python $Manta --normalBam $nbam --tumorBam $tbam --referenceFasta $reference --runDir $MantaDir/result --callRegions $outdir/TR/$trname.gz --generateEvidenceBam --outputContig && \\\n";
    }
    $MantaCmd.="$pythonBin/python $MantaDir/result/runWorkflow.py -m local -j 20 && \\\n";
	$MantaCmd.="rm -rf $MantaDir/result/workspace && \\\n";
    print MANTA $MantaCmd;
    print MANTA echostring("$MantaDir/shell/Manta.$samplePair.sh");
    close MANTA;
    
    open MANTAANN,">$MantaDir/shell/Manta.$samplePair.ann.sh" || die $!;
    my $AnnotationCmd="#!/bin/bash\necho ==========start at : `date` ==========\n";
    $AnnotationCmd.="check=`zcat $MantaDir/result/results/variants/somaticSV.vcf.gz|grep -v \"#\"` \n";
    $AnnotationCmd.="if [ \"\$check\" ];then export ANNOTSV=$BIN/AnnotSV/v2.0/ && \\\n";
    $AnnotationCmd.="export PATH=/hwfssz1/ST_CANCER/POL/USER/lifuqiang/tools/anaconda2/bin:\$PATH  && \\\n";
    $AnnotationCmd.="$AnnotSV -bedtools $bedtools_2_26_0 -genomeBuild GRCh38 -SVinputFile $MantaDir/result/results/variants/somaticSV.vcf.gz -SVinputInfo 1 -svtBEDcol 4 -outputDir $MantaDir/result/results/variants/ -outputFile $MantaDir/result/results/variants/$samplePair.Manta.sv.hg38_multianno.tsv >& $MantaDir/shell/Manta.$samplePair.ann.sh.log && \\\n";
    $AnnotationCmd.="export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH && \\\n";
    $AnnotationCmd.="$bgzip -f $MantaDir/result/results/variants/$samplePair.Manta.sv.hg38_multianno.tsv \n";
    $AnnotationCmd.="else echo no SV detected > $MantaDir/shell/Manta.$samplePair.ann.sh.log;fi && \\\n";
    $AnnotationCmd.=echostring("$MantaDir/shell/Manta.$samplePair.ann.sh");
    print MANTAANN $AnnotationCmd;
    close MANTAANN;    
    
}

sub Strelka2Calling{
    my ($normal,$tumor,$nbam,$tbam)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $Strelka2Dir="$outdir/somatic/$samplePair/Strelka2";
    my $Num=$Strelka2_TRsplitNum-1;
    `mkdir -p $Strelka2Dir/shell`;
    #`mkdir -p $Strelka2Dir/raw`;
    `mkdir -p $Strelka2Dir/result`;
    #calling snv/indel
    if ($seqType =~ /target/i && $Strelka2_TRsplitNum eq 1){
        `mkdir -p $Strelka2Dir/raw`;            
        open CFGSH,">$Strelka2Dir/shell/Strelka2.$samplePair.sh";
        my $ST2CONFIG="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
        $ST2CONFIG.="export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:/opt/python/lib:\$LD_LIBRARY_PATH && \\\n";
        $ST2CONFIG.="rm -fr $Strelka2Dir/raw/* && \\\n";
        my $add_para="";
        if ($Use_MantaAndStrelka2=~/true/i){
            $add_para="--indelCandidates $outdir/somatic/$samplePair/Manta/result/results/variants/candidateSmallIndels.vcf.gz";
        }  
        $ST2CONFIG.="$pythonBin/python $Strelka2 --normalBam $nbam --tumorBam $tbam --ref $reference --exome --callRegions $outdir/TR/$trname.gz  $add_para --runDir $Strelka2Dir/raw && \\\n";
        $ST2CONFIG.="$pythonBin/python $Strelka2Dir/raw/runWorkflow.py -m local -j 4 && \\\n"; 
        $ST2CONFIG.="export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH && \\\n";       
        $ST2CONFIG.="zcat $Strelka2Dir/raw/results/variants/somatic.snvs.vcf.gz|$bgzip > $Strelka2Dir/raw/$samplePair.Strelka2.snv.vcf.gz && \\\n";
        $ST2CONFIG.="zcat $Strelka2Dir/raw/results/variants/somatic.indels.vcf.gz|$bgzip > $Strelka2Dir/raw/$samplePair.Strelka2.indel.vcf.gz && \\\n";
        print CFGSH $ST2CONFIG;
        print CFGSH echostring("$Strelka2Dir/shell/Strelka2.$samplePair.sh");
        close CFGSH;
    }
    my $all_vcfs="";
    my $all_indels="";  
    my $add_para="";
    if(($seqType =~ /wgs/i) || ($seqType =~ /target/i && $Strelka2_TRsplitNum > 1)){
        foreach my $n(0 .. $Num){
            my $bed="$outdir/TR/split_intervals$Strelka2_TRsplitNum/$n-scattered.intervals.bed.gz"; 
            `mkdir -p $Strelka2Dir/raw/RAW$n`;
            open CFGSH,">$outdir/somatic/$samplePair/Strelka2/shell/Strelka2.$samplePair.$n.sh";
            my $ST2CONFIG="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
            $ST2CONFIG.="export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:/opt/python/lib:\$LD_LIBRARY_PATH && \\\n";
            $ST2CONFIG.="rm -rf $Strelka2Dir/raw/RAW$n/* && \\\n";
            if ($Use_MantaAndStrelka2=~/true/i){
                $add_para="--indelCandidates $outdir/somatic/$samplePair/Manta/result/results/variants/candidateSmallIndels.vcf.gz";
            }
            if($seqType =~ /target/i){
                $ST2CONFIG.="$pythonBin/python $Strelka2 --normalBam $nbam --tumorBam $tbam --ref $reference --exome --callRegions $bed $add_para --runDir $Strelka2Dir/raw/RAW$n && \\\n";
            }
            elsif($seqType =~ /wgs/i){
                $ST2CONFIG.="$pythonBin/python $Strelka2 --normalBam $nbam --tumorBam $tbam --ref $reference --callRegions $bed $add_para --runDir $Strelka2Dir/raw/RAW$n && \\\n";
            }           
            $ST2CONFIG.="$pythonBin/python $Strelka2Dir/raw/RAW$n/runWorkflow.py -m local -j 4 && \\\n";
            print CFGSH $ST2CONFIG;
            print CFGSH echostring("$Strelka2Dir/shell/Strelka2.$samplePair.$n.sh");
            $all_vcfs.="I=$Strelka2Dir/raw/RAW$n/results/variants/somatic.snvs.vcf.gz ";  
            $all_indels.="I=$Strelka2Dir/raw/RAW$n/results/variants/somatic.indels.vcf.gz ";
        }
    }
    
    open ST2PRO,">$Strelka2Dir/shell/Strelka2.$samplePair.process.sh" or die $!;
    my $ProcessCmd="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    if(($seqType =~ /wgs/i) || ($seqType =~ /target/i && $Strelka2_TRsplitNum > 1)){
        $ProcessCmd.="$java8 -jar $Picard SortVcf SD=$refdict $all_vcfs O=$Strelka2Dir/raw/$samplePair.Strelka2.snv.vcf.gz && \\\n";
        $ProcessCmd.="$java8 -jar $Picard SortVcf SD=$refdict $all_indels O=$Strelka2Dir/raw/$samplePair.Strelka2.indel.vcf.gz && \\\n";
    }
    $ProcessCmd.="export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH && \\\n";                
    $ProcessCmd.="zcat $Strelka2Dir/raw/$samplePair.Strelka2.snv.vcf.gz | awk '\$1~/^#/ || \$7~/PASS/ {print}' | $bgzip > $Strelka2Dir/result/$samplePair.Strelka2.snv.vcf.gz && \\\n"; 
    $ProcessCmd.="zcat $Strelka2Dir/raw/$samplePair.Strelka2.indel.vcf.gz | awk '\$1~/^#/ || \$7~/PASS/ {print}' | $bgzip > $Strelka2Dir/result/$samplePair.Strelka2.indel.vcf.gz && \\\n"; 
    print ST2PRO $ProcessCmd;
    print ST2PRO echostring("$Strelka2Dir/shell/Strelka2.$samplePair.process.sh");
    close ST2PRO;    
    #annovar 
    &SnvIndelAnn("Strelka2",$samplePair,"snv,indel");
}

## calling somatic snv/indel by Lancet
sub LancetCalling{
    my ($normal,$tumor,$nbam,$tbam)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $LancetDir="$outdir/somatic/$samplePair/Lancet";
    `mkdir -p $LancetDir/shell/`;
    `mkdir -p $LancetDir/raw`;
    `mkdir -p $LancetDir/result`;
    ##calling snv/indel
    open LANCET, ">$LancetDir/shell/Lancet.$samplePair.sh";
    my $LancetCmd="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    $LancetCmd.="export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH && \\\n";
    $LancetCmd.="$Lancet --tumor $tbam --normal $nbam --ref $reference --bed $outdir/somatic/$samplePair/called/$samplePair.targetForLancet.bed --num-thread 4|$bgzip > $LancetDir/raw/Lancet.$samplePair.raw.vcf.gz && \\\n";
    $LancetCmd.="zcat $LancetDir/raw/Lancet.$samplePair.raw.vcf.gz|awk '\$1~/^#/ || \$7~/PASS/ {print}' | $bgzip > $LancetDir/result/$samplePair.Lancet.snv_shortindel.vcf.gz && \\\n";
    print LANCET $LancetCmd;
    print LANCET echostring("$LancetDir/shell/Lancet.$samplePair.sh");
    close LANCET;        
    
    #annovar 
    &SnvIndelAnn("Lancet",$samplePair,"snv_shortindel");
}

## calling somatic indel by Platypus
sub PlatypusCalling{
    my ($normal,$tumor,$nbam,$tbam)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $PlatypusDir="$outdir/somatic/$samplePair/Platypus";
    my $Num=$Platypus_TRsplitNum-1;
    `mkdir -p $PlatypusDir`;
    `mkdir -p $PlatypusDir/shell`;
    `mkdir -p $PlatypusDir/raw`;
    `mkdir -p $PlatypusDir/result`;

    #calling sindel
    if ($seqType =~ /target/i && $Platypus_TRsplitNum eq 1){
        open SINDELSH, ">$PlatypusDir/shell/Platypus.$samplePair.sh";
        print SINDELSH "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
        print SINDELSH "export PATH=$pythonBin:\$PATH\n";
        print SINDELSH "export LD_LIBRARY_PATH=$ldlib:\$LD_LIBRARY_PATH\n";
        print SINDELSH "export PYTHONPATH=$pythonPath:\${PYTHONPATH-}\n";
        print SINDELSH "$pythonBin/python $Platypus/Platypus.py callVariants $PlatypusPara --refFile=$reference --bamFiles=$nbam,$tbam --regions=$outdir/TR/$trname --output=$PlatypusDir/raw/$samplePair.Platypus.vcf && \\\n";
        print SINDELSH echostring("$PlatypusDir/shell/Platypus.$samplePair.sh");
        close SINDELSH;
    }
    
    my $all_vcfs;
    if(($seqType =~ /wgs/i) || ($seqType =~ /target/i && $Platypus_TRsplitNum > 1)){
        foreach my $n(0 .. $Num){
            open SINDELSH, ">$PlatypusDir/shell/Platypus.$samplePair.$n.sh" or die $!;
            my $bed="$outdir/TR/split_intervals$Platypus_TRsplitNum/$n-scattered.intervals.bed"; 
            print SINDELSH "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
            print SINDELSH "export PATH=$pythonBin:\$PATH\n";
            print SINDELSH "export LD_LIBRARY_PATH=$ldlib:\$LD_LIBRARY_PATH\n";
            print SINDELSH "export PYTHONPATH=$pythonPath:\${PYTHONPATH-}\n";
            print SINDELSH "$pythonBin/python $Platypus/Platypus.py callVariants $PlatypusPara --refFile=$reference --bamFiles=$nbam,$tbam --regions=$bed --output=$PlatypusDir/raw/$samplePair.$n.Platypus.vcf && \\\n";
            print SINDELSH "export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH && \\\n";
            print SINDELSH "$bgzip -f $PlatypusDir/raw/$samplePair.$n.Platypus.vcf && \\\n";
            print SINDELSH echostring("$PlatypusDir/shell/Platypus.$samplePair.$n.sh"); 
            close SINDELSH;            
            $all_vcfs.="I=$PlatypusDir/raw/$samplePair.$n.Platypus.vcf.gz ";
        }            
    }
    
    open PLTPRO, ">$PlatypusDir/shell/Platypus.$samplePair.process.sh" or die $!;
    print PLTPRO "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    print PLTPRO "export PATH=$pythonBin:\$PATH\n";
    print PLTPRO "export LD_LIBRARY_PATH=$ldlib:\$LD_LIBRARY_PATH\n";
    print PLTPRO "export PYTHONPATH=$pythonPath:\${PYTHONPATH-}\n";
    if(($seqType =~ /wgs/i) || ($seqType =~ /target/i && $Platypus_TRsplitNum > 1)){
        print PLTPRO "$java8 -jar $Picard SortVcf SD=$refdict $all_vcfs O=$PlatypusDir/raw/$samplePair.Platypus.vcf && \\\n";  
    }   
	  print PLTPRO "$pythonBin/python $Platypus/somaticMutationDetector.py $PlatypussmtStrictPara --inputVCF $PlatypusDir/raw/$samplePair.Platypus.vcf --outputVCF $PlatypusDir/raw/$samplePair.Platypus.somatic.vcf --tumourSample $tumor --normalSample $normal && \\\n";
    print PLTPRO "export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH && \\\n";
    print PLTPRO "$bgzip -f $PlatypusDir/raw/$samplePair.Platypus.somatic.vcf && \\\n";
    print PLTPRO "$bgzip -f $PlatypusDir/raw/$samplePair.Platypus.vcf && \\\n";
    print PLTPRO "$tabix -f -p vcf $PlatypusDir/raw/$samplePair.Platypus.somatic.vcf.gz && \\\n";
    print PLTPRO "/share/app/perl-5.22.0/bin/perl $PlatypusSomaticIndelFilter $PlatypusDir/raw/$samplePair.Platypus.somatic.vcf.gz INDEL $PlatypusDir/result/$samplePair.Platypus.indel.vcf.gz $normal 0.02 $tumor 0.001 && \\\n";  
    print PLTPRO echostring("$PlatypusDir/shell/Platypus.$samplePair.process.sh");
    close PLTPRO;

    #annovar
    &SnvIndelAnn("Platypus",$samplePair,"indel");
}

## calling somatic CNV by FACETS
sub RunFACETS{
    my ($normal,$tumor,$nbam,$tbam)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $FACETSDir="$outdir/somatic/$samplePair/FACETS";
    
    `mkdir -p $FACETSDir`;
    `mkdir -p $FACETSDir/shell`;
    `mkdir -p $FACETSDir/result`;
	
    open FSH, ">$FACETSDir/shell/FACETS.$samplePair.pileup.sh" || die $!;
    print FSH "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";	
    print FSH "export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH && \\\n";
    print FSH "if [ -e \"$FACETSDir/result/$samplePair.FACETS.pileup.gz\" ];then rm -rf $FACETSDir/result/$samplePair.FACETS.pileup.gz;fi && \\\n";
    if($seqType =~ /target/i){
        print FSH "$BIN/snp-pileup -g -q15 -Q20 -P100 -r25,0 $Dbsnpcommon $FACETSDir/result/$samplePair.FACETS.pileup.gz $nbam $tbam && \\\n";
    }
    elsif($seqType =~ /wgs/i){
        print FSH "$BIN/snp-pileup -g -q15 -Q20 -P100 -r10,0 $Dbsnpcommon $FACETSDir/result/$samplePair.FACETS.pileup.gz $nbam $tbam && \\\n";
    }
    print FSH echostring("$FACETSDir/shell/FACETS.$samplePair.pileup.sh");
    close FSH;
	
    open FSH, ">$FACETSDir/shell/FACETS.$samplePair.calling.sh" || die $!;
    print FSH "#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    print FSH "export R_LIBS=/hwfssz1/ST_CANCER/POL/USER/lifuqiang/tools/library/R-3.3.2:/share/app/R-3.3.2/lib64/R/library && \\\n";
    print FSH "export LD_LIBRARY_PATH=/share/app/gcc-4.9.3/lib64:\$LD_LIBRARY_PATH && \\\n";
    print FSH "/share/app/R-3.3.2/bin/Rscript $BIN/facets_hg38.R 9922 $FACETSDir/result/$samplePair.FACETS.pileup.gz $FACETSDir/result/$samplePair.FACETS.purity.txt $FACETSDir/result/$samplePair.FACETS.cnv.txt $FACETSDir/result/$samplePair.FACETS.cnv.pdf && \\\n";
    print FSH "gzip -f $FACETSDir/result/$samplePair.FACETS.cnv.txt && \\\n";
    print FSH echostring("$FACETSDir/shell/FACETS.$samplePair.calling.sh");
    close FSH;
	
    #Annotate CNV by AnnotSV
    open FACETSANN,">$FACETSDir/shell/FACETS.$samplePair.ann.sh" || die $!;
    my $AnnotationCmd="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    $AnnotationCmd.="zcat $FACETSDir/result/$samplePair.FACETS.cnv.txt.gz|awk '{if (\$1 ==\"chrom\") print \$1\"\\t\"\$10\"\\t\"\$11\"\\t\"\$3\"\\t\"\$5\"\\t\"\$12\"\\t\"\$13\"\\t\"\$14;else if (\$1!=\"chrom\") print \"chr\"\$1\"\\t\"\$10\"\\t\"\$11\"\\t\"\$3\"\\t\"\$5\"\\t\"\$12\"\\t\"\$13\"\\t\"\$14}' > $FACETSDir/result/$samplePair.FACETS.cnv.bed && \\\n";
    $AnnotationCmd.="if [ -s $FACETSDir/result/$samplePair.FACETS.cnv.bed ];then export ANNOTSV=$BIN/AnnotSV/v2.0/ && \\\n";
    $AnnotationCmd.="export PATH=/hwfssz1/ST_CANCER/POL/USER/lifuqiang/tools/anaconda2/bin:\$PATH  && \\\n";
    $AnnotationCmd.="$AnnotSV -bedtools $bedtools_2_26_0 -genomeBuild GRCh38 -SVinputFile $FACETSDir/result/$samplePair.FACETS.cnv.bed -SVinputInfo 1 -svtBEDcol 4 -outputDir $FACETSDir/result -outputFile $FACETSDir/result/$samplePair.FACETS.cnv.hg38_multianno.tsv >& $FACETSDir/shell/FACETS.$samplePair.ann.sh.log && \\\n";
    $AnnotationCmd.="export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH && \\\n";
    $AnnotationCmd.="$bgzip -f $FACETSDir/result/$samplePair.FACETS.cnv.hg38_multianno.tsv \n";
    $AnnotationCmd.="else echo no CNA detected > $FACETSDir/shell/FACETS.$samplePair.ann.sh.log;fi && \\\n";
    $AnnotationCmd.=echostring("$FACETSDir/shell/FACETS.$samplePair.ann.sh");
    print FACETSANN $AnnotationCmd;
    close FACETSANN; 
}

sub runSvaba{
    my ($normal,$tumor,$nbam,$tbam)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $SvabaDir="$outdir/somatic/$samplePair/Svaba";
    `mkdir -p $SvabaDir`;
    `mkdir -p $SvabaDir/shell`;
    `mkdir -p $SvabaDir/result`;
    open SvABA,">$SvabaDir/shell/Svaba.$samplePair.sh" || die $!;
    my $SvabaCmd="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    $SvabaCmd.="$Svaba run -z -t $tbam -n $nbam -D $Dbsnp -p 20 -a $SvabaDir/result/$samplePair -G $reference && \\\n";
    print SvABA $SvabaCmd;
    print SvABA echostring("$SvabaDir/shell/Svaba.$samplePair.sh");
    close SvABA;
    
    open SVABAANN,">$SvabaDir/shell/Svaba.$samplePair.ann.sh" || die $!;
    my $AnnotationCmd="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    $AnnotationCmd.="check=`zcat $SvabaDir/result/$samplePair.svaba.somatic.sv.vcf.gz|grep -v \"#\"` \n";
    $AnnotationCmd.="if [ \"\$check\" ];then export ANNOTSV=$BIN/AnnotSV/v2.0/ && \\\n";
    $AnnotationCmd.="export PATH=/hwfssz1/ST_CANCER/POL/USER/lifuqiang/tools/anaconda2/bin:\$PATH  && \\\n";
    $AnnotationCmd.="$AnnotSV -bedtools $bedtools_2_26_0 -genomeBuild GRCh38 -SVinputFile $SvabaDir/result/$samplePair.svaba.somatic.sv.vcf.gz -SVinputInfo 1 -svtBEDcol 4 -outputDir $SvabaDir/result -outputFile $SvabaDir/result/$samplePair.Svaba.sv.hg38_multianno.tsv >& $SvabaDir/shell/Svaba.$samplePair.ann.sh.log && \\\n";
    $AnnotationCmd.="export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH && \\\n";
    $AnnotationCmd.="$bgzip -f $SvabaDir/result/$samplePair.Svaba.sv.hg38_multianno.tsv \n";
    $AnnotationCmd.="else echo no SV detected > $SvabaDir/shell/Svaba.$samplePair.ann.sh.log;fi && \\\n";    
    $AnnotationCmd.=echostring("$SvabaDir/shell/Svaba.$samplePair.ann.sh");
    print SVABAANN $AnnotationCmd;
    close SVABAANN;    
}

#&mergeBed;
sub mergeBed{
    my ($normal,$tumor,$nbam,$tbam,$soft,$type)=@_;
    my $samplePair="$normal-VS-$tumor";
    `mkdir -p $outdir/somatic/$samplePair/called`;
    `mkdir -p $outdir/somatic/$samplePair/called/shell`;
    my @TP=split/,/,$type;
    open(GBSH,">$outdir/somatic/$samplePair/$soft/shell/$soft.$samplePair.getbed.sh") or die "$!\n";
    my $GetBedSH="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    $GetBedSH.="rm -fr $outdir/somatic/$samplePair/$soft/result/$samplePair.$soft.called.bed && \\\n";
    $GetBedSH.="export LD_LIBRARY_PATH=/share/app/gcc-4.9.3/lib64:/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH && \\\n";
    foreach my $T(@TP){     
        $GetBedSH.="$Vt decompose -s $outdir/somatic/$samplePair/$soft/result/$samplePair.$soft.$T.vcf.gz | $bcftools norm -f $reference| $bcftools query -f '\%CHROM\\t\%POS\\t\%END\\t%REF\\t%ALT\\n'|awk '{if (length(\$4)==1 && length(\$5)==1) print \$1\"\\t\"\$2\"\\t\"\$3;else print \$1\"\\t\"\$2-10\"\\t\"\$3+10 }' >> $outdir/somatic/$samplePair/$soft/result/$samplePair.$soft.called.bed && \\\n"; 
    }       
        ###
    $GetBedSH.="cat $outdir/somatic/$samplePair/$soft/result/$samplePair.$soft.called.bed >> $outdir/somatic/$samplePair/called/$samplePair.NonLancet.called.bed && \\\n" if ($soft !~ /Lancet/);        
    print GBSH $GetBedSH;
    print GBSH echostring("$outdir/somatic/$samplePair/$soft/shell/$soft.$samplePair.getbed.sh"); 
    ##merge_bed   
    open(MR,">$outdir/somatic/$samplePair/called/shell/mergeNonLancetBed.$samplePair.sh") or die "$!\n";
    my $MergeBedSH="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";  
    $MergeBedSH.="export LD_LIBRARY_PATH=/share/app/gcc-4.9.3/lib64:/share/app/libz/zlib-1.2.11:\$LD_LIBRARY_PATH && \\\n";  
    $MergeBedSH.="cat $outdir/somatic/$samplePair/called/$samplePair.NonLancet.called.bed|sort -n -k 1.4,1.5 -k 2 -k 3 -u|$bedtools intersect -a - -b $DATAPATH/chrmosomes.bed > $outdir/somatic/$samplePair/called/$samplePair.targetForLancet.bed && \\\n";
    print MR $MergeBedSH;
    print MR echostring("$outdir/somatic/$samplePair/called/shell/mergeNonLancetBed.$samplePair.sh");
    close MR;
}
   
sub mergeMulTools{
    my ($normal,$tumor,$nbam,$tbam)=@_;
    my $add_para="";
    my $samplePair="$normal-VS-$tumor";
    `mkdir -p $outdir/somatic/$samplePair/called/multTools`;    
    foreach my $soft('MuSE','Strelka2','MuTect2','SomaticSniper','Lancet','Svaba'){
        if ($softwares{$soft}=~/true/i){
            if ($soft=~/MuSE|SomaticSniper/){
                $add_para.="-$soft $outdir/somatic/$samplePair/$soft/result/$samplePair.$soft.snv.vcf.gz ";
            }
            elsif($soft=~/MuTect2|Lancet/){
                $add_para.="-$soft $outdir/somatic/$samplePair/$soft/result/$samplePair.$soft.snv_shortindel.vcf.gz ";
            }
            elsif($soft=~/Strelka2/){
                $add_para.="-Strelka2_snv $outdir/somatic/$samplePair/$soft/result/$samplePair.$soft.snv.vcf.gz -Strelka2_indel $outdir/somatic/$samplePair/$soft/result/$samplePair.$soft.indel.vcf.gz ";
            }
            elsif($soft=~/Svaba/){
                $add_para.="-Svaba_indel $outdir/somatic/$samplePair/$soft/result/$samplePair.svaba.somatic.indel.vcf.gz ";
            }
        }
    }
    open(MTT,">$outdir/somatic/$samplePair/called/shell/multTools.$samplePair.sh") or die $!;
    my $MergeMutSH="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    $MergeMutSH.="export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:/hwfssz1/ST_CANCER/POL/SHARE/tools/bam-readcount/v0.8.0/lib:/share/app/gcc-4.9.3/lib64:\$LD_LIBRARY_PATH && \\\n";
    $MergeMutSH.="/share/app/perl-5.22.0/bin/perl $BIN/multTools_somatic_snvs_indels_ensemble-v1.0.pl -tName $tumor -nName $normal --minSnvTools 2 -minIndelTools 2 -seqType $seqType -ref $reference -tBAM $tbam -nBAM $nbam $add_para -outDir $outdir/somatic/$samplePair/called/multTools && \\\n";
    print MTT $MergeMutSH;
    print MTT echostring("$outdir/somatic/$samplePair/called/shell/multTools.$samplePair.sh");
    close MTT;
    
    #annovar 
    open ANNO, ">$outdir/somatic/$samplePair/called/shell/multTools.$samplePair.ann.sh" || die $!;
    my $AnnotationCmd="#!/bin/bash\nset -euo pipefail\necho ==========start at : `date` ==========\n";
    die "$annovar does not exist!\n" if (!-e "$annovar");
    die "$annovardb does not exist!\n" if (!-e "$annovardb"); 
    $AnnotationCmd.="/share/app/perl-5.22.0/bin/perl $annovar/table_annovar.pl $outdir/somatic/$samplePair/called/multTools/$samplePair.merged.norm.final.vcf.gz $annovardb -buildver hg38 -out $outdir/somatic/$samplePair/called/multTools/$samplePair.merged.norm.final -remove -protocol refGene,ensGene,knownGene,cytoBand,genomicSuperDups,gwasCatalog,phastConsElements100way,wgRna,avsnp150,clinvar_20190305,cosmic89,dbnsfp35c,exac03nonpsych,exac03nontcga,gnomad_exome,gnomad_genome,intervar_20180118 -operation g,g,g,r,r,r,r,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput && \\\n";
    $AnnotationCmd.="gzip -f $outdir/somatic/$samplePair/called/multTools/$samplePair.merged.norm.final.hg38_multianno.txt && \\\n";
    $AnnotationCmd.="rm -rf $outdir/somatic/$samplePair/called/multTools/$samplePair.merged.norm.final.avinput $outdir/somatic/$samplePair/called/multTools/$samplePair.merged.norm.final.hg38_multianno.vcf && \\\n";
    $AnnotationCmd.=echostring("$outdir/somatic/$samplePair/called/shell/multTools.$samplePair.ann.sh");
    print ANNO $AnnotationCmd;
    close ANNO;
}    
    
sub edgeList{
    `mkdir -p $outdir/shell_run/`;
    my $annMem="6G:1cpu";
    my $SentieonMem="100G:$SentieonCPU"."cpu:sentieon.q:$priority"."_sentieon";
    open EDGE,">$outdir/shell_run/edge.$subProjectName.list" or die $!;
	
    my %flagsh;
    
    for my $samp (sort keys %samples){
        my $dependsh;
        for my $lib (sort keys %{$samples{$samp}}){
            for my $lane_barcode (sort keys %{$samples{$samp}{$lib}}){
                print EDGE "$outdir/$samp/$lane_barcode/shell/clean.sh:6G:4cpu $outdir/$samp/$lane_barcode/shell/fqcheck.sh:3G:1cpu\n";
                print EDGE "$outdir/$samp/$lane_barcode/shell/clean.sh:6G:4cpu $outdir/$samp/$lane_barcode/shell/Sentieon.bwamem.sh:$SentieonMem\n";
                if($oneSentieon=~/false/i){print EDGE "$outdir/$samp/$lane_barcode/shell/Sentieon.bwamem.sh:$SentieonMem $outdir/$samp/shell/Sentieon.mkdup.$samp.sh:$SentieonMem\n";}
                elsif($oneSentieon=~/true/i){print EDGE "$outdir/$samp/$lane_barcode/shell/Sentieon.bwamem.sh:$SentieonMem $outdir/$samp/shell/Sentieon.mkdup.realign.bqsr.qc.$samp.sh:$SentieonMem\n";}
                if($rmCleanReads =~ /true/i){
                    print EDGE "$outdir/$samp/$lane_barcode/shell/fqcheck.sh:3G:1cpu $outdir/$samp/shell/rmCleanReads.$samp.sh:1G:1cpu\n";
                }
            }
		}
        if($oneSentieon=~/false/i){
            print EDGE "$outdir/$samp/shell/Sentieon.mkdup.$samp.sh:$SentieonMem $outdir/$samp/shell/Sentieon.realn.$samp.sh:$SentieonMem\n";
            print EDGE "$outdir/$samp/shell/Sentieon.realn.$samp.sh:$SentieonMem $outdir/$samp/shell/Sentieon.bqsr.$samp.sh:$SentieonMem\n";
            print EDGE "$outdir/$samp/shell/Sentieon.bqsr.$samp.sh:$SentieonMem $outdir/$samp/shell/Sentieon.qc.$samp.sh:$SentieonMem\n";

            $dependsh = "$outdir/$samp/shell/Sentieon.qc.$samp.sh:$SentieonMem";
        }
        elsif($oneSentieon=~/true/i){
            $dependsh = "$outdir/$samp/shell/Sentieon.mkdup.realign.bqsr.qc.$samp.sh:$SentieonMem";
        }

        if($seqType =~ /target/i){
			print EDGE "$dependsh $outdir/$samp/shell/covdep.$samp.sh:2G:1cpu\n";
            if($bam2cram =~ /true/i){print EDGE "$outdir/$samp/shell/covdep.$samp.sh:2G:1cpu $outdir/$samp/shell/Sentieon.bam2cram.$samp.sh:$SentieonMem\n"}
		}

        if($bam2cram =~ /true/i){print EDGE "$dependsh $outdir/$samp/shell/Sentieon.bam2cram.$samp.sh:$SentieonMem\n";}
        print EDGE "$dependsh $outdir/$samp/shell/variant_Ann.$samp.sh:$annMem\n";
		
        print EDGE "$dependsh $outdir/$samp/shell/bam_readcount.$samp.sh:1G:1cpu\n";
        if($bam2cram =~ /true/i){print EDGE "$outdir/$samp/shell/bam_readcount.$samp.sh:1G:1cpu $outdir/$samp/shell/Sentieon.bam2cram.$samp.sh:$SentieonMem\n";}
        print EDGE "$dependsh $outdir/$samp/QC/shell/ContEst.$samp.sh:7G:1cpu\n";
        if($bam2cram =~ /true/i){print EDGE "$outdir/$samp/QC/shell/ContEst.$samp.sh:7G:1cpu $outdir/$samp/shell/Sentieon.bam2cram.$samp.sh:$SentieonMem\n";}

=header		
        print EDGE "$dependsh $outdir/$samp/QC/shell/flagstat.$samp.sh:3G:1cpu\n";
        print EDGE "$outdir/$samp/QC/shell/flagstat.$samp.sh:3G:1cpu $outdir/$samp/QC/shell/summary.$samp.sh:1G:1cpu\n";
		
        print EDGE "$dependsh $outdir/$samp/QC/shell/coverageQc.$samp.sh:3G:1cpu\n";
        print EDGE "$outdir/$samp/QC/shell/coverageQc.$samp.sh:3G:1cpu $outdir/$samp/QC/shell/summary.$samp.sh:1G:1cpu\n";
		
        print EDGE "$dependsh $outdir/$samp/QC/shell/genomeCoverage.$samp.sh:3G:4cpu\n";
        print EDGE "$outdir/$samp/QC/shell/genomeCoverage.$samp.sh:3G:4cpu $outdir/$samp/QC/shell/summary.$samp.sh:1G:1cpu\n";
		
        print EDGE "$dependsh $outdir/$samp/QC/shell/bam_stats.$samp.sh:3G:1cpu\n";
        print EDGE "$outdir/$samp/QC/shell/bam_stats.$samp.sh:3G:1cpu $outdir/$samp/QC/shell/summary.$samp.sh:1G:1cpu\n";

        if($bam2cram =~ /true/i){print EDGE "$outdir/$samp/QC/shell/summary.$samp.sh:1G:1cpu $outdir/$samp/shell/Sentieon.bam2cram.$samp.sh:$SentieonMem\n";}
=cut

        if($rmCleanReads =~ /true/i){
            print EDGE "$dependsh $outdir/$samp/shell/rmCleanReads.$samp.sh:1G:1cpu\n";
            if($bam2cram =~ /true/i){print EDGE "$outdir/$samp/shell/rmCleanReads.$samp.sh:1G:1cpu $outdir/$samp/shell/Sentieon.bam2cram.$samp.sh:$SentieonMem\n";}
        }
    }
      
    foreach my $samp(sort keys %ALL_SAMP){
        my $dependsh;
        if($somatic && exists $hash_sampair{$samp}){
            if(exists $samples{$samp}){
                if($oneSentieon=~/false/i){$dependsh = "$outdir/$samp/shell/Sentieon.bqsr.$samp.sh:$SentieonMem";}
                elsif($oneSentieon=~/true/i){$dependsh = "$outdir/$samp/shell/Sentieon.mkdup.realign.bqsr.qc.$samp.sh:$SentieonMem";}
            }
            elsif(exists $cram2bam{$samp}){
                $dependsh = $cram2bam{$samp}{"S"}.":$SentieonMem";
            }
            else{
				$dependsh = "$outdir/TR/shell/pipeline_with_existed_bam.sh:1G:1cpu";
			}
			
			foreach my $TRsplitNum(sort keys %TRs){
				print EDGE "$outdir/TR/shell/split_intervals.$TRsplitNum.sh:1G:1cpu $dependsh\n";
			}
            for my $samp2 (sort keys %{$hash_sampair{$samp}}){
                $flagsh{$samp}{$samp2} += 1;
                $flagsh{$samp2}{$samp} += 1;
                my $samplePair = $hash_sampair{$samp}{$samp2};
				
                print EDGE "$dependsh $outdir/somatic/$samplePair/QC/somaticQC.sh:5G:1cpu\n";
                if($bam2cram =~ /true/i){print EDGE "$outdir/somatic/$samplePair/QC/somaticQC.sh:5G:1cpu $outdir/$samp/shell/Sentieon.bam2cram.$samp.sh:$SentieonMem\n";}
                if(exists $cram2bam{$samp}){print EDGE "$outdir/somatic/$samplePair/QC/somaticQC.sh:5G:1cpu $cram2bam{$samp}{E}:1G:1cpu\n";}
                if($Use_TNcaller =~ /true/i){
                    print EDGE "$dependsh $outdir/somatic/$samplePair/Sentieon/shell/Sentieon.TNcaller.$samplePair.sh:$SentieonMem\n";
                    if($flagsh{$samp}{$samp2}==2){
                        if($bam2cram =~ /true/i){
                            if(exists $samples{$samp}){print EDGE "$outdir/somatic/$samplePair/Sentieon/shell/Sentieon.TNcaller.$samplePair.sh:$SentieonMem $outdir/$samp/shell/Sentieon.bam2cram.$samp.sh:$SentieonMem\n";}
                            if(exists $samples{$samp2}){print EDGE "$outdir/somatic/$samplePair/Sentieon/shell/Sentieon.TNcaller.$samplePair.sh:$SentieonMem $outdir/$samp2/shell/Sentieon.bam2cram.$samp2.sh:$SentieonMem\n";}
                        }
                    }
                    if(exists $cram2bam{$samp}){print EDGE "$outdir/somatic/$samplePair/Sentieon/shell/Sentieon.TNcaller.$samplePair.sh:$SentieonMem $cram2bam{$samp}{E}:1G:1cpu\n";}
                }
                                
                if($Use_MuTect =~ /true/i){ ## snv calling by MuTect1
                    if ($seqType =~ /target/i && $MuTect_TRsplitNum eq 1){
                        print EDGE "$dependsh $outdir/somatic/$samplePair/MuTect/shell/MuTect.$samplePair.sh:25G:2cpu\n";
                        if($flagsh{$samp}{$samp2}==2){
                            print EDGE "$outdir/somatic/$samplePair/MuTect/shell/MuTect.$samplePair.sh:25G:2cpu $outdir/somatic/$samplePair/MuTect/shell/MuTect.$samplePair.process.sh:$annMem\n";
                        }
                    }
                    elsif(($seqType =~ /wgs/i) || ($seqType =~ /target/i && $MuTect_TRsplitNum > 1)){
                        my $Num=$MuTect_TRsplitNum-1;
                        foreach my $n(0 .. $Num){
                            print EDGE "$dependsh $outdir/somatic/$samplePair/MuTect/shell/MuTect.$samplePair.$n.sh:25G:2cpu\n";
                            if($flagsh{$samp}{$samp2}==2){
                                print EDGE "$outdir/somatic/$samplePair/MuTect/shell/MuTect.$samplePair.$n.sh:25G:2cpu $outdir/somatic/$samplePair/MuTect/shell/MuTect.$samplePair.process.sh:$annMem\n";
                            }
                        }
                    }
                    if($flagsh{$samp}{$samp2}==2){
                        print EDGE "$outdir/somatic/$samplePair/MuTect/shell/MuTect.$samplePair.process.sh:$annMem $outdir/somatic/$samplePair/MuTect/shell/MuTect.$samplePair.ann.sh:$annMem\n";
                        if($bam2cram =~ /true/i){
                            if(exists $samples{$samp}){print EDGE "$outdir/somatic/$samplePair/MuTect/shell/MuTect.$samplePair.ann.sh:$annMem $outdir/$samp/shell/Sentieon.bam2cram.$samp.sh:$SentieonMem\n";}
                            if(exists $samples{$samp2}){print EDGE "$outdir/somatic/$samplePair/MuTect/shell/MuTect.$samplePair.ann.sh:$annMem $outdir/$samp2/shell/Sentieon.bam2cram.$samp2.sh:$SentieonMem\n";}
                        }                            
                    }
                    if(exists $cram2bam{$samp}){print EDGE "$outdir/somatic/$samplePair/MuTect/shell/MuTect.$samplePair.ann.sh:$annMem $cram2bam{$samp}{E}:1G:1cpu\n";}
                }                               
                if($Use_MuTect2 =~ /true/i){ ## snv calling by MuTect2 
                    if ($seqType =~ /target/i && $MuTect2_TRsplitNum eq 1){
                        print EDGE "$dependsh $outdir/somatic/$samplePair/MuTect2/shell/MuTect2.$samplePair.sh:10G:1cpu\n";
                        if($flagsh{$samp}{$samp2}==2){print EDGE "$outdir/somatic/$samplePair/MuTect2/shell/MuTect2.$samplePair.sh:10G:1cpu $outdir/somatic/$samplePair/MuTect2/shell/MuTect2.$samplePair.process.sh:$annMem\n";}
                    }
                    elsif(($seqType =~ /wgs/i) || ($seqType =~ /target/i && $MuTect2_TRsplitNum > 1)){
                        my $Num=$MuTect2_TRsplitNum-1;
                        foreach my $n(0 .. $Num){
                            print EDGE "$dependsh $outdir/somatic/$samplePair/MuTect2/shell/MuTect2.$samplePair.$n.sh:10G:1cpu\n";
                            if($flagsh{$samp}{$samp2}==2){print EDGE "$outdir/somatic/$samplePair/MuTect2/shell/MuTect2.$samplePair.$n.sh:10G:1cpu $outdir/somatic/$samplePair/MuTect2/shell/MuTect2.$samplePair.process.sh:$annMem\n";
                            }
                        }
                    }
                    if($flagsh{$samp}{$samp2}==2){
                        print EDGE "$outdir/somatic/$samplePair/MuTect2/shell/MuTect2.$samplePair.process.sh:$annMem $outdir/somatic/$samplePair/MuTect2/shell/MuTect2.$samplePair.getbed.sh:2G:1cpu\n";
                        print EDGE "$outdir/somatic/$samplePair/MuTect2/shell/MuTect2.$samplePair.getbed.sh:2G:1cpu $outdir/somatic/$samplePair/called/shell/mergeNonLancetBed.$samplePair.sh:2G:1cpu\n";
                        print EDGE "$outdir/somatic/$samplePair/MuTect2/shell/MuTect2.$samplePair.process.sh:$annMem $outdir/somatic/$samplePair/MuTect2/shell/MuTect2.$samplePair.ann.sh:$annMem\n";
                        print EDGE "$outdir/somatic/$samplePair/MuTect2/shell/MuTect2.$samplePair.process.sh:$annMem $outdir/somatic/$samplePair/called/shell/multTools.$samplePair.sh:20G:2cpu\n";
                        print EDGE "$outdir/somatic/$samplePair/called/shell/multTools.$samplePair.sh:20G:2cpu $outdir/somatic/$samplePair/called/shell/multTools.$samplePair.ann.sh:$annMem\n";                            
                        if($bam2cram =~ /true/i){
                            if(exists $samples{$samp}){
                                print EDGE "$outdir/somatic/$samplePair/MuTect2/shell/MuTect2.$samplePair.ann.sh:$annMem $outdir/$samp/shell/Sentieon.bam2cram.$samp.sh:$SentieonMem\n";
                                print EDGE "$outdir/somatic/$samplePair/called/shell/multTools.$samplePair.ann.sh:$annMem $outdir/$samp/shell/Sentieon.bam2cram.$samp.sh:$SentieonMem\n";
                            }
                            if(exists $samples{$samp2}){
                                print EDGE "$outdir/somatic/$samplePair/MuTect2/shell/MuTect2.$samplePair.ann.sh:$annMem $outdir/$samp2/shell/Sentieon.bam2cram.$samp2.sh:$SentieonMem\n";
                                print EDGE "$outdir/somatic/$samplePair/called/shell/multTools.$samplePair.ann.sh:$annMem $outdir/$samp2/shell/Sentieon.bam2cram.$samp2.sh:$SentieonMem\n";
                            }
                        }
                    }
                    if(exists $cram2bam{$samp}){print EDGE "$outdir/somatic/$samplePair/MuTect2/shell/MuTect2.$samplePair.ann.sh:$annMem $cram2bam{$samp}{E}:1G:1cpu\n";}
                }                
                if($Use_MuSE =~ /true/i){ ## snv calling by MuSE
                    if ($seqType =~ /target/i && $MuSE_TRsplitNum eq 1){
                        print EDGE "$dependsh $outdir/somatic/$samplePair/MuSE/shell/MuSE.$samplePair.sh:4G:1cpu\n";
                        if($flagsh{$samp}{$samp2}==2){print EDGE "$outdir/somatic/$samplePair/MuSE/shell/MuSE.$samplePair.sh:4G:1cpu $outdir/somatic/$samplePair/MuSE/shell/MuSE.$samplePair.process.sh:$annMem\n";}
                    }
                    elsif(($seqType =~ /wgs/i) || ($seqType =~ /target/i && $MuSE_TRsplitNum > 1)){
                        my $Num=$MuSE_TRsplitNum-1;
                        foreach my $n(0 .. $Num){
                            print EDGE "$dependsh $outdir/somatic/$samplePair/MuSE/shell/MuSE.$samplePair.$n.sh:4G:1cpu\n";
                            if($flagsh{$samp}{$samp2}==2){print EDGE "$outdir/somatic/$samplePair/MuSE/shell/MuSE.$samplePair.$n.sh:4G:1cpu $outdir/somatic/$samplePair/MuSE/shell/MuSE.$samplePair.process.sh:$annMem\n";
                            }
                        }
                    }
                    if($flagsh{$samp}{$samp2}==2){
                        print EDGE "$outdir/somatic/$samplePair/MuSE/shell/MuSE.$samplePair.process.sh:$annMem $outdir/somatic/$samplePair/MuSE/shell/MuSE.$samplePair.getbed.sh:2G:1cpu\n";
                        print EDGE "$outdir/somatic/$samplePair/MuSE/shell/MuSE.$samplePair.getbed.sh:2G:1cpu $outdir/somatic/$samplePair/called/shell/mergeNonLancetBed.$samplePair.sh:2G:1cpu\n";
                        print EDGE "$outdir/somatic/$samplePair/MuSE/shell/MuSE.$samplePair.process.sh:$annMem $outdir/somatic/$samplePair/MuSE/shell/MuSE.$samplePair.ann.sh:$annMem\n";
                        print EDGE "$outdir/somatic/$samplePair/MuSE/shell/MuSE.$samplePair.process.sh:$annMem $outdir/somatic/$samplePair/called/shell/multTools.$samplePair.sh:20G:2cpu\n";
                        print EDGE "$outdir/somatic/$samplePair/called/shell/multTools.$samplePair.sh:20G:2cpu $outdir/somatic/$samplePair/called/shell/multTools.$samplePair.ann.sh:$annMem\n";
                        if($bam2cram =~ /true/i){
                            if(exists $samples{$samp}){
                                print EDGE "$outdir/somatic/$samplePair/MuSE/shell/MuSE.$samplePair.ann.sh:$annMem $outdir/$samp/shell/Sentieon.bam2cram.$samp.sh:$SentieonMem\n";
                                print EDGE "$outdir/somatic/$samplePair/called/shell/multTools.$samplePair.ann.sh:$annMem $outdir/$samp/shell/Sentieon.bam2cram.$samp.sh:$SentieonMem\n";
                            }
                            if(exists $samples{$samp2}){
                                print EDGE "$outdir/somatic/$samplePair/MuSE/shell/MuSE.$samplePair.ann.sh:$annMem $outdir/$samp2/shell/Sentieon.bam2cram.$samp2.sh:$SentieonMem\n";
                                print EDGE "$outdir/somatic/$samplePair/called/shell/multTools.$samplePair.ann.sh:$annMem $outdir/$samp2/shell/Sentieon.bam2cram.$samp2.sh:$SentieonMem\n";
                            }
                        }                   
                    }
                    if(exists $cram2bam{$samp}){print EDGE "$outdir/somatic/$samplePair/MuSE/shell/MuSE.$samplePair.ann.sh:$annMem $cram2bam{$samp}{E}:1G:1cpu\n";}
                }                
                if($Use_SomaticSniper =~ /true/i){ ## snv calling by SomaticSniper 
                    print EDGE "$dependsh $outdir/somatic/$samplePair/SomaticSniper/shell/SomaticSniper.$samplePair.sh:4G:1cpu\n";
                    if($flagsh{$samp}{$samp2}==2){
                        if ($seqType =~ /target/i && $SomaticSniper_TRsplitNum eq 1){
                            print EDGE "$outdir/somatic/$samplePair/SomaticSniper/shell/SomaticSniper.$samplePair.sh:4G:1cpu $outdir/somatic/$samplePair/SomaticSniper/shell/SomaticSniper.$samplePair.mpileup.sh:2G:1cpu\n";
                            print EDGE "$outdir/somatic/$samplePair/SomaticSniper/shell/SomaticSniper.$samplePair.mpileup.sh:2G:1cpu $outdir/somatic/$samplePair/SomaticSniper/shell/SomaticSniper.$samplePair.filter.sh:2G:1cpu\n";
                        }
                        elsif(($seqType =~ /wgs/i) || ($seqType =~ /target/i && $SomaticSniper_TRsplitNum > 1)){
                            my $Num=$SomaticSniper_TRsplitNum-1;
                            foreach my $n(0..$Num){
                                print EDGE "$outdir/somatic/$samplePair/SomaticSniper/shell/SomaticSniper.$samplePair.sh:4G:1cpu $outdir/somatic/$samplePair/SomaticSniper/shell/SomaticSniper.$samplePair.$n.mpileup.sh:2G:1cpu\n";
                                print EDGE "$outdir/somatic/$samplePair/SomaticSniper/shell/SomaticSniper.$samplePair.$n.mpileup.sh:2G:1cpu $outdir/somatic/$samplePair/SomaticSniper/shell/SomaticSniper.$samplePair.filter.sh:2G:1cpu\n";
                            }
                        }                                                
                        print EDGE "$outdir/somatic/$samplePair/SomaticSniper/shell/SomaticSniper.$samplePair.filter.sh:2G:1cpu $outdir/somatic/$samplePair/SomaticSniper/shell/SomaticSniper.$samplePair.getbed.sh:2G:1cpu\n";
                        print EDGE "$outdir/somatic/$samplePair/SomaticSniper/shell/SomaticSniper.$samplePair.getbed.sh:2G:1cpu $outdir/somatic/$samplePair/called/shell/mergeNonLancetBed.$samplePair.sh:2G:1cpu\n";
                        print EDGE "$outdir/somatic/$samplePair/SomaticSniper/shell/SomaticSniper.$samplePair.filter.sh:2G:1cpu $outdir/somatic/$samplePair/SomaticSniper/shell/SomaticSniper.$samplePair.ann.sh:$annMem\n";
                        print EDGE "$outdir/somatic/$samplePair/SomaticSniper/shell/SomaticSniper.$samplePair.filter.sh:2G:1cpu $outdir/somatic/$samplePair/called/shell/multTools.$samplePair.sh:20G:2cpu\n";
                        print EDGE "$outdir/somatic/$samplePair/called/shell/multTools.$samplePair.sh:20G:2cpu $outdir/somatic/$samplePair/called/shell/multTools.$samplePair.ann.sh:$annMem\n";
                        if($bam2cram =~ /true/i){
                            if(exists $samples{$samp}){
                                print EDGE "$outdir/somatic/$samplePair/SomaticSniper/shell/SomaticSniper.$samplePair.ann.sh:$annMem $outdir/$samp/shell/Sentieon.bam2cram.$samp.sh:$SentieonMem\n";
                                print EDGE "$outdir/somatic/$samplePair/called/shell/multTools.$samplePair.ann.sh:$annMem $outdir/$samp/shell/Sentieon.bam2cram.$samp.sh:$SentieonMem\n";
                            }
                            if(exists $samples{$samp2}){
                                print EDGE "$outdir/somatic/$samplePair/SomaticSniper/shell/SomaticSniper.$samplePair.ann.sh:$annMem $outdir/$samp2/shell/Sentieon.bam2cram.$samp2.sh:$SentieonMem\n";
                                print EDGE "$outdir/somatic/$samplePair/called/shell/multTools.$samplePair.ann.sh:$annMem $outdir/$samp2/shell/Sentieon.bam2cram.$samp2.sh:$SentieonMem\n";
                            }
                        }
                    }
                    if(exists $cram2bam{$samp}){print EDGE "$outdir/somatic/$samplePair/SomaticSniper/shell/SomaticSniper.$samplePair.ann.sh:$annMem $cram2bam{$samp}{E}:1G:1cpu\n";}
                }               
                if($Use_Strelka2 =~ /true/i){ ## snv calling by Strelka2
                    if ($seqType =~ /target/i && $Strelka2_TRsplitNum eq 1){
                        print EDGE "$dependsh $outdir/somatic/$samplePair/Strelka2/shell/Strelka2.$samplePair.sh:10G:4cpu\n";
                        if($flagsh{$samp}{$samp2}==2){
                            print EDGE "$outdir/somatic/$samplePair/Manta/shell/Manta.$samplePair.sh:45G:8cpu $outdir/somatic/$samplePair/Strelka2/shell/Strelka2.$samplePair.sh:10G:4cpu\n" if ($Use_MantaAndStrelka2 =~ /true/i);
                            print EDGE "$outdir/somatic/$samplePair/Strelka2/shell/Strelka2.$samplePair.sh:10G:4cpu $outdir/somatic/$samplePair/Strelka2/shell/Strelka2.$samplePair.process.sh:$annMem\n";
                        }
                    }
                    elsif(($seqType =~ /wgs/i) || ($seqType =~ /target/i && $Strelka2_TRsplitNum > 1)){
                        my $Num=$Strelka2_TRsplitNum-1;                    
                        foreach my $n(0 .. $Num){
                            print EDGE "$dependsh $outdir/somatic/$samplePair/Strelka2/shell/Strelka2.$samplePair.$n.sh:10G:4cpu\n";
                            if($flagsh{$samp}{$samp2}==2){print EDGE "$outdir/somatic/$samplePair/Manta/shell/Manta.$samplePair.sh:45G:8cpu $outdir/somatic/$samplePair/Strelka2/shell/Strelka2.$samplePair.$n.sh:10G:4cpu\n" if ($Use_MantaAndStrelka2 =~ /true/i);print EDGE "$outdir/somatic/$samplePair/Strelka2/shell/Strelka2.$samplePair.$n.sh:10G:4cpu $outdir/somatic/$samplePair/Strelka2/shell/Strelka2.$samplePair.process.sh:$annMem\n";
                            }
                        }
                    }
                    if($flagsh{$samp}{$samp2}==2){
                        print EDGE "$outdir/somatic/$samplePair/Strelka2/shell/Strelka2.$samplePair.process.sh:$annMem $outdir/somatic/$samplePair/Strelka2/shell/Strelka2.$samplePair.getbed.sh:2G:1cpu\n";
                        print EDGE "$outdir/somatic/$samplePair/Strelka2/shell/Strelka2.$samplePair.getbed.sh:2G:1cpu $outdir/somatic/$samplePair/called/shell/mergeNonLancetBed.$samplePair.sh:2G:1cpu\n";
                        print EDGE "$outdir/somatic/$samplePair/Strelka2/shell/Strelka2.$samplePair.process.sh:$annMem $outdir/somatic/$samplePair/Strelka2/shell/Strelka2.$samplePair.ann.sh:$annMem\n";
                        print EDGE "$outdir/somatic/$samplePair/Strelka2/shell/Strelka2.$samplePair.process.sh:$annMem $outdir/somatic/$samplePair/called/shell/multTools.$samplePair.sh:20G:2cpu\n";
                        print EDGE "$outdir/somatic/$samplePair/called/shell/multTools.$samplePair.sh:20G:2cpu $outdir/somatic/$samplePair/called/shell/multTools.$samplePair.ann.sh:$annMem\n";
                        if($bam2cram =~ /true/i){
                            if(exists $samples{$samp}){
                                print EDGE "$outdir/somatic/$samplePair/Strelka2/shell/Strelka2.$samplePair.ann.sh:$annMem $outdir/$samp/shell/Sentieon.bam2cram.$samp.sh:$SentieonMem\n";
                                print EDGE "$outdir/somatic/$samplePair/called/shell/multTools.$samplePair.ann.sh:$annMem $outdir/$samp/shell/Sentieon.bam2cram.$samp.sh:$SentieonMem\n";
                            }
                            if(exists $samples{$samp2}){
                                print EDGE "$outdir/somatic/$samplePair/Strelka2/shell/Strelka2.$samplePair.ann.sh:$annMem $outdir/$samp2/shell/Sentieon.bam2cram.$samp2.sh:$SentieonMem\n";
                                print EDGE "$outdir/somatic/$samplePair/called/shell/multTools.$samplePair.ann.sh:$annMem $outdir/$samp2/shell/Sentieon.bam2cram.$samp2.sh:$SentieonMem\n";    
                            }
                        }                    
                    }
                    if(exists $cram2bam{$samp}){print EDGE "$outdir/somatic/$samplePair/Strelka2/shell/Strelka2.$samplePair.ann.sh:$annMem $cram2bam{$samp}{E}:1G:1cpu\n";}
                }               
                if($Use_Lancet =~ /true/i){
                    print EDGE "$dependsh $outdir/somatic/$samplePair/Lancet/shell/Lancet.$samplePair.sh:4G:1cpu\n";
                    if($flagsh{$samp}{$samp2}==1){print EDGE "$outdir/somatic/$samplePair/called/shell/mergeNonLancetBed.$samplePair.sh:2G:1cpu $outdir/somatic/$samplePair/Lancet/shell/Lancet.$samplePair.sh:10G:4cpu\n";}
                    if($flagsh{$samp}{$samp2}==2){
                        print EDGE "$outdir/somatic/$samplePair/Lancet/shell/Lancet.$samplePair.sh:10G:4cpu $outdir/somatic/$samplePair/Lancet/shell/Lancet.$samplePair.getbed.sh:2G:1cpu\n";
                        print EDGE "$outdir/somatic/$samplePair/Lancet/shell/Lancet.$samplePair.sh:10G:4cpu $outdir/somatic/$samplePair/Lancet/shell/Lancet.$samplePair.ann.sh:$annMem\n";
                        print EDGE "$outdir/somatic/$samplePair/Lancet/shell/Lancet.$samplePair.sh:10G:4cpu $outdir/somatic/$samplePair/called/shell/multTools.$samplePair.sh:20G:2cpu\n";
                        print EDGE "$outdir/somatic/$samplePair/called/shell/multTools.$samplePair.sh:20G:2cpu $outdir/somatic/$samplePair/called/shell/multTools.$samplePair.ann.sh:$annMem\n";
                        if($bam2cram =~ /true/i){
                            if(exists $samples{$samp}){
                                print EDGE "$outdir/somatic/$samplePair/Lancet/shell/Lancet.$samplePair.ann.sh:$annMem $outdir/$samp/shell/Sentieon.bam2cram.$samp.sh:$SentieonMem\n";
                                print EDGE "$outdir/somatic/$samplePair/called/shell/multTools.$samplePair.ann.sh:$annMem $outdir/$samp/shell/Sentieon.bam2cram.$samp.sh:$SentieonMem\n";
                            }
                            if(exists $samples{$samp2}){
                                print EDGE "$outdir/somatic/$samplePair/Lancet/shell/Lancet.$samplePair.ann.sh:$annMem $outdir/$samp2/shell/Sentieon.bam2cram.$samp2.sh:$SentieonMem\n";
                                print EDGE "$outdir/somatic/$samplePair/called/shell/multTools.$samplePair.ann.sh:$annMem $outdir/$samp2/shell/Sentieon.bam2cram.$samp2.sh:$SentieonMem\n";
                            }
                        }
                    }
                    if(exists $cram2bam{$samp}){print EDGE "$outdir/somatic/$samplePair/Lancet/shell/Lancet.$samplePair.ann.sh:$annMem $cram2bam{$samp}{E}:1G:1cpu\n";}
                }                                
                
                if($Use_Platypus =~ /true/i){ ## somatic indel calling by platypus
                    if ($seqType =~ /target/i && $Platypus_TRsplitNum eq 1){
                        print EDGE "$dependsh $outdir/somatic/$samplePair/Platypus/shell/Platypus.$samplePair.sh:10G:4cpu\n";
                        if($flagsh{$samp}{$samp2}==2){print EDGE "$outdir/somatic/$samplePair/Platypus/shell/Platypus.$samplePair.sh:10G:4cpu $outdir/somatic/$samplePair/Platypus/shell/Platypus.$samplePair.process.sh:$annMem\n";}
                    }
                    elsif(($seqType =~ /wgs/i) || ($seqType =~ /target/i && $Platypus_TRsplitNum > 1)){
                        my $Num=$Platypus_TRsplitNum-1;
                        foreach my $n(0 .. $Num){
                            print EDGE "$dependsh $outdir/somatic/$samplePair/Platypus/shell/Platypus.$samplePair.$n.sh:10G:4cpu\n";
                            if($flagsh{$samp}{$samp2}==2){print EDGE "$outdir/somatic/$samplePair/Platypus/shell/Platypus.$samplePair.$n.sh:10G:4cpu $outdir/somatic/$samplePair/Platypus/shell/Platypus.$samplePair.process.sh:$annMem\n";
                            }
                        }
                    }
                    if($flagsh{$samp}{$samp2}==2){
                        print EDGE "$outdir/somatic/$samplePair/Platypus/shell/Platypus.$samplePair.process.sh:$annMem $outdir/somatic/$samplePair/Platypus/shell/Platypus.$samplePair.ann.sh:$annMem\n";
                        if($bam2cram =~ /true/i){
                            if(exists $samples{$samp}){print EDGE "$outdir/somatic/$samplePair/Platypus/shell/Platypus.$samplePair.ann.sh:$annMem $outdir/$samp/shell/Sentieon.bam2cram.$samp.sh:$SentieonMem\n";}
                            if(exists $samples{$samp2}){print EDGE "$outdir/somatic/$samplePair/Platypus/shell/Platypus.$samplePair.ann.sh:$annMem $outdir/$samp2/shell/Sentieon.bam2cram.$samp2.sh:$SentieonMem\n";}
                        }                    
                    }
                    if(exists $cram2bam{$samp}){print EDGE "$outdir/somatic/$samplePair/Platypus/shell/Platypus.$samplePair.ann.sh:$annMem $cram2bam{$samp}{E}:1G:1cpu\n";}
                }                
                if($Use_FACETS =~ /true/i){ ## somatic cnv calling by FACETS
                    print EDGE "$dependsh $outdir/somatic/$samplePair/FACETS/shell/FACETS.$samplePair.pileup.sh:3G:1cpu\n";
                    if($flagsh{$samp}{$samp2}==2){
                        print EDGE "$outdir/somatic/$samplePair/FACETS/shell/FACETS.$samplePair.pileup.sh:3G:1cpu $outdir/somatic/$samplePair/FACETS/shell/FACETS.$samplePair.calling.sh:3G:1cpu\n";
                        print EDGE "$outdir/somatic/$samplePair/FACETS/shell/FACETS.$samplePair.calling.sh:3G:1cpu $outdir/somatic/$samplePair/FACETS/shell/FACETS.$samplePair.ann.sh:$annMem\n";
                        if($bam2cram =~ /true/i){
                            if(exists $samples{$samp}){print EDGE "$outdir/somatic/$samplePair/FACETS/shell/FACETS.$samplePair.ann.sh:$annMem $outdir/$samp/shell/Sentieon.bam2cram.$samp.sh:$SentieonMem\n";}
                            if(exists $samples{$samp2}){print EDGE "$outdir/somatic/$samplePair/FACETS/shell/FACETS.$samplePair.ann.sh:$annMem $outdir/$samp2/shell/Sentieon.bam2cram.$samp2.sh:$SentieonMem\n";}
                        }
                    }
                    if(exists $cram2bam{$samp}){print EDGE "$outdir/somatic/$samplePair/FACETS/shell/FACETS.$samplePair.ann.sh:$annMem $cram2bam{$samp}{E}:1G:1cpu\n";}
                }
                if($Use_Manta =~ /true/i){ ## somatic SVs and indels calling by Manta
                    print EDGE "$dependsh $outdir/somatic/$samplePair/Manta/shell/Manta.$samplePair.sh:45G:8cpu\n";
                    if($flagsh{$samp}{$samp2}==2){
                        print EDGE "$outdir/somatic/$samplePair/Manta/shell/Manta.$samplePair.sh:45G:8cpu $outdir/somatic/$samplePair/Manta/shell/Manta.$samplePair.ann.sh:$annMem\n";
                            if($bam2cram =~ /true/i){
                                if(exists $samples{$samp}){print EDGE "$outdir/somatic/$samplePair/Manta/shell/Manta.$samplePair.ann.sh:$annMem $outdir/$samp/shell/Sentieon.bam2cram.$samp.sh:$SentieonMem\n";}
                                if(exists $samples{$samp2}){print EDGE "$outdir/somatic/$samplePair/Manta/shell/Manta.$samplePair.ann.sh:$annMem $outdir/$samp2/shell/Sentieon.bam2cram.$samp2.sh:$SentieonMem\n";}
                            }
                        }
                    if(exists $cram2bam{$samp}){print EDGE "$outdir/somatic/$samplePair/Manta/shell/Manta.$samplePair.ann.sh:$annMem $cram2bam{$samp}{E}:1G:1cpu\n";}
                }
                    
                if($Use_Svaba =~ /true/i){ ## somatic SVs and indels calling by Svaba
                    print EDGE "$dependsh $outdir/somatic/$samplePair/Svaba/shell/Svaba.$samplePair.sh:45G:20cpu\n";
                    if($flagsh{$samp}{$samp2}==2){
                        print EDGE "$outdir/somatic/$samplePair/Svaba/shell/Svaba.$samplePair.sh:45G:20cpu $outdir/somatic/$samplePair/Svaba/shell/Svaba.$samplePair.ann.sh:$annMem\n";
                        print EDGE "$outdir/somatic/$samplePair/Svaba/shell/Svaba.$samplePair.sh:45G:20cpu $outdir/somatic/$samplePair/called/shell/multTools.$samplePair.sh:20G:2cpu\n";
                        print EDGE "$outdir/somatic/$samplePair/called/shell/multTools.$samplePair.sh:20G:2cpu $outdir/somatic/$samplePair/called/shell/multTools.$samplePair.ann.sh:$annMem\n";
                        if($bam2cram =~ /true/i){
                            if(exists $samples{$samp}){
                                print EDGE "$outdir/somatic/$samplePair/Svaba/shell/Svaba.$samplePair.ann.sh:$annMem $outdir/$samp/shell/Sentieon.bam2cram.$samp.sh:$SentieonMem\n";
                                print EDGE "$outdir/somatic/$samplePair/called/shell/multTools.$samplePair.ann.sh:$annMem $outdir/$samp/shell/Sentieon.bam2cram.$samp.sh:$SentieonMem\n";
                            }
                            if(exists $samples{$samp2}){
                                print EDGE "$outdir/somatic/$samplePair/Svaba/shell/Svaba.$samplePair.ann.sh:$annMem $outdir/$samp2/shell/Sentieon.bam2cram.$samp2.sh:$SentieonMem\n";
                                print EDGE "$outdir/somatic/$samplePair/called/shell/multTools.$samplePair.ann.sh:$annMem $outdir/$samp2/shell/Sentieon.bam2cram.$samp2.sh:$SentieonMem\n";
                            }
                        }
                    }
                    if(exists $cram2bam{$samp}){print EDGE "$outdir/somatic/$samplePair/Svaba/shell/Svaba.$samplePair.ann.sh:$annMem $cram2bam{$samp}{E}:1G:1cpu\n";}
                }
                if(exists $cram2bam{$samp}){print EDGE "$outdir/somatic/$samplePair/called/shell/multTools.$samplePair.ann.sh:$annMem $cram2bam{$samp}{E}:1G:1cpu\n";}                 
            }
        }
    }
    close EDGE;

    open RUNSH,">$outdir/shell_run/run.$subProjectName.sh" or die $!;
    print RUNSH "/hwfssz1/ST_CANCER/POL/SHARE/tools/PyMonitor/v1.6/monitor2.7 cron -m 4\n";
    print RUNSH "/hwfssz1/ST_CANCER/POL/SHARE/tools/PyMonitor/v1.6/monitor2.7 cron -m 5\n";
    print RUNSH "/hwfssz1/ST_CANCER/POL/SHARE/tools/PyMonitor/v1.6/monitor2.7 taskmonitor -q $queue -P $priority -p $subProjectName -i $outdir/shell_run/edge.$subProjectName.list\n";
    close RUNSH;
}
