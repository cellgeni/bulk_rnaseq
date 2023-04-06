#!/bin/bash 

## all-in-one processing script that should work as individual submission as well as part of job array
## all necessary tools are packaged in Singularity
## everything here is for paired-end reads, single end would require many option changes 

SIF="/nfs/cellgeni/singularity/images/starsolo_2-7-10a-alpha-220818_samtools_1-15-1_seqtk-1-13_bbmap_38-97_RSEM-1-3-3.sif"
CMD="singularity run --nv --bind /nfs,/lustre,/software $SIF"

FQDIR=$1
TAG=$2
CPUS=16

FQDIR=`readlink -f $FQDIR`
SREF=/nfs/cellgeni/STAR/human/Gencode_v19_full/index
RREF=/nfs/cellgeni/STAR/human/Gencode_v19_full/Gencode_v19_rsem 
GTF=/nfs/cellgeni/STAR/human/Gencode_v19_full/gencode.v19.chr_patch_hapl_scaff.annotation.gtf

## step 1 - adapter and quality trimming, with some polyA 
## version of bbduk geared towards bulk RNA-seq
## 1) everything is dumped into one directory, usually /work 
## 2) if there are multiple fastq files per sample, they will be input together 
## 3) resulting reads are not archived and should be processed faster

ADAPTERS=/software/cellgeni/bbmap/resources/adapters.fa

R1=""
R2=""
if [[ `find $FQDIR/* | grep $TAG | grep "_1\.fastq"` != "" ]]
then 
  R1=`find $FQDIR/* | grep $TAG | grep "_1\.fastq" | sort | tr '\n' ',' | sed "s/,$//g"`
  R2=`find $FQDIR/* | grep $TAG | grep "_2\.fastq" | sort | tr '\n' ',' | sed "s/,$//g"`
elif [[ `find $FQDIR/* | grep $TAG | grep "R1\.fastq"` != "" ]]
then
  R1=`find $FQDIR/* | grep $TAG | grep "R1\.fastq" | sort | tr '\n' ',' | sed "s/,$//g"`
  R2=`find $FQDIR/* | grep $TAG | grep "R2\.fastq" | sort | tr '\n' ',' | sed "s/,$//g"`
elif [[ `find $FQDIR/* | grep $TAG | grep "_R1_.*\.fastq"` != "" ]]
then
  R1=`find $FQDIR/* | grep $TAG | grep "_R1_" | sort | tr '\n' ',' | sed "s/,$//g"`
  R2=`find $FQDIR/* | grep $TAG | grep "_R2_" | sort | tr '\n' ',' | sed "s/,$//g"`
else 
  >&2 echo "ERROR: No appropriate fastq files were found! Please check file formatting, and check if you have set the right FQDIR."
  exit 1
fi 

$CMD bbduk.sh -Xmx100G in1=$R1 in2=$R2 out1=bbduk_fastqs/$TAG.bbduk.R1.fastq out2=bbduk_fastqs/$TAG.bbduk.R2.fastq ref=$ADAPTERS trimpolya=10 ktrim=r k=23 mink=11 hdist=1 tpe tbo &> $TAG.bbduk.log

## redefine FQDIR and R1/R2
FQDIR=`readlink -f bbduk_fastqs` 
R1=$FQDIR/bbduk_fastqs/$TAG.bbduk.R1.fastq
R2=$FQDIR/bbduk_fastqs/$TAG.bbduk.R2.fastq

## step 2 - align reads with STAR to genome/transcriptome

mkdir $TAG && cd $TAG

## ENCODE options for RNA-seq alignment
$CMD STAR --runThreadN $CPUS --genomeDir $SREF --readFilesIn $R1 $R2 --outFilterMultimapNmax 20 \
     --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
     --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMheaderCommentFile COfile.txt \
     --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outFilterType BySJout \
     --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --sjdbScore 1 \
     --outBAMsortingBinsN 500 --limitBAMsortRAM 60000000000 --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes All

$CMD samtools index -@$CPUS Aligned.sortedByCoord.out.bam

## step 3 - evaluate strand-specificity 

STRAND=""
FCSTR=""
UNS=`grep "^ENS" ReadsPerGene.out.tab | awk '{sum+=$2} END {printf "%d\n",sum}'`
FWD=`grep "^ENS" ReadsPerGene.out.tab | awk '{sum+=$3} END {printf "%d\n",sum}'`
REV=`grep "^ENS" ReadsPerGene.out.tab | awk '{sum+=$4} END {printf "%d\n",sum}'`

RTO1=$((FWD*100/UNS))
RTO2=$((REV*100/UNS))

echo "Strand-specificity evaluation, read counts: UNS = $UNS, FWD = $FWD, REV = $REV; FWD/UNS = $RTO1%, REV/UNS = $RTO2%"

if (( $RTO1 > 40 && $RTO1 < 60 && $RTO2 > 40 && $RTO2 < 60 )) 
then
  STRAND="none"
  FCSTR="0"
elif (( $RTO1 > 80 && $RTO2 < 20 )) 
then
  STRAND="forward"
  FCSTR="1"
elif (( $RTO2 > 80 && $RTO1 < 20 )) 
then
  STRAND="reverse"
  FCSTR="2"
else 
  >&2 echo "ERROR: Strand-specificity could not be determined! please examine the logs and metrics carefully."
fi 

echo "Strand-specificity was determined to be: STRAND = $STRAND"

## step 4 - quantify the expression using RSEM and featureCounts

$CMD rsem-calculate-expression --paired-end -p $CPUS --bam --estimate-rspd --seed 12345 --no-bam-output --strandedness $STRAND Aligned.toTranscriptome.out.bam $RREF $TAG.rsem &> $TAG.rsem.log
featureCounts -p -T $CPUS -t gene -g gene_id -s $FCSTR -a $GTF -o $TAG.feature_counts.tsv Aligned.sortedByCoord.out.bam &> $TAG.fcounts.log
