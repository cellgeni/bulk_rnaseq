#!/bin/bash 

## all-in-one processing script that should work as individual submission as well as part of job array
## all necessary tools are packaged in Singularity
## everything here is for paired-end reads, single end would require many option changes 

SIF="/nfs/cellgeni/singularity/images/star_bulk_v1.sif"
CMD="singularity run --nv --bind /nfs,/lustre,/software $SIF"

FQDIR=$1
TAG=$2
CPUS=16

if [[ $FQDIR == "" || $TAG == "" ]]
then
  >&2 echo "Usage: ./star_bulk.sh <fastq_dir> <sample_id>"
  >&2 echo "(make sure you set the correct reference/index variables below)"
  exit 1
fi

FQDIR=`readlink -f $FQDIR`
SREF=/nfs/cellgeni/STAR/human/2020A-full/index
RREF=/nfs/cellgeni/STAR/human/2020A-full/GRCh38_v32_rsem
SALM=/nfs/cellgeni/STAR/human/2020A-full/GRCh38_v32_salmon
GTF=/nfs/cellgeni/STAR/human/2020A-full/GRCh38_v32_modified.gtf
ADAPTERS=/software/cellgeni/bbmap/resources/adapters.fa

## step 1 - adapter and quality trimming, with some polyA 
## version of bbduk geared towards bulk RNA-seq
## 1) everything is dumped into one directory, usually /work 
## 2) if there are multiple fastq files per sample, they will be input together 
## 3) resulting reads are not archived and should be processed faster


R1=""
R2=""
if [[ `find $FQDIR/* | grep -P "$TAG[\/\._]" | grep "_1\.f.*q"` != "" ]]
then 
  R1=`find $FQDIR/* | grep -P "$TAG[\/\._]" | grep "_1\.f.*q" | sort | tr '\n' ' ' | sed "s/ $//g"`
  R2=`find $FQDIR/* | grep -P "$TAG[\/\._]" | grep "_2\.f.*q" | sort | tr '\n' ' ' | sed "s/ $//g"`
elif [[ `find $FQDIR/* | grep -P "$TAG[\/\._]" | grep "R1\.f.*q"` != "" ]]
then
  R1=`find $FQDIR/* | grep -P "$TAG[\/\._]" | grep "R1\.f.*q" | sort | tr '\n' ' ' | sed "s/ $//g"`
  R2=`find $FQDIR/* | grep -P "$TAG[\/\._]" | grep "R2\.f.*q" | sort | tr '\n' ' ' | sed "s/ $//g"`
elif [[ `find $FQDIR/* | grep -P "$TAG[\/\._]" | grep "_R1_.*\.f.*q"` != "" ]]
then
  R1=`find $FQDIR/* | grep -P "$TAG[\/\._]" | grep "_R1_" | sort | tr '\n' ' ' | sed "s/ $//g"`
  R2=`find $FQDIR/* | grep -P "$TAG[\/\._]" | grep "_R2_" | sort | tr '\n' ' ' | sed "s/ $//g"`
else 
  >&2 echo "ERROR: No appropriate fastq files were found! Please check file formatting, and check if you have set the right FQDIR."
  exit 1
fi
NF=`echo $R1 | grep -o " " | wc -l`

if [[ ! -d bbduk_fastqs ]]
then 
  mkdir bbduk_fastqs
fi  

if (( $NF == 0 )) 
then 
  $CMD bbduk.sh -Xmx100G in1=$R1 in2=$R2 out1=bbduk_fastqs/$TAG.bbduk.R1.fastq out2=bbduk_fastqs/$TAG.bbduk.R2.fastq ref=$ADAPTERS trimpolya=10 ktrim=r k=23 mink=11 hdist=1 tpe tbo &> $TAG.bbduk.log
else
  >&2 echo "WARNING: Multiple fastq files found for sample $TAG. Processing them sequentially.." 
  a=( $R1 ) 
  b=( $R2 ) 
  for i in `seq 0 $NF`
  do
    $CMD bbduk.sh -Xmx100G in1=${a[$i]} in2=${b[$i]} out1=bbduk_fastqs/$TAG.bbduk.$i.R1.fastq out2=bbduk_fastqs/$TAG.bbduk.$i.R2.fastq ref=$ADAPTERS trimpolya=10 ktrim=r k=23 mink=11 hdist=1 tpe tbo &> $TAG.bbduk.$i.log
  done 
fi 


## redefine FQDIR and R1/R2 (will work for multiple files too)
FQDIR=`readlink -f bbduk_fastqs` 
R1=`find $FQDIR/* | grep "$TAG\.bbduk" | grep "\.R1\.fastq" | sort | tr '\n' ',' | sed "s/,$//g"`
R2=`find $FQDIR/* | grep "$TAG\.bbduk" | grep "\.R2\.fastq" | sort | tr '\n' ',' | sed "s/,$//g"`

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

echo "Strand-specificity evaluation, read counts: UNS = $UNS, FWD = $FWD, REV = $REV; FWD/UNS = $RTO1%, REV/UNS = $RTO2%" > strand_info.txt

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
  >&2 echo "WARNING: Strand-specificity could not be determined! please examine the logs and metrics carefully."
  STRAND="none"
  FCSTR="0"
fi 

echo "Strand-specificity was determined to be: STRAND = $STRAND" >> strand_info.txt

## step 4 - quantify the expression using RSEM, featureCounts, and Salmon (x2)

$CMD rsem-calculate-expression --paired-end -p $CPUS --bam --estimate-rspd --seed 12345 --no-bam-output --strandedness $STRAND Aligned.toTranscriptome.out.bam $RREF $TAG.rsem &> $TAG.rsem.log

## featureCounts sorting always breaks something, so we'll have to settle for this solution: select only unique and concordantly mapped reads, and name-sort them. 
$CMD samtools view -@$CPUS -F8 -d NH:1 -O BAM Aligned.sortedByCoord.out.bam | samtools sort -@$CPUS -n -O BAM - > Uniq_paired_namesorted.bam
$CMD featureCounts -p --countReadPairs --donotsort -T $CPUS -t exon -g gene_id -s $FCSTR -a $GTF -o $TAG.feature_counts.tsv Uniq_paired_namesorted.bam &> $TAG.fcounts.log

## finally, use Salmon to do the same thing as RSEM (salmon_aln), and selective-map to transcriptome (salmon_reads):
## two lines below are because Salmon expects a space in case of multiple fastq files, while STAR/bowtie2/hisat2 etc use comma
R1=`echo $R1 | tr ',' ' '`
R2=`echo $R2 | tr ',' ' '`
$CMD salmon quant -p $CPUS -g $GTF -i $SALM -l A -1 $R1 -2 $R2 --validateMappings -o $TAG.salmon_reads &> $TAG.salmon_reads.log
$CMD salmon quant -p $CPUS -g $GTF -t $RREF.idx.fa -l A -a Aligned.toTranscriptome.out.bam -o $TAG.salmon_aln &> $TAG.salmon_aln.log
