#!/bin/bash 

## all Singularity, no bs

SIF="/nfs/cellgeni/singularity/images/starsolo_2-7-10a-alpha-220818_samtools_1-15-1_seqtk-1-13_bbmap_38-97_RSEM-1-3-3.sif"
CMD="singularity run --nv --bind /nfs,/lustre,/software $SIF"

FQDIR=$1
TAG=$2
CPUS=16

FQDIR=`readlink -f $FQDIR`
SREF=/nfs/cellgeni/STAR/human/Gencode_v19_full/index
RREF=/nfs/cellgeni/STAR/human/Gencode_v19_full/Gencode_v19_rsem 
GTF=/nfs/cellgeni/STAR/human/Gencode_v19_full/gencode.v19.chr_patch_hapl_scaff.annotation.gtf

## for single-end reads, change both STAR and RSEM options
## for strand-specific processing, change --forward-prob in rsem-calculate-expression (1 is FR, 0 is RF, 0.5 is non-strand-specific). 

## assume we ran bbduk on thangs before, so really there's only 1 (or 2 for PE) fastq file that's not gzipped 
mkdir $TAG && cd $TAG

R1=$FQDIR/bbduk_fastqs/$TAG.bbduk.R1.fastq
R2=$FQDIR/bbduk_fastqs/$TAG.bbduk.R2.fastq

## ENCODE options for RNA-seq alignment
$CMD STAR --runThreadN $CPUS --genomeDir $SREF --readFilesIn $R1 $R2 --outFilterMultimapNmax 20 \
     --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
     --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMheaderCommentFile COfile.txt \
     --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outFilterType BySJout \
     --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --sjdbScore 1 \
     --outBAMsortingBinsN 500 --limitBAMsortRAM 60000000000 --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes All

$CMD samtools index -@$CPUS Aligned.sortedByCoord.out.bam

STRAND=""
FCSTR=""
UNS=`grep "^ENS" ReadsPerGene.out.tab | awk '{sum+=$2} END {printf "%d\n",sum}'`
FWD=`grep "^ENS" ReadsPerGene.out.tab | awk '{sum+=$3} END {printf "%d\n",sum}'`
REV=`grep "^ENS" ReadsPerGene.out.tab | awk '{sum+=$4} END {printf "%d\n",sum}'`

RTO1=$((FWD*100/UNS))
RTO2=$((REV*100/UNS))

echo "Strand-specificity evaluation, read counts: UNS = $UNS, FWD = $FWD, REV = $REV; FWD/UNS = $RTO1 %, REV/UNS = $RTO2 %"

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

$CMD rsem-calculate-expression --paired-end -p $CPUS --bam --estimate-rspd --seed 12345 --no-bam-output --strandedness $STRAND Aligned.toTranscriptome.out.bam $RREF $TAG.rsem &> $TAG.rsem.log
featureCounts -p -T $CPUS -t gene -g gene_id -s $FCSTR -a $GTF -o $TAG.feature_counts.tsv Aligned.sortedByCoord.out.bam &> $TAG.fcounts.log
