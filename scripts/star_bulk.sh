#!/bin/bash

set -euo pipefail

## v0.3
## my in-house bulk RNA-seq pipeline streamlined to work with many counter + STAR-Fusion
## bbduk.sh is used for read trimming always - check that adapter list is up to date
## esp if mapped read length is very different from actual read length
## all necessary tools are packaged in Singularity
## should work equally well for both single- and paired-end, and with .gz/.bz2/.fastq

if [[ $# -ne 2 ]]
then
  >&2 echo "Usage: ./star_bulk.sh <fastq_dir> <sample_id>"
  >&2 echo "(make sure you set the correct reference/index variables below)"
  exit 1
fi

FQDIR=$1
TAG=$2
CPUS=16
RAM=100
FQDIR=`readlink -f $FQDIR`

SIF="/nfs/cellgeni/singularity/images/star_bulk_v1.sif"
CMD="singularity run --nv --bind /nfs,/lustre,/software $SIF"
SREF="/nfs/cellgeni/STAR/human/2020A-full/index"
RREF="/nfs/cellgeni/STAR/human/2020A-full/GRCh38_v32_rsem"
SALM="/nfs/cellgeni/STAR/human/2020A-full/GRCh38_v32_salmon"
GTF="/nfs/cellgeni/STAR/human/2020A-full/GRCh38_v32_modified.gtf"
ADAPTERS="/software/cellgen/cellgeni/bbmap/resources/adapters.fa"

## step 1 - adapter and quality trimming, with some polyA trimming too  
## 1) everything is dumped into one directory, usually /work 
## 2) if there are multiple fastq files per sample, they will be input together 
## 3) resulting reads are not archived for faster processing downstream


R1=""
R2=""
if [[ `find $FQDIR/* | grep -P "\/$TAG[\/\._]" | grep "_1\.f.*q" || true` != "" ]]
then 
    R1=`find $FQDIR/* | grep -P "\/$TAG[\/\._]" | grep "_1\.f.*q" | sort | tr '\n' ' ' | sed "s/ $//g"`
    R2=`find $FQDIR/* | grep -P "\/$TAG[\/\._]" | grep "_2\.f.*q" | sort | tr '\n' ' ' | sed "s/ $//g" || true`
elif [[ `find $FQDIR/* | grep -P "\/$TAG[\/\._]" | grep "R1\.f.*q" || true` != "" ]]
then
    R1=`find $FQDIR/* | grep -P "\/$TAG[\/\._]" | grep "R1\.f.*q" | sort | tr '\n' ' ' | sed "s/ $//g"`
    R2=`find $FQDIR/* | grep -P "\/$TAG[\/\._]" | grep "R2\.f.*q" | sort | tr '\n' ' ' | sed "s/ $//g" || true`
elif [[ `find $FQDIR/* | grep -P "\/$TAG[\/\._]" | grep "_R1_.*\.f.*q" || true` != "" ]]
then
    R1=`find $FQDIR/* | grep -P "\/$TAG[\/\._]" | grep "_R1_" | sort | tr '\n' ' ' | sed "s/ $//g"`
    R2=`find $FQDIR/* | grep -P "\/$TAG[\/\._]" | grep "_R2_" | sort | tr '\n' ' ' | sed "s/ $//g" || true`
else 
    >&2 echo "ERROR: No appropriate fastq files were found! Please check file formatting, and check if you have set the right FQDIR."
    exit 1
fi
NF=`echo $R1 | grep -o " " | wc -l || true`

mkdir -p bbduk_fastqs

if [[ $R2 != "" ]]
then
    if (( $NF == 0 )) 
    then
        echo "One set of PAIRED-END fastq files found for sample $TAG. Running bbduk.sh read trimming .." 
        $CMD bbduk.sh -Xmx${RAM}G in1=$R1 in2=$R2 out1=bbduk_fastqs/$TAG.bbduk.R1.fastq out2=bbduk_fastqs/$TAG.bbduk.R2.fastq ref=$ADAPTERS trimpolya=10 ktrim=r k=23 mink=11 hdist=1 tpe tbo &> $TAG.bbduk.log
    else
        echo "Multiple sets of PAIRED-END fastq files found for sample $TAG. Running bbduk.sh read trimming sequentially .." 
        a=( $R1 ) 
        b=( $R2 ) 
        for i in `seq 0 $NF`
        do
            $CMD bbduk.sh -Xmx${RAM}G in1=${a[$i]} in2=${b[$i]} out1=bbduk_fastqs/$TAG.bbduk.$i.R1.fastq out2=bbduk_fastqs/$TAG.bbduk.$i.R2.fastq ref=$ADAPTERS trimpolya=10 ktrim=r k=23 mink=11 hdist=1 tpe tbo &> $TAG.bbduk.$i.log
        done 
    fi 
else 
    if (( $NF == 0 )) 
    then
        echo "One set of SINGLE-END fastq files found for sample $TAG. Running bbduk.sh read trimming .." 
        $CMD bbduk.sh -Xmx${RAM}G in=$R1 out=bbduk_fastqs/$TAG.bbduk.R1.fastq ref=$ADAPTERS trimpolya=10 ktrim=r k=23 mink=11 hdist=1 &> $TAG.bbduk.log
    else
        echo "Multiple sets of SINGLE-END fastq files found for sample $TAG. Running bbduk.sh read trimming sequentially .." 
        a=( $R1 ) 
        for i in `seq 0 $NF`
        do
            $CMD bbduk.sh -Xmx${RAM}G in=${a[$i]} out=bbduk_fastqs/$TAG.bbduk.$i.R1.fastq ref=$ADAPTERS trimpolya=10 ktrim=r k=23 mink=11 hdist=1 &> $TAG.bbduk.$i.log
        done 
    fi 
fi

## redefine FQDIR and R1/R2 (will work for multiple files too)
FQDIR=`readlink -f bbduk_fastqs` 
R1=`find $FQDIR/* | grep "$TAG\.bbduk" | grep "\.R1\.fastq" | sort | tr '\n' ',' | sed "s/,$//g"`
R2=`find $FQDIR/* | grep "$TAG\.bbduk" | grep "\.R2\.fastq" | sort | tr '\n' ',' | sed "s/,$//g" || true`
PAIRED="true"
if [[ $R2 == "" ]]
then
    PAIRED="false" 
fi

## step 2 - align reads with STAR to genome/transcriptome

mkdir $TAG && cd $TAG
mv ../$TAG.bbduk.log bbduk.log

## ENCODE options for RNA-seq alignment
if [[ $PAIRED == "true" ]]
then
    echo "Running STAR mapping for paired-end reads .." 
    $CMD STAR --runThreadN $CPUS --genomeDir $SREF --readFilesIn $R1 $R2 --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMheaderCommentFile COfile.txt \
        --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outFilterType BySJout \
        --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --sjdbScore 1 \
        --outBAMsortingBinsN 500 --limitBAMsortRAM ${RAM}000000000 --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes All

    $CMD samtools index -@$CPUS Aligned.sortedByCoord.out.bam
    
    ## step 3 - evaluate strand-specificity 
    
    STRAND=""
    FCSTR=""
    UNS=`grep "^ENS" ReadsPerGene.out.tab | awk '{sum+=$2} END {printf "%d\n",sum}'`
    FWD=`grep "^ENS" ReadsPerGene.out.tab | awk '{sum+=$3} END {printf "%d\n",sum}'`
    REV=`grep "^ENS" ReadsPerGene.out.tab | awk '{sum+=$4} END {printf "%d\n",sum}'`
    
    RTO1=$((FWD*100/UNS))
    RTO2=$((REV*100/UNS))
    
    echo "Strand-specificity evaluation, read counts: UNS = $UNS, FWD = $FWD, REV = $REV; FWD/UNS = $RTO1%, REV/UNS = $RTO2%" | tee -a strand.txt
    
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
    
    echo "Strand-specificity was determined to be: STRAND = $STRAND" | tee -a strand.txt
    
    ## step 4 - quantify the expression using RSEM, featureCounts, and Salmon (x2)
    echo "Running read counting using RSEM on a transcriptome BAM .." 
    $CMD rsem-calculate-expression --paired-end -p $CPUS --bam --estimate-rspd --seed 12345 --no-bam-output --strandedness $STRAND Aligned.toTranscriptome.out.bam $RREF rsem &> rsem.log
    
    ## featureCounts sorting always breaks something, so we'll have to settle for this solution: select only unique and concordantly mapped reads, and name-sort them. 
    echo "Making name-sorted paired-end BAM file and running read counting using featureCounts .." 
    $CMD samtools view -@$CPUS -F8 -d NH:1 -O BAM Aligned.sortedByCoord.out.bam | samtools sort -@$CPUS -n -O BAM - > Uniq_paired_namesorted.bam
    $CMD featureCounts -p --countReadPairs --donotsort -T $CPUS -t exon -g gene_id -s $FCSTR -a $GTF -o feature_counts.tsv Uniq_paired_namesorted.bam &> fcounts.log
    
    ## finally, use Salmon to do the same thing as RSEM (salmon_aln), and selective-map to transcriptome (salmon_reads):
    ## two lines below are because Salmon expects a space in case of multiple fastq files, while STAR/bowtie2/hisat2 etc use comma
    R1=`echo $R1 | tr ',' ' '`
    R2=`echo $R2 | tr ',' ' '`
    echo "Running read-based counting using Salmon with full decoy reference .." 
    $CMD salmon quant -p $CPUS -g $GTF -i $SALM -l A -1 $R1 -2 $R2 --validateMappings -o salmon_reads &> salmon_reads.log
    echo "Running read counting using Salmon on a transcriptome BAM .." 
    $CMD salmon quant -p $CPUS -g $GTF -t $RREF.idx.fa -l A -a Aligned.toTranscriptome.out.bam -o salmon_aln &> salmon_aln.log
else
    echo "Running STAR mapping for single-end reads .." 
    $CMD STAR --runThreadN $CPUS --genomeDir $SREF --readFilesIn $R1 --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMheaderCommentFile COfile.txt \
        --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outFilterType BySJout \
        --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --sjdbScore 1 \
        --outBAMsortingBinsN 500 --limitBAMsortRAM ${RAM}000000000 --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes All

    $CMD samtools index -@$CPUS Aligned.sortedByCoord.out.bam
    
    ## step 3 - evaluate strand-specificity 
    
    STRAND=""
    FCSTR=""
    UNS=`grep "^ENS" ReadsPerGene.out.tab | awk '{sum+=$2} END {printf "%d\n",sum}'`
    FWD=`grep "^ENS" ReadsPerGene.out.tab | awk '{sum+=$3} END {printf "%d\n",sum}'`
    REV=`grep "^ENS" ReadsPerGene.out.tab | awk '{sum+=$4} END {printf "%d\n",sum}'`
    
    RTO1=$((FWD*100/UNS))
    RTO2=$((REV*100/UNS))
    
    echo "Strand-specificity evaluation, read counts: UNS = $UNS, FWD = $FWD, REV = $REV; FWD/UNS = $RTO1%, REV/UNS = $RTO2%" | tee -a strand.txt
    
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
    
    echo "Strand-specificity was determined to be: STRAND = $STRAND" | tee -a strand.txt
    
    ## step 4 - quantify the expression using RSEM, featureCounts, and Salmon (x2)
    
    echo "Running read counting using RSEM on a transcriptome BAM (single-end) .." 
    $CMD rsem-calculate-expression -p $CPUS --bam --estimate-rspd --seed 12345 --no-bam-output --strandedness $STRAND Aligned.toTranscriptome.out.bam $RREF rsem &> rsem.log
    
    ## featureCounts sorting always breaks something, so we'll have to settle for this solution: select only unique and concordantly mapped reads, and name-sort them. 
    echo "Running single-end read counting using featureCounts .." 
    $CMD featureCounts -T $CPUS -t exon -g gene_id -s $FCSTR -a $GTF -o feature_counts.tsv Aligned.sortedByCoord.out.bam &> fcounts.log
    
    ## finally, use Salmon to do the same thing as RSEM (salmon_aln), and selective-map to transcriptome (salmon_reads):
    ## two lines below are because Salmon expects a space in case of multiple fastq files, while STAR/bowtie2/hisat2 etc use comma
    R1=`echo $R1 | tr ',' ' '`
    echo "Running read-based counting using Salmon with full decoy reference (single-end) .." 
    $CMD salmon quant -p $CPUS -g $GTF -i $SALM -l A -r $R1 --validateMappings -o salmon_reads &> salmon_reads.log
    echo "Running read counting using Salmon on a transcriptome BAM (single-end) .." 
    $CMD salmon quant -p $CPUS -g $GTF -t $RREF.idx.fa -l A -a Aligned.toTranscriptome.out.bam -o salmon_aln &> salmon_aln.log
fi 

echo "ALL DONE!"
