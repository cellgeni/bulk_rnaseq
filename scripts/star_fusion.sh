#!/bin/bash 

## all Singularity, no bs

SIF="/nfs/cellgeni/singularity/images/starfusion-1.12.sif"
CMD="singularity run --nv --bind /nfs,/lustre,/software $SIF"

FQDIR=$1
TAG=$2
CPUS=16

FQDIR=`readlink -f $FQDIR`
FREF="/nfs/cellgeni/STAR-Fusion/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir"

mkdir ${TAG}_fusion && cd ${TAG}_fusion

R1=$FQDIR/bbduk_fastqs/$TAG.bbduk.R1.fastq
R2=$FQDIR/bbduk_fastqs/$TAG.bbduk.R2.fastq

$CMD STAR-Fusion --left_fq $R1 --right_fq $R2 --genome_lib_dir $FREF --FusionInspector validate --denovo_reconstruct --examine_coding_effect
