#!/bin/bash 

set -euo pipefail
## all Singularity, no bs
## paired-end only for now I think 
## this is normally run alongside star_bulk.sh so you'd have your normal expression tables too 

SIF="/nfs/cellgeni/singularity/images/starfusion-1.12.sif"
CMD="singularity run --nv --bind /nfs,/lustre,/software $SIF"

FQDIR=$1
TAG=$2
CPUS=16

FQDIR=`readlink -f $FQDIR`
FREF="/nfs/cellgeni/refdata_starfusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir" 

mkdir -p ${TAG}_fusion && cd ${TAG}_fusion

R1=$FQDIR/bbduk_fastqs/$TAG.bbduk.R1.fastq
R2=$FQDIR/bbduk_fastqs/$TAG.bbduk.R2.fastq

if [[ ! -s $R1 || ! -s $R2 ]]
then
    >&2 echo "ERROR: read 1 ($R1) or read 2 ($R2) cannot be found!" 
    exit 1
fi

$CMD STAR-Fusion --left_fq $R1 --right_fq $R2 --tmpdir `pwd` --CPU $CPUS --genome_lib_dir $FREF --FusionInspector validate --denovo_reconstruct --examine_coding_effect &> $TAG.star-fusion.log 
