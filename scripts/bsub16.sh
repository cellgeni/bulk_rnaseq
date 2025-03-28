#!/bin/bash

SCRIPT=$1
ARG=${@:2}
AA=( "$@" )
## assume the last passed argument is a sample ID or something equivalent
TAG=${AA[-1]}  
## in case the last argument is a path, we only want the sample ID
TAG=`basename $TAG` 

GROUP=`id -gn`
CPUS=16
RAM=128000 
QUE="normal"
WDIR=`pwd`

#################

bsub -G $GROUP -n$CPUS -R"span[hosts=1] select[mem>$RAM] rusage[mem=$RAM]" -M$RAM -o $WDIR/$TAG.%J.bsub.log -e $WDIR/$TAG.%J.bsub.err -q $QUE $WDIR/$SCRIPT $ARG
