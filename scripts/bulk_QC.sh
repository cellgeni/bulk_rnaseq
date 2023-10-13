#!/bin/bash 

## check that STAR temporary dir is removed for all samples, and that archived unmapped reads are created
>&2 echo "Checking that all STAR jobs went to completion .." 
for i in *
do
  if [[ -d $i && -s $i/Log.out ]]
  then
    if [[ -d $i/_STARtmp ]]
    then
      >&2 echo "WARNING: Sample $i did not run to completion: _STARtmp is still present!" 
    fi
  fi 
done

echo -e "Sample\tRead_len\tMapped_len\tN_reads\tN_uniq\tN_multi\tPct_uniq\tPct_M1\tPct_M2\tPct_U1\tPct_U2\tPct_U3\tSTAR_U\tFcounts\tRSEM\tSalmon_aln\tSalmon_reads"

for i in *
do
  if [[ -d $i && -s $i/Log.out ]]
  then 
    RL=`grep "Average input read length " $i/Log.final.out | awk '{print $6}'`	
    ML=`grep "Average mapped length " $i/Log.final.out | awk '{print $5}'`	
    N1=`grep "Number of input reads " $i/Log.final.out | awk '{print $6}'`	
    N2=`grep "Uniquely mapped reads number " $i/Log.final.out | awk '{print $6}'`
    N3=`grep "Number of reads mapped to multiple loci " $i/Log.final.out | awk '{print $9}'`
    
    P1=`grep "Uniquely mapped reads %" $i/Log.final.out | awk -F "|" '{print $2}' | awk '{print $1}'`	
    P2=`grep "% of reads mapped to multiple loci" $i/Log.final.out | awk -F "|" '{print $2}' | awk '{print $1}'`	
    P3=`grep "% of reads mapped to too many loci" $i/Log.final.out | awk -F "|" '{print $2}' | awk '{print $1}'`	
    P4=`grep "% of reads unmapped: too many mismatches" $i/Log.final.out | awk -F "|" '{print $2}' | awk '{print $1}'`	
    P5=`grep "% of reads unmapped: too short" $i/Log.final.out | awk -F "|" '{print $2}' | awk '{print $1}'`	
    P6=`grep "% of reads unmapped: other" $i/Log.final.out | awk -F "|" '{print $2}' | awk '{print $1}'`	
    
    NS=`grep "^ENS" $i/ReadsPerGene.out.tab | awk '{sum+=$2} END {printf "%d\n",sum}'`
    NF=`grep "^ENS" $i/$i.feature_counts.tsv | awk '{sum+=$7} END {printf "%d\n",sum}'`
    NR=`grep "^ENS" $i/$i.rsem.genes.results | awk '{sum+=$5} END {printf "%d\n",sum}'`
    SA=`grep "^ENS" $i/$i.salmon_aln/quant.genes.sf | awk '{sum+=$5} END {printf "%d\n",sum}'`
    SR=`grep "^ENS" $i/$i.salmon_reads/quant.genes.sf | awk '{sum+=$5} END {printf "%d\n",sum}'`
    echo -e "$i\t$RL\t$ML\t$N1\t$N2\t$N3\t$P1\t$P2\t$P3\t$P4\t$P5\t$P6\t$NS\t$NF\t$NR\t$SA\t$SR"
  fi
done 
