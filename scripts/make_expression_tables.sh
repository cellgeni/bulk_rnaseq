#!/bin/bash 

## v2 - with gene names
## also, let's not use horrific gene lengths given by featureCounts 

## assume we have a featureCounts output done against the same GTF file 
echo "Finding the appropriate GTF file and extracting ENSG-to-name relationships.."
FCTSV=`find * | grep "\.feature_counts.tsv$" | head -n1`
GTF=`head -n1 $FCTSV | tr '\"' '\n' | grep "gtf$"`
echo "GTF file found: $GTF. Continuing.. "
perl -ne 'print "$1\t$2\n" if (m/\tgene\t.*gene_id \"(.*?)\";.*gene_name \"(.*?)\";/)' $GTF | sort -k1,1 > ensg_to_name.tsv
N=`cat ensg_to_name.tsv | wc -l`
echo "$N genes and their names extracted. Continuing.." 
echo 

for i in *
do
  if [[ -d $i && -s $i/Log.final.out ]]
  then
    echo "Processing sample $i.."
    if [[ ! -f names.tmp ]]
    then
      ## run a quick failsafe thing 
      NR=`grep -c "^ENS" $i/$i.rsem.genes.results`
      NS=`grep -c "^ENS" $i/$i.salmon_reads/quant.genes.sf`
      if (( $N != $NR || $N != $NS ))
      then
        >&2 echo "ERROR: Number of genes in RSEM and/or Salmon output does not match the number of genes in the GTF file ($GTF)! Please investigate."
        exit 1
      fi 
    
      ## make a table of 4 cols: ENSEMBL ID, Gene symbol, gene length, and gene effective length. 
      echo -e "ensemblID\tgeneSymbol\tgeneLength\teffLength" > names.tmp
      grep "^ENS" $i/$i.rsem.genes.results | sort -k1,1 | cut -f3,4 >> sorted_lengths.tmp
      paste ensg_to_name.tsv sorted_lengths.tmp >> names.tmp
      rm sorted_lengths.tmp
    fi 
    echo $i > $i.fc.tmp
    echo $i > $i.r1.tmp
    echo $i > $i.r2.tmp
    echo $i > $i.r3.tmp
    echo $i > $i.sa1.tmp
    echo $i > $i.sa2.tmp
    echo $i > $i.sr1.tmp
    echo $i > $i.sr2.tmp
    grep "^ENS" $i/$i.feature_counts.tsv | sort -k1,1 | cut -f7 >> $i.fc.tmp
    grep "^ENS" $i/$i.rsem.genes.results | sort -k1,1 | cut -f5 >> $i.r1.tmp
    grep "^ENS" $i/$i.rsem.genes.results | sort -k1,1 | cut -f6 >> $i.r2.tmp
    grep "^ENS" $i/$i.rsem.genes.results | sort -k1,1 | cut -f7 >> $i.r3.tmp
    grep "^ENS" $i/$i.salmon_aln/quant.genes.sf | sort -k1,1 | cut -f5 >> $i.sa1.tmp
    grep "^ENS" $i/$i.salmon_aln/quant.genes.sf | sort -k1,1 | cut -f4 >> $i.sa2.tmp
    grep "^ENS" $i/$i.salmon_reads/quant.genes.sf | sort -k1,1 | cut -f5 >> $i.sr1.tmp
    grep "^ENS" $i/$i.salmon_reads/quant.genes.sf | sort -k1,1 | cut -f4 >> $i.sr2.tmp
  fi
done 

paste names.tmp *.fc.tmp  > featureCounts.counts.tsv
paste names.tmp *.r1.tmp  > rsem.counts.tsv
paste names.tmp *.r2.tmp  > rsem.TPM.tsv
paste names.tmp *.r3.tmp  > rsem.FPKM.tsv
paste names.tmp *.sa1.tmp > salmon_aln.counts.tsv
paste names.tmp *.sa2.tmp > salmon_aln.TPM.tsv
paste names.tmp *.sr1.tmp > salmon_reads.counts.tsv
paste names.tmp *.sr2.tmp > salmon_reads.TPM.tsv
rm names.tmp *.fc.tmp *.r?.tmp *.s??.tmp 

echo "ALL DONE!" 
  
