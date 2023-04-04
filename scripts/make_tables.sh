#!/bin/bash 

## v2 - with gene names
## also, let's not use horrific gene lengths given by featureCounts 

GTF=/nfs/cellgeni/STAR/human/Gencode_v19_full/gencode.v19.chr_patch_hapl_scaff.annotation.gtf
perl -ne 'print "$1\t$2\n" if (m/\tgene\t.*gene_id \"(.*?)\";.*gene_name \"(.*?)\";/)' $GTF | sort -k1,1 > ensg_to_name.tsv

for i in *
do
  if [[ -d $i && -f $i/$i.feature_counts.tsv ]]
  then
    echo "Processing sample $i.."
    if [[ ! -f 000_names.tmp ]]
    then
      echo -e "ensemblID\tgeneSymbol\tgeneLength\teffLength" > 000_names.tmp
      awk 'NR>1' $i/$i.rsem.genes.results | sort -k1,1 | cut -f3,4 >> sorted_lengths.tmp
      paste ensg_to_name.tsv sorted_lengths.tmp >> 000_names.tmp
      rm sorted_lengths.tmp
    fi 
    echo $i > $i.fc.tmp
    echo $i > $i.r1.tmp
    echo $i > $i.r2.tmp
    echo $i > $i.r3.tmp
    awk 'NR>2' $i/$i.feature_counts.tsv | sort -k1,1 | cut -f7 >> $i.fc.tmp
    awk 'NR>1' $i/$i.rsem.genes.results | sort -k1,1 | cut -f5 >> $i.r1.tmp
    awk 'NR>1' $i/$i.rsem.genes.results | sort -k1,1 | cut -f6 >> $i.r2.tmp
    awk 'NR>1' $i/$i.rsem.genes.results | sort -k1,1 | cut -f7 >> $i.r3.tmp
  fi
done 

paste 000_names.tmp *.fc.tmp > featureCounts.counts.tsv
paste 000_names.tmp *.r1.tmp > rsem.counts.tsv
paste 000_names.tmp *.r2.tmp > rsem.TPM.tsv
paste 000_names.tmp *.r3.tmp > rsem.FPKM.tsv
#rm *tmp

echo "ALL DONE!" 
  
