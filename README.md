# Wrapper scripts for bulk RNA-seq pipeline

These are the updated scripts used for CellGenIT for uniform processing of bulk RNA-seq. All listed methods but `Salmon` use [STAR](https://github.com/alexdobin/STAR) aligner to align reads to the reference genome and/or transcriptome. Following the alignment, quantification is done using the following options: 

  - built-in STAR read counter on genomic BAM;
  - featureCounts on genomic BAM;
  - RSEM on transcriptomic BAM for proper EM-based multimapper counting; 
  - Salmon in _read mapping_ mode with full genome decoy, for [general awesomeness](https://salmon.readthedocs.io/en/latest/salmon.html);
  - Salmon in _alignment_ mode for debugging and RSEM sanity check.  

Additionally, [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion) could be used to find and annotate fusion transcripts. Usually this is done on cancer samples.  

## Software installation

### Tool versions

Current tool versions are listed below: 

```
star_version=2.7.10a_alpha_220818
samtools_version=1.15.1
bbmap_version=38.97
rsem_version=1.3.3
subread_version=2.0.2
salmon_version=1.10.0
```

In case `STAR-Fusion` is required, this should be done using a separate run (starting from `BBduk`-trimmed fastq files), and using a separate `Singularity` image (see below). Currently `STAR-Fusion` v1.12 is used in this pipeline. 

### Docker/Singularity containers

Everything CellGenIT runs routinely is done with `Singularity` images for reproducibility. This pipepline uses 2 images: one for general bulk RNA-seq, and another one for `STAR-Fusion`. Provided _Dockerfile_ is used for bulk RNA-seq processing. 
`STAR-Fusion` image was built using this command (latest Docker here is digest `9bf8e66c0c86`, `STAR-Fustion` v1.12): 

`singularity build starfusion-1.12.sif docker://trinityctat/starfusion:latest` 

Scripts provided in `/scripts` are Farm-specific; singularity images referred in the scripts should be available to all Sanger users. 

## Reference genome and annotation

Unlike `STARsolo` processing, for which we often exactly mimic the reference provided by `Cell Ranger`, bulk RNA-seq might require many different versions of reference genome and annotation. For the purpose of scRNA-seq to bulk RNA-seq mapping, it's best to use the exactly matching reference (e.g., `/nfs/cellgeni/STAR/human/2020A` or `/nfs/cellgeni/STAR/human/2020A-full`). When `STAR-Fusion` analysis is required, it's convenient to use pre-made reference bundles, which currently are available for Gencode v19, v22, and v37. 

### Genome and GTF source and pre-processing 

Normally, assembly and GTF files provided by [Gencode]() are used to create reference indexes. The following things need to be noted about the reference genome and GTF preparation:

  - Whenever available, [primary assembly](https://www.ncbi.nlm.nih.gov/grc/help/definitions/) and a corresponding GTF file is used;
  - Genome and GTF files should be pre-processed to get rid of ENSEMBL gene/transcript versions in the downstream results, similarly to how it's done in 10x reference preparation instructions; 
  - The second copy of pseudo-autosomal genes should be removed from the GTF. 
  
### Index generation for `STAR`, `RSEM`, and `Salmon`

`STAR` index is generated using the following command (you'll need about 64G of RAM): 

`STAR --runThreadN 16 --runMode genomeGenerate --genomeDir index --genomeFastaFiles $FA --sjdbGTFfile $GTF`

`RSEM` transcriptome reference is made from the GTF file and genome fasta; e.g. for the `2020A-full` reference, GRCh38_v32_modified.fa (genome) and GRCh38_v32_modified.gtf (GTF), the command would look like this:

`rsem-prepare-reference --gtf GRCh38_v32_modified.gtf GRCh38_v32_modified.fa Gencode_v32_rsem`

`Salmon` full genome decay reference is made as follows (last step is rather resource-demanding, run it on 16 cores & 128 Gb RAM): 

```bash
grep "^>" GRCh38_v32_modified.fa | sed 's/>//g' > decoys.txt
cat GRCh38_v32_rsem.idx.fa GRCh38_v32_modified.fa | gzip > GRCh38_v32_gentrome.fa.gz
salmon index -t GRCh38_v32_gentrome.fa.gz -d decoys.txt -p 16 -i GRCh38_v32_salmon -k 25 --keepDuplicates --no-clip
```

For `STAR-Fusion`, [this link](https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/) can be used to obtain the pre-made reference files. Currently, Gencode v19 is available for genome assembly GRCh37, and Gencode v22 and v37 are available for GRCh38. There are also some files for T2T already, which is kewl. 

### Reference location on Farm

All **CellGenIT** pre-made `STAR` references are located in `/nfs/cellgeni/STAR/`. `Salmon`, `RSEM`, and other indexes can usually be found in the same directory. 
 
## Resource requirements 

By default, all processing is done using 16 CPUs and 128 Gb of RAM. Using the latest settings, you _should_ be able to process most 10x experiments with 64 Gb of RAM even when you need a sorted BAM output. Farm typically has ~8 Gb of RAM per core, so 8 CPUs/64 Gb RAM, or 16 CPUs/128 Gb RAM is probably optimal. 

## Mapping reads with `STAR` using ENCODE presets

By default, all reads are mapped to the genome and transcriptome using ENCODE options (described in details in STAR manual). Note that transcriptome mapping does not allow soft clipping, which causes severe problems if your reads have untrimmed adaptors or UMIs. RSEM also does not allow any indels to be present in the (transcriptomic) alignment. The output of `bulk_QC.sh` should be informative as to whether transcriptome mapping was successful.

`STAR` quantifies the reads using the 3 strand-specificity presets (unstranded, forward, reverse); these are available in the file named `ReadsPerGene.out.tab`. Our scripts use the output to evaluate the strand-specificity of each sample. 

### Counting the multimapping reads

While `STAR` and `featureCounts` do not count multimapping reads by default, `RSEM` with a transcriptomic mapping is a well-established solution that fills this gap. `Salmon` does so as well, but using a somewhat different Bayesian algorithm and with other improvements.

### Using `Salmon` with full genome decoy

`Salmon` can use _selective alignment_ algorithm, in which whole genome is basically used a decoy to avoid spurious mapping of intronic reads to transcripts, often observed with classical pseudo-mapping. [This paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8) contains all the necessary details. All in all, `Salmon` in this mode should probably be considered the most reliable and accurate estimate of expression. 

## Quick evaluation of multiple bulk RNA-seq runs

If you've used these scripts to process multiple RNA-seq samples, you can get a quick look at the results by copying `bulk_QC.sh` script from this repo to the directory with `STAR` output folders, and running

```bash
./bulk_QC.sh | column -t 
```
The script will provide a rough overview of mapping and read counting statistics by each tool described above. 
