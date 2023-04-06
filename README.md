# Wrapper scripts for bulk RNA-seq pipeline

These are the updated scripts used for CellGenIT for uniform processing of bulk RNA-seq. All listed methods use [STAR](https://github.com/alexdobin/STAR) aligner to align reads to the reference genome and transcriptome. Following the alignment, quantification is done using 3 different options: 1) built-in STAR read counter on genomic BAM; 2) featureCounts (see below for exact options) on genomic BAM; 3) RSEM on transcriptomic BAM for proper EM-based multimapper counting; 4) Salmon in alignment mode with full genome decoy, for EM-based multimapper counting and perhaps slightly fewer annoyances related to STAR mapping to transcriptome and RSEM counting. 

Additionally, [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion) could be used to find and annotate fusion transcripts, usually from cancer samples.  

## Software installation

### Tool versions

`STAR` of version 2.7.9a or above is recommended (2.7.10a is the latest and greatest, as of August'22). RSEM, featureCounts, and Salmon are now stable and are not frequently updated. BBduk is used for adapter and quality trimming. 

In case `STAR-Fusion` is required, this should be done using a separate run (starting from bbduk-trimmed fastq files), and using a separate `singularity` image. 

### Docker/Singularity containers

Everything CellGenIT runs routinely is done with `singularity` images for reproducibility. 
 
Using Dockerfile which builds with specific versions of all the tools used. Once the container is built you can do `cat /versions.txt` and it will show the versions of the tools in the container. Scripts provided in `/scripts` are Farm-specific; singularity images referred in the scripts should be available to all Sanger users. 

`STAR-Fusion` image was built using this command (latest Docker here is digest `9bf8e66c0c86`, `STAR-Fustion` v1.12): 

`singularity build starfusion-1.12.sif docker://trinityctat/starfusion:latest` 

## Reference genome and annotation

Unlike `STARsolo` processing, for which we often exactly mimic the reference provided by `Cell Ranger`, bulk RNA-seq might require many different versions of reference genome and annotation. For the purpose of scRNA-seq to bulk RNA-seq mapping, it's best to use the exactly matching reference (e.g., `/nfs/cellgeni/STAR/human/2020A` or `/nfs/cellgeni/STAR/human/2020A-full`). When `STAR-Fusion` analysis is required, it's convenient 

Once you have downloaded the needed fasta and GTF files, and applied the necessary filtering (see the 10x genomics link above), run the following command using 16 CPUs/64 Gbs of RAM: 

`STAR --runThreadN 16 --runMode genomeGenerate --genomeDir STAR --genomeFastaFiles $FA --sjdbGTFfile $GTF`

`RSEM` reference is made from the GTF file and genome fasta; e.g. for the `2020A-full` reference, GRCh38_v32_modified.fa (genome) and GRCh38_v32_modified.gtf (GTF), the command would look like this:

`rsem-prepare-reference --gtf GRCh38_v32_modified.gtf GRCh38_v32_modified.fa Gencode_v32_rsem`

`Salmon` full genome decay reference is made as follows (last step is rather resource-demanding, run it on 16 cores & 128 Gb RAM): 

```bash
grep "^>" GRCh38_v32_modified.fa | sed 's/>//g' > decoys.txt
cat GRCh38_v32_rsem.idx.fa GRCh38_v32_modified.fa | gzip > GRCh38_v32_gentrome.fa.gz
salmon index -t GRCh38_v32_gentrome.fa.gz -d decoys.txt -p 16 -i GRCh38_v32_salmon -k 25 --keepDuplicates --no-clip
```

For `STAR-Fusion`, [this link](https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/) can be used to obtain the pre-made reference files. Currently, Gencode v19 is available for genome assembly GRCh37, and Gencode v22 and v37 are available for GRCh38. There are also some files for T2T already, which is kewl. 

### Pre-made reference location

All **CellGenIT** pre-made `STAR` references are located in `/nfs/cellgeni/STAR/`. `Salmon`, `RSEM`, and other indexes can usually be found in the same directory. 
 
## Resource requirements 

By default, all processing is done using 16 CPUs and 128 Gb of RAM. Using the latest settings, you _should_ be able to process most 10x experiments with 64 Gb of RAM even when you need a sorted BAM output. Farm typically has ~8 Gb of RAM per core, so 8 CPUs/64 Gb RAM, or 16 CPUs/128 Gb RAM is probably optimal. 

### Mapping reads with STAR using ENCODE presets

By default, all reads are mapped to the genome and transcriptome using ENCODE options (described in details in STAR manual). 

Note that transcriptome mapping does not allow soft clipping, which causes severe problems if your reads have untrimmed adaptors or UMIs. The output of `bulk_QC.sh` should be informative as to whether transcriptome mapping was successful.   

### Counting the multimapping reads

## Quick evaluation of multiple bulk RNA-seq runs

If you've used these scripts to process multiple RNA-seq samples, you can get a quick look at the results by copying `bulk_QC.sh` script from this repo to the directory with `STAR` output folders, and running

```bash
./bulk_QC.sh | column -t 
```
The script will provide a rough overview of mapping and read counting statistics. 
