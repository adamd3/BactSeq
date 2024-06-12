# BactSeq

![BactSeq](https://github.com/adamd3/BactSeq/actions/workflows/ci.yml/badge.svg)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**BactSeq** is a Nextflow pipeline for performing bacterial RNA-Seq analysis.

## News and updates

- 12/06/24: Dockerfile updated to install environment from yml. There were issues with installing a couple of R packages in the previous version.
- 06/10/23: Docker and Singularity support now working again. The pipeline will automatically pull the image from Docker Hub with default settings.
- 28/09/23: Please see example contrasts table and functional enrichment file below in README.

## Pipeline summary

The pipeline will perform the following steps:

1. Trim adaptors from reads ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Align reads to reference genome ([`BWA-MEM`](https://github.com/lh3/bwa/))
4. Size-factor scaling and gene length (RPKM) scaling of counts (Using median-of-ratios method from [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and TMM from [`edgeR`](http://bioconductor.org/packages/release/bioc/html/edgeR.html))
5. Principal component analysis (PCA) of normalised expression values
6. Differential gene expression ([`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)) (optional)
7. Functional enrichment of differentially expressed genes ([`topGO`](https://bioconductor.org/packages/release/bioc/html/topGO.html)) (optional)

## Installation

You will need to install [`Nextflow`](https://www.nextflow.io/) (version 21.10.3+).

```
Usage:
nextflow run BactSeq --data_dir [dir] --sample_file [file] --ref_genome [file] --ref_ann [file] -profile conda [other_options]

Mandatory arguments:
  --data_dir [file]               Path to directory containing FastQ files.
  --ref_genome [file]             Path to FASTA file containing reference genome sequence (bwa) or multi-FASTA file containing coding gene sequences (kallisto).
  --sample_file [file]            Path to file containing sample information.
  -profile [str]                  Configuration profile to use.
                                  Available: conda, docker, singularity.

Other options:
  --aligner [str]                 (Pseudo-)aligner to be used. Options: `bwa`, `kallisto`. Default = bwa.
  --ref_ann [file]                Path to GFF file containing reference genome annotation. Required only if bwa aligner used.
  --cont_tabl [file]              Path to tsv file containing contrasts to be performed for differential expression.
  --fragment_len [str]            Estimated average fragment length for kallisto transcript quantification (only required for single-end reads). Default = 150.
  --fragment_sd [str]             Estimated standard deviation of fragment length for kallisto transcript quantification (only required for single-end reads). Default = 20.
  --func_file [file]              Path to GFF3-format file containing functional annotations.
  --l2fc_thresh [str]             Absolute log2(FoldChange) threshold for identifying differentially expressed genes. Default = 1.
  --outdir [file]                 The output directory where the results will be saved (Default: './results').
  --paired [str]                  Data are paired-end.
  --p_thresh [str]                Adjusted p-value threshold for identifying differentially expressed genes. Default = 0.05.
  --skip_trimming [bool]          Do not trim adaptors from FastQ files.
  --strandedness [str]            Is data stranded? Options: `unstranded`, `forward`, `reverse`. Default = reverse.
  -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
  -resume                         Re-start the pipeline if it has been previously run.
```

Explanation of parameters:

- `ref_genome`: genome sequence for mapping reads.
- `ref_ann`: annotation of genes/features in the reference genome.
- `sample_file`: TSV file containing sample information (see below)
- `data_dir`: path to directory containing FASTQ files.
- `paired`: data are paired-end (default is to assume single-end)
- `strandedness`: is data stranded? Options: `unstranded`, `forward`, `reverse`. Default = `reverse`.
- `cont_tabl`: (optional) table of contrasts to be performed for differential expression.
- `func_file`: (optional) functional annotation file - if provided, functional enrichment of DE genes will be performed.
- `p_thresh`: adjusted p-value threshold for identifying differentially expressed genes. Default = 0.05.
- `l2fc_thresh`: absolute log2(FoldChange) threshold for identifying differentially expressed genes. Default = 1.
- `skip_trimming`: do not trim adaptors from reads.
- `outdir`: the output directory where the results will be saved (Default: `./results`).

## Required inputs

- **Note**: See the [test data](https://github.com/adamd3/BactSeq/tree/main/test_data) folder for example mandatory inputs for a minimal run. Also see examples below.
- **Genome sequence**: FASTA file containing the genome sequence. Can be retrieved from NCBI.
- **Gene annotation file**: GFF file containing the genome annotation. Can be retrieved from NCBI.
- **Sample file**: TSV file containing sample information. Must contain the following columns:

  - `sample`: sample ID
  - `file_name`: name of the FASTQ file.
  - `group`: grouping factor for differential expression and exploratory plots.
  - `rep_no`: repeat number (if more than one sample per group).

  Example:

  If data are single-end, leave the `file2` column blank.

  ```console
  sample	file1   file2	group	rep_no
  AS_1	SRX1607051_T1.fastq.gz	    Artificial_Sputum	1
  AS_2	SRX1607052_T1.fastq.gz	    Artificial_Sputum	2
  AS_3	SRX1607053_T1.fastq.gz	    Artificial_Sputum	3
  MB_1	SRX1607054_T1.fastq.gz	    Middlebrook	1
  MB_2	SRX1607055_T1.fastq.gz	    Middlebrook	2
  MB_3	SRX1607056_T1.fastq.gz	    Middlebrook	3
  ER_1	SRX1607060_T1.fastq.gz	    Erythromycin	1
  ER_2	SRX1607061_T1.fastq.gz	    Erythromycin	2
  ER_3	SRX1607062_T1.fastq.gz	    Erythromycin	3
  KN_1	SRX1607066_T1.fastq.gz	    Kanamycin	1
  KN_2	SRX1607067_T1.fastq.gz	    Kanamycin	2
  KN_3	SRX1607068_T1.fastq.gz	    Kanamycin	3
  ```

## Optional inputs

- **Contrasts table**: TSV file containing contrasts to be performed to identify differentially expressed genes.
  Contains 2 columns, representing the groups (as defined in Samples file) to be contrasted.

  Example:

  ```console
  Condition1  Condition2
  Artificial_Sputum Middlebrook
  Artificial_Sputum Kanamycin
  Artificial_Sputum Erythromycin
  Middlebrook Erythromycin
  Middlebrook Kanamycin
  Erythromycin  Kanamycin
  ```

- **Functional annotation file**: CSV file containing functional categories for genes. Enrichment testing will be performed
  on results from differential gene expression contrasts. First column contains the gene ID (must match the gene IDs in `locus_tag` of the GFF annotation file); second column contains the functional groups (GO terms).

  Example:

  ```console
  MAB_0013c,"GO:0003674,GO:0003824,GO:0008150,GO:0008152,GO:0016407,GO:0016740,GO:0016746,GO:0016747"
  MAB_0018c,"GO:0003674,GO:0003824,GO:0005575,GO:0008150,GO:0008152,GO:0008168,GO:0016020,GO:0016021,GO:0016740,GO:0016741,GO:0031224,GO:0032259,GO:0044425"

  ```

## Output

1. **trim_galore** directory containing adaptor-trimmed RNA-Seq files and FastQC results.
2. **read_counts** directory containing:
   1. `ref_gene_df.tsv`: table of genes in the annotation.
   2. `gene_counts.tsv`: raw read counts per gene.
   3. `deseq_counts.tsv`: size factor-scaled, log2 counts matrix, normalised using DESeq2.
   4. `cpm_counts.tsv`: size factor-scaled, log2 counts per million (CPM) matrix, normalised using edgeR.
   5. `rpkm_counts.tsv`: size factor scaled and gene length-scaled counts, expressed as reads per kilobase per million mapped reads (RPKM).
3. **PCA_samples** directory containing principal component analysis results.
4. **diff_expr** directory containing differential expression results.
5. **func_enrich** directory containing functional enrichment results (optional).
