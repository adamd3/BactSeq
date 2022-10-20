# BactSeq

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**BactSeq** is a Nextflow pipeline for performing bacterial RNA-Seq analysis.

## Pipeline summary

The pipeline will perform the following steps:

1. Trim adaptors from reads ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Align reads to reference genome ([`BWA-MEM`](https://github.com/lh3/bwa/))
4. Size-factor scaling and gene length (RPKM) scaling of counts (TMM from [`edgeR`](http://bioconductor.org/packages/release/bioc/html/edgeR.html))
5. Principal component analysis (PCA) of normalised expression values
6. Differential gene expression ([`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)) (optional)
6. Functional enrichment of differentially expressed genes ([`topGO`](https://bioconductor.org/packages/release/bioc/html/topGO.html)) (optional)


## Installation

You will need to install [`Nextflow`](https://www.nextflow.io/) (version 21.10.3+).

You can run the pipeline as follows:

    nextflow run /home/adam/BactSeq \
        --data_dir [/path/to/data_dir] \
        --sample_file [sample_file.tsv] \
        --ref_genome [genome.fasta] --ref_ann [genome_annot.gff3] \
        --func_file [func_file.csv] \
        -profile docker -resume

You can run with [`Docker`](https://www.docker.com/) or [`Singularity`](https://sylabs.io/guides/3.5/user-guide/introduction.html) by specifying ` -profile docker` or ` -profile singularity`, respectively.

The `-resume` parameter will re-start the pipeline if it has been previously run.

Explanation of parameters:
- `ref_genome`: genome sequence for mapping reads.
- `ref_ann`: annotation of genes/features in the reference genome.
- `sample_file`: TSV file containing sample information (see below)
- `data_dir`: path to directory containing FASTQ files.
- `paired`: data are paired-end (default is to assume single-end)
- `strandedness`: is data stranded? Options: `unstranded`, `forward`, `reverse`. Default = reverse.
- `cont_tabl`: (optional) table of contrasts to be performed for differential expression.
- `func_file`: (optional) functional annotation file - if provided, functional enrichment of DE genes will be performed.
- `p_thresh`: adjusted p-value threshold for identifying differentially expressed genes. Default = 0.05.
- `l2fc_thresh`: absolute log2(FoldChange) threshold for identifying differentially expressed genes. Default = 1.
- `skip_trimming`: do not trim adaptors from reads.
- `outdir`: the output directory where the results will be saved (Default: `./results`).


## Required input

- __Genome sequence__: FASTA file containing the genome sequence. Can be retrieved from NCBI.
- __Gene annotation file__: GFF file containing the genome annotation. Can be retrieved from NCBI.
- __Sample file__: TSV file containing sample information. Must contain the following columns:
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

## Output

1. __trim_galore__ directory containing adaptor-trimmed RNA-Seq files and FastQC results.
2. __read_counts__ directory containing:
    1. `ref_gene_df.tsv`: table of genes in the annotation.
    2. `gene_counts.tsv`: raw read counts per gene.
    3. `cpm_counts.tsv`: size factor scaled counts per million (CPM).
    4. `rpkm_counts.tsv`: size factor scaled and gene length-scaled counts, expressed as reads per kilobase per million mapped reads (RPKM).
3. __PCA_samples__ directory containing principal component analysis results.
4. __diff_expr__ directory containing differential expression results.
5. __func_enrich__ directory containing functional enrichment results (optional).
