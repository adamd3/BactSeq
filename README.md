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
6. Differential gene expression ([`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))
6. Functional enrichment of differentially expressed genes ([`fgsea`](https://bioconductor.org/packages/release/bioc/html/fgsea.html)) (optional)


## Installation

You will need to install [`Nextflow`](https://www.nextflow.io/) (version 21.10.3+).

You can run the pipeline as follows:

    nextflow run /home/adam/BactSeq \
        --data_dir [/path/to/data_dir] \
        --sample_file [sample_file.tsv] \
        --ref_genome [genome.fasta] --ref_ann [genome_annot.gff3] \
        -profile docker -resume

You can run with [`Docker`](https://www.docker.com/) or [`Singularity`](https://sylabs.io/guides/3.5/user-guide/introduction.html) by specifying ` -profile docker` or ` -profile singularity`, respectively.

The `-resume` parameter will re-start the pipeline if it has been previously run.

Explanation of parameters:
- `ref_genome`: genome sequence for mapping reads.
- `ref_ann`: annotation of genes/features in the reference genome.
- `sample_file`: TSV file containing sample information (see below)
- `data_dir`: path to directory containing FASTQ files.

Other available parameters:
- `func_file`: (optional) functional annotation file - if provided, functional enrichment of DE genes will be performed.
- `skip_trimming`: do not trim adaptors from reads.
- `outdir`: the output directory where the results will be saved (Default: `./results`).


## Required input

- __Genome sequence__: FASTA file containing the genome sequence. Can be retrieved from `https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/`.
- __Gene annotation file__: GFF file containing the genome annotation. Can be retrieved from `https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/`.
- __Sample file__: TSV file containing sample information. Must contain the following columns:
  - `sample`: sample ID
  - `file_name`: name of the FASTQ file.
  - `group`: grouping factor for differential expression and exploratory plots.
  - `rep_no`: repeat number (if more than one sample per group).

  Example:

    ```console
    sample	file_name	group	rep_no
    ZG301975	SRX5123742_T1.fastq.gz	ST_313	1
    ZG205864	SRX5123741_T1.fastq.gz	ST_111	1
    ZG302367	SRX5123744_T1.fastq.gz	ST_274	1
    ZG302359	SRX5123743_T1.fastq.gz	ST_244	1
    PSAE1647	SRX5123713_T1.fastq.gz	ST_313	1
    PSAE1649	SRX5123714_T1.fastq.gz	ST_313	1
    MHH17441	SRX5123695_T1.fastq.gz	ST_235	1
    PSAE1975	SRX5123726_T1.fastq.gz	ST_395	1
    PSAE1745	SRX5123719_T1.fastq.gz	ST_111	1
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
