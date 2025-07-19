# BactSeq

![BactSeq](https://github.com/adamd3/BactSeq/actions/workflows/ci.yml/badge.svg)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**BactSeq** is a Nextflow pipeline for performing bacterial RNA-Seq analysis. The pipeline supports both BWA alignment and Kallisto pseudo-alignment approaches, with optional differential expression and functional enrichment analysis.

## Quick Start

```bash
# Run with test data
nextflow run BactSeq -profile test,docker

# Run with your own data
nextflow run BactSeq \
  --data_dir /path/to/fastq/files \
  --sample_file samples.tsv \
  --ref_genome genome.fasta \
  --ref_ann genome.gff3 \
  -profile docker
```

## Pipeline Summary

The pipeline performs the following steps:

### Core Analysis

1. **Quality Control & Trimming**

   - Trim adaptors from reads ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
   - Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))

2. **Read Alignment** (choose one)

   - **BWA**: Align reads to reference genome ([`BWA-MEM`](https://github.com/lh3/bwa/))
   - **Kallisto**: Pseudo-align reads to coding sequences ([`Kallisto`](https://pachterlab.github.io/kallisto/))

3. **Expression Quantification & Normalization**

   - Size-factor scaling using [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
   - TMM normalization using [`edgeR`](http://bioconductor.org/packages/release/bioc/html/edgeR.html)
   - RPKM scaling for gene length normalization

4. **Exploratory Analysis**
   - Principal component analysis (PCA) of normalized expression values
   - Sample clustering and visualization

### Optional Analysis

5. **Differential Expression** ([`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))

   - Pairwise comparisons based on provided contrasts
   - Volcano plots and summary statistics

6. **Functional Enrichment** ([`topGO`](https://bioconductor.org/packages/release/bioc/html/topGO.html))
   - GO term enrichment analysis of differentially expressed genes
   - Enrichment plots and gene lists

## Installation

You will need to install [`Nextflow`](https://www.nextflow.io/) (version 21.10.3+).

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash

# Or via conda
conda install -c bioconda nextflow
```

## Usage

### Basic Usage

```bash
nextflow run BactSeq \
  --data_dir /path/to/fastq/files \
  --sample_file samples.tsv \
  --ref_genome genome.fasta \
  --ref_ann genome.gff3 \
  -profile docker
```

### Parameters

#### Required Parameters

| Parameter       | Description                                                                                       |
| --------------- | ------------------------------------------------------------------------------------------------- |
| `--data_dir`    | Path to directory containing FastQ files                                                          |
| `--sample_file` | Path to file containing sample information                                                        |
| `--ref_genome`  | Path to FASTA file containing reference genome sequence (BWA) or coding gene sequences (Kallisto) |
| `-profile`      | Configuration profile to use: `conda`, `docker`, `singularity`                                    |

#### Optional Parameters

| Parameter         | Default     | Description                                                                |
| ----------------- | ----------- | -------------------------------------------------------------------------- |
| `--aligner`       | `bwa`       | Aligner to use: `bwa`, `kallisto`                                          |
| `--ref_ann`       | -           | Path to GFF file containing reference genome annotation (required for BWA) |
| `--contrast_file` | -           | Path to TSV file containing contrasts for differential expression          |
| `--func_file`     | -           | Path to functional annotation file for enrichment analysis                 |
| `--strandedness`  | `reverse`   | Data strandedness: `unstranded`, `forward`, `reverse`                      |
| `--p_thresh`      | `0.05`      | Adjusted p-value threshold for differential expression                     |
| `--l2fc_thresh`   | `1`         | Absolute log2(FoldChange) threshold for differential expression            |
| `--fragment_len`  | `150`       | Average fragment length for Kallisto (single-end only)                     |
| `--fragment_sd`   | `20`        | Fragment length standard deviation for Kallisto (single-end only)          |
| `--skip_trimming` | `false`     | Skip adapter trimming                                                      |
| `--outdir`        | `./results` | Output directory for results                                               |

#### Standard Nextflow Parameters

| Parameter   | Description                        |
| ----------- | ---------------------------------- |
| `-resume`   | Resume a previous run              |
| `-name`     | Name for the pipeline run          |
| `-work-dir` | Work directory for temporary files |

## Input Files

> **Note**: See the [test data](https://github.com/adamd3/BactSeq/tree/main/test_data) folder for example inputs.

### Required Files

1. **Sample Sheet** (`samples.tsv`)

   - TSV file containing sample information with the following columns:
     - `sample`: Sample ID
     - `file1`: Name of R1 FastQ file
     - `file2`: Name of R2 FastQ file (leave blank for single-end)
     - `group`: Grouping factor for differential expression
     - `rep_no`: Replicate number

   **Example:**

   ```tsv
   sample	file1	file2	group	rep_no
   AS_1	SRX1607051_T1.fastq.gz		Artificial_Sputum	1
   AS_2	SRX1607052_T1.fastq.gz		Artificial_Sputum	2
   MB_1	SRX1607054_T1.fastq.gz		Middlebrook	1
   MB_2	SRX1607055_T1.fastq.gz		Middlebrook	2
   ```

2. **Reference Genome** (`genome.fasta`)

   - FASTA file containing the reference genome sequence
   - Can be downloaded from NCBI RefSeq

3. **Gene Annotation** (`genome.gff3`) - _Required for BWA only_
   - GFF3 file containing gene annotations
   - Can be downloaded from NCBI RefSeq

### Optional Files

1. **Contrasts Table** (`contrasts.tsv`) - _For differential expression_

   - TSV file with 2 columns defining comparisons to perform
   - Column names: `Condition1`, `Condition2`

   **Example:**

   ```tsv
   Condition1	Condition2
   Artificial_Sputum	Middlebrook
   Artificial_Sputum	Kanamycin
   Middlebrook	Kanamycin
   ```

2. **Functional Annotation File** (`functional_annotation.csv`) - _For enrichment analysis_

   - CSV file containing GO terms for genes
   - Column 1: Gene ID (must match `locus_tag` in GFF)
   - Column 2: GO terms (comma-separated)

   **Example:**

   ```csv
   Gene,GO_terms
   MAB_0001,"GO:0005737,GO:0008150"
   MAB_0002,"GO:0016020,GO:0006810"
   ```

## Output

The pipeline generates the following output directories:

| Directory        | Contents                                           |
| ---------------- | -------------------------------------------------- |
| `trim_galore/`   | Adapter-trimmed FastQ files and FastQC reports     |
| `read_counts/`   | Raw and normalized gene count matrices             |
| `PCA_samples/`   | Principal component analysis plots and coordinates |
| `diff_expr/`     | Differential expression results and volcano plots  |
| `func_enrich/`   | Functional enrichment analysis results             |
| `pipeline_info/` | Pipeline execution reports and metadata            |

### Key Output Files

**Count Matrices:**

- `gene_counts.tsv`: Raw read counts per gene
- `deseq_counts.tsv`: DESeq2 normalized counts (log2 transformed)
- `cpm_counts.tsv`: Counts per million (CPM) normalized
- `rpkm_counts.tsv`: RPKM normalized counts

**Analysis Results:**

- `DGE_*.tsv`: Differential expression results for each contrast
- `volcano_plot_*.png`: Volcano plots for each contrast
- `pca_grouped.png`: PCA plot colored by sample groups
- `*_enrich.tsv`: GO enrichment results (if functional annotation provided)

## Citation

If you use BactSeq in your research, please cite:

```bibtex
@software{bactseq,
  author = {Adam Dinan},
  title = {BactSeq: A Nextflow pipeline for bacterial RNA-seq analysis},
  url = {https://github.com/adamd3/BactSeq},
  version = {dev},
  year = {2024}
}
```

## Changelog

### Recent Updates

- **2025-07-19**: Docker image is now compatible with nextflow 25.04.0+
- **2023-10-06**: Docker and Singularity support restored
- **2023-09-28**: Added example contrasts table and functional enrichment file
