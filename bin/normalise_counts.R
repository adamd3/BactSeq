#!/usr/bin/env Rscript

library(optparse)
library(edgeR)
library(readr)
library(tibble)
library(DESeq2)

option_list <- list(
    make_option(c("-t", "--log_transform"),
        type = "character", default = NULL,
        help = "log transform the counts? default = FALSE",
        metavar = "character"
    ),
    make_option(c("-o", "--outdir"),
        type = "character", default = NULL,
        help = "output directory for results", metavar = "character"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

log <- if (opt$log_transform == "TRUE") TRUE else FALSE
outdir <- opt$outdir

## Read data
counts_tab <- read_tsv("gene_counts_pc.tsv")
ref_gene_tab <- read_tsv("ref_gene_df.tsv")

gene_names <- counts_tab[["feature_id"]]
counts_tab <- counts_tab %>%
    select(-feature_id) %>%
    as.data.frame()
rownames(counts_tab) <- gene_names

ref_tab_sub <- ref_gene_tab
non_rRNA_counts <- counts_tab

## ensure that ref gene annotations order matches counts table
ref_tab_sub <- ref_tab_sub[match(
    rownames(non_rRNA_counts), ref_tab_sub$locus_tag
), ]

## DESeq2 scaled counts
colData <- tibble(sample_name = colnames(non_rRNA_counts))

dds <- DESeqDataSetFromMatrix(
    countData = round(non_rRNA_counts),
    colData = colData, design = ~1
)

dds <- estimateSizeFactors(dds)
deseq_norm <- counts(dds, normalized = TRUE)
deseq_df <- as_tibble(log2(deseq_norm + 1)) %>%
    rownames_to_column("feature_id") %>%
    select(feature_id, all_of(colData$sample_name))

write_tsv(deseq_df, file.path(outdir, "deseq_counts.tsv"))

## counts per million (cpm) - normalised for lib size but not gene length
y <- DGEList(counts = non_rRNA_counts)
y <- calcNormFactors(y, method = "TMM")
cpm_df <- edgeR::cpm(y, log = log) %>%
    as_tibble(rownames = "feature_id")

write_tsv(cpm_df, file.path(outdir, "cpm_counts.tsv"))

## reads per kilobase per million (rpkm) - normalised for lib size + gene length
y <- DGEList(
    counts = non_rRNA_counts,
    genes = tibble(gene.length = as.numeric(ref_tab_sub$gene_length))
)
y <- calcNormFactors(y)
rpkm_df <- edgeR::rpkm(y, log = log) %>%
    as_tibble(rownames = "feature_id")

write_tsv(rpkm_df, file.path(outdir, "rpkm_counts.tsv"))
