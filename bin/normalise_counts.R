#!/usr/bin/env Rscript

library(optparse)
library(edgeR)
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
counts_tab <- read.csv(
    "gene_counts_pc.tsv",
    header = TRUE, na.strings = c("", "NA"), sep = "\t",
    stringsAsFactors = FALSE
)
ref_gene_tab <- read.csv(
    "ref_gene_df.tsv",
    header = TRUE, na.strings = c("", "NA"), sep = "\t",
    stringsAsFactors = FALSE
)

gene_names <- counts_tab[["feature_id"]]
counts_tab[["feature_id"]] <- NULL
counts_tab <- as.data.frame(sapply(counts_tab, as.numeric))
rownames(counts_tab) <- gene_names

# ## remove rRNA genes
# ref_tab_sub <- ref_gene_tab[!ref_gene_tab$biotype == "rRNA", ]
# non_rRNA_counts <- counts_tab[rownames(counts_tab) %in% ref_tab_sub$locus_tag, ]
ref_tab_sub <- ref_gene_tab
non_rRNA_counts <- counts_tab

## ensure that ref gene annotations order matches counts table
ref_tab_sub <- ref_tab_sub[match(
    rownames(non_rRNA_counts), ref_tab_sub$locus_tag
), ]

## DESeq2 scaled counts
colData <- data.frame(sample_name = colnames(non_rRNA_counts))

dds <- DESeqDataSetFromMatrix(
    countData = round(non_rRNA_counts),
    colData = colData, design = ~1
)

dds <- estimateSizeFactors(dds)
deseq_norm <- counts(dds, normalized = TRUE)
deseq_df <- as_tibble(as.data.frame(log2(deseq_norm + 1)))
deseq_df$feature_id <- rownames(deseq_df)
deseq_df <- deseq_df[c("feature_id", colData$sample_name)]

write.table(
    deseq_df, file.path(outdir, "deseq_counts.tsv"),
    col.names = TRUE,
    row.names = FALSE, sep = "\t", quote = FALSE
)

## counts per million (cpm) - normalised for lib size but not gene length
y <- DGEList(counts = non_rRNA_counts)
y <- calcNormFactors(y, method = "TMM")
cpm_df <- as.data.frame(edgeR::cpm(y, log = log))
cpm_df <- tibble::rownames_to_column(as.data.frame(cpm_df), "feature_id")


write.table(
    cpm_df, file.path(outdir, "cpm_counts.tsv"),
    col.names = TRUE,
    row.names = FALSE, sep = "\t", quote = FALSE
)

## reads per kilobase per million (rpkm) - normalised for lib size + gene length
y <- DGEList(
    counts = non_rRNA_counts,
    genes = data.frame(gene.length = as.numeric(ref_tab_sub$gene_length))
)
y <- calcNormFactors(y)
rpkm_df <- as.data.frame(edgeR::rpkm(y, log = log))

rpkm_df <- tibble::rownames_to_column(as.data.frame(rpkm_df), "feature_id")

write.table(
    rpkm_df, file.path(outdir, "rpkm_counts.tsv"),
    col.names = TRUE,
    row.names = FALSE, sep = "\t", quote = FALSE
)
