#!/usr/bin/env Rscript

## load / install packages
if (!require("optparse")){
    install.packages("optparse")
}
if (!require("edgeR")){
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("edgeR")
}

library(optparse)
library(edgeR)

option_list <- list(
    make_option(c("-t", "--log_transform"), type="character", default=NULL,
        help="log transform the counts? default = FALSE", metavar="character"),
    make_option(c("-o", "--outdir"), type="character", default=NULL,
        help="output directory for results", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

log <- if(opt$log_transform == "TRUE") TRUE else FALSE
outdir <- opt$outdir

## Read data
counts_tab <- read.csv(
    "gene_counts.tsv", header = TRUE, na.strings=c("","NA"), sep = "\t",
    stringsAsFactors = FALSE
)
ref_gene_tab <- read.csv(
    "ref_gene_df.tsv", header = TRUE, na.strings=c("","NA"), sep = "\t",
    stringsAsFactors = FALSE
)

## remove rRNA genes
ref_tab_sub <- ref_gene_tab[!ref_gene_tab$biotype=="rRNA",]
non_rRNA_counts <- counts_tab[rownames(counts_tab) %in% ref_tab_sub$locus_tag,]


## ensure that ref gene annotations order matches counts table
ref_tab_sub <- ref_tab_sub[match(rownames(non_rRNA_counts),ref_tab_sub$locus_tag),]

## counts per million (cpm) - normalised for lib size but not gene length
y <- DGEList(counts = non_rRNA_counts)
y <- calcNormFactors(y, method = "TMM")
cpm_df <- as.data.frame(edgeR::cpm(y, log = log))

write.table(
    cpm_df, file.path(outdir,"cpm_counts.tsv"), col.names = TRUE,
    row.names = TRUE, sep = "\t", quote = FALSE
)

## reads per kilobase per million (rpkm) - normalised for lib size + gene length
y <- DGEList(
    counts = non_rRNA_counts,
    genes = data.frame(gene.length = as.numeric(ref_tab_sub$gene_length))
)
y <- calcNormFactors(y)
rpkm_df <- as.data.frame(edgeR::rpkm(y, log = log))

write.table(
    rpkm_df, file.path(outdir,"rpkm_counts.tsv"), col.names = TRUE,
    row.names = TRUE, sep = "\t", quote = FALSE
)
