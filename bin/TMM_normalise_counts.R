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

option_list <- list(
    make_option(c("-t", "--log_transform"), type="character", default=NULL,
        help="log transform the counts? default = FALSE", metavar="character"),
    make_option(c("-o", "--outf"), type="character", default=NULL,
        help="output file for results", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

log <- if(opt$log_transform == "TRUE") TRUE else FALSE
outf <- opt$outf

## Read data
counts_tab <- read.csv(
    "gene_counts.tsv", header = TRUE, na.strings=c("","NA"), sep = "\t",
    stringsAsFactors = FALSE
)

y <- DGEList(counts = counts_tab)
y <- calcNormFactors(y, method = "TMM")
libSizes <- y$samples$lib.size
res_df <- as.data.frame(cpm(
    y, log = log, lib.size = (libSizes)*(y$samples$norm.factors)
))


write.table(
    res_df, outf, col.names = TRUE, row.names = TRUE,
    sep = "\t", quote = FALSE
)
