#!/usr/bin/env Rscript

## load / install packages
if (!require("optparse")){
    install.packages("optparse")
}
if (!require("ggplot2")){
    install.packages("ggplot2")
}
if (!require("RColorBrewer")){
    install.packages("RColorBrewer")
}
if (!require("DESeq2")){
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("DESeq2")
}
if (!require("EnhancedVolcano")){
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("EnhancedVolcano")
}


option_list <- list(
    make_option(c("-o", "--outdir"), type="character", default=NULL,
        help="output directory for files", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

counts_f <- "cpm_counts.tsv"
meta_f <- "sample_metadata.tsv"
outdir <- opt$outdir


large_disc_pal <- brewer.pal.info[brewer.pal.info$category == 'qual',]
colpal_large <- unlist(
    mapply(brewer.pal, large_disc_pal$maxcolors, rownames(large_disc_pal)))
colpal_large[c(5:8)] <- colpal_large[c(70:73)] ## replace to avoid colour clashes


##------------------------------------------------------------------------------
## Read and process data
##------------------------------------------------------------------------------
counts_tab <- read.csv(
    "gene_counts.tsv", header = TRUE, na.strings=c("","NA"), sep = "\t",
    stringsAsFactors = FALSE
)
meta_tab <- read.table(
    meta_f, header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

## factorise group column
meta_tab$group <- as.factor(as.character(meta_tab$group))

## order rows to match counts columns
meta_tab <- meta_tab[match(colnames(counts_tab),meta_tab$sample),]



##------------------------------------------------------------------------------
## Differential gene expression
##------------------------------------------------------------------------------
## get all possible combinations of group column:
comb_list <- combn(levels(meta_tab$group), 2, simplify = FALSE)

comb_names <- sapply(comb_list, function(x){
    contrast_name <- paste0(x[1],"_",x[2])
    contrast_name
})
names(comb_list) <- comb_names

## make DESeq2 object
dds <- DESeqDataSetFromMatrix(
    countData = counts_tab, colData = meta_tab,
    design = ~ group
)
dds <- DESeq(dds)

## perform pairwise contrasts of groups
contrast_list <- lapply(comb_list, function(x){
    gp_1 <- (comb_list[[x]])[1]
    gp_2 <- (comb_list[[x]])[2]
    res <- lfcShrink(dds, contrast = c("group", gp_1, gp_2), type = "normal")
    res
})


## export tables of genes with log2FC and p-values:
lapply(seq_along(contrast_list), function(x){
    contrast_name <- names(contrast_list)[x]
    write.table(
        contrast_list[x],
        file = file.path(outdir, paste0("DGE_", contrast_name, ".tsv")),
        quote = FALSE, sep = "\t",
        row.names = TRUE, col.names = TRUE
    )
})
