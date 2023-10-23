#!/usr/bin/env Rscript

library(optparse)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(RSQLite)
library(plyr)
library(tibble)

if (!require("EnhancedVolcano")) {
    if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    BiocManager::install("EnhancedVolcano")
    library(EnhancedVolcano)
}

option_list <- list(
    make_option(c("-p", "--p_threshold"),
        type = "double", default = 0.05,
        help = "adjusted p-value threshold (default = 0.05)",
        metavar = "character"
    ),
    make_option(c("-l", "--log2fc_threshold"),
        type = "double", default = 1,
        help = "absolute log2FoldChange threshold (default = 1)",
        metavar = "character"
    ),
    make_option(c("-o", "--outdir"),
        type = "character", default = NULL,
        help = "output directory for files", metavar = "character"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

counts_f <- "gene_counts_pc.tsv"
meta_f <- "sample_metadata.tsv"
cont_tab_f <- "contrast_table.tsv"

p_thresh <- opt$p_threshold
l2fc_thresh <- opt$log2fc_threshold
outdir <- opt$outdir




## ------------------------------------------------------------------------------
## Read and process data
## ------------------------------------------------------------------------------
counts_tab <- read.csv(
    counts_f,
    header = TRUE, na.strings = c("", "NA"), sep = "\t",
    stringsAsFactors = FALSE
)
meta_tab <- read.table(
    meta_f,
    header = TRUE, sep = "\t", stringsAsFactors = FALSE
)
contrast_tab <- read.table(
    cont_tab_f,
    header = TRUE, sep = "\t", stringsAsFactors = FALSE
)
contrast_tab <- data.frame(
    rbind(apply(contrast_tab, 2, function(x) gsub("\\s+", "", x)))
)


gene_names <- counts_tab[["feature_id"]]
counts_tab[["feature_id"]] <- NULL
counts_tab <- as.data.frame(sapply(counts_tab, as.numeric))
rownames(counts_tab) <- gene_names


## factorise group column
meta_tab$group <- as.factor(as.character(meta_tab$group))

## order rows to match counts columns
meta_tab <- meta_tab[match(colnames(counts_tab), meta_tab$sample), ]


## ------------------------------------------------------------------------------
## Differential gene expression
## ------------------------------------------------------------------------------

## make list of contrasts to be performed from contrast table
comb_list <- lapply(1:nrow(contrast_tab), function(idx) {
    c(contrast_tab[idx, 1], contrast_tab[idx, 2])
})
comb_names <- lapply(1:nrow(contrast_tab), function(idx) {
    contrast_name <- paste0(contrast_tab[idx, 1], "_", contrast_tab[idx, 2])
    contrast_name
})
names(comb_list) <- comb_names

# ## previous version - get all possible combinations of group column:
# comb_list <- combn(levels(meta_tab$group), 2, simplify = FALSE)
#
# comb_names <- sapply(comb_list, function(x){
#     contrast_name <- paste0(x[1],"_",x[2])
#     contrast_name
# })
# names(comb_list) <- comb_names

## make DESeq2 object
dds <- DESeqDataSetFromMatrix(
    countData = round(counts_tab), colData = meta_tab,
    design = ~group
)
dds <- DESeq(dds)

## perform pairwise contrasts of groups
contrast_list <- lapply(comb_list, function(x) {
    gp_1 <- x[1]
    gp_2 <- x[2]
    res <- lfcShrink(dds, contrast = c("group", gp_1, gp_2), type = "normal")
    res
})

## export tables of genes with log2FC and p-values:
lapply(seq_along(contrast_list), function(x) {
    contrast_name <- names(contrast_list)[x]
    res_df <- tibble::rownames_to_column(as.data.frame(contrast_list[x]), "feature_id")
    colnames(res_df) <- gsub(contrast_name, "", colnames(res_df))
    colnames(res_df) <- gsub("\\.", "", colnames(res_df))
    write.table(
        res_df,
        file = file.path(outdir, paste0("DGE_", contrast_name, ".tsv")),
        quote = FALSE, sep = "\t",
        row.names = FALSE, col.names = TRUE
    )
})



## ------------------------------------------------------------------------------
## Volcano plots
## ------------------------------------------------------------------------------
## choose the limits for the x- and y-axes using the log2FoldChanges and pvalues
lfc_list <- lapply(contrast_list, function(x) {
    subset(x, padj < p_thresh)$log2FoldChange
})
xmin <- round_any(min(unlist(lfc_list)), 0.1, floor)
xmax <- round_any(max(unlist(lfc_list)), 0.1, ceiling)

if (abs(xmax) > abs(xmin)) {
    xmin <- xmax * -1
} else if (abs(xmin) > abs(xmax)) {
    xmax <- xmin * -1
}

ymin <- 0
pval_list <- lapply(contrast_list, function(x) {
    subset(x, padj < p_thresh)$padj
})
ymax <- round_any(max(unlist(pval_list)), 0.5, ceiling)

set1pal <- brewer.pal(9, "Set1")


lapply(seq_along(contrast_list), function(x) {
    ymin <- 0
    ymax <- max(-log10(
        subset(contrast_list[[x]], padj < p_thresh & log2FoldChange > l2fc_thresh)$padj
    ))
    xmin <- round_any(
        min(unlist(subset(
            contrast_list[[x]], padj < p_thresh & log2FoldChange > l2fc_thresh
        )$log2FoldChange)),
        0.1, floor
    )
    xmax <- round_any(
        max(unlist(subset(
            contrast_list[[x]], padj < p_thresh & log2FoldChange > l2fc_thresh
        )$log2FoldChange)),
        0.1, ceiling
    )
    if (abs(xmax) > abs(xmin)) {
        xmin <- xmax * -1
    } else if (abs(xmin) > abs(xmax)) {
        xmax <- xmin * -1
    }

    ## show gene labels where the log2FoldChange exceeds a certain threshold:
    keepLabs <- rownames(
        subset(contrast_list[[x]], padj < p_thresh & (abs(log2FoldChange) > xmax * 0.9))
    )
    keepLabs <- c(keepLabs, rownames(
        subset(contrast_list[[x]], padj < p_thresh & (abs(-log10(padj)) > ymax * 0.5))
    ))

    ## Custom colour scheme
    keyvals <- ifelse(
        contrast_list[[x]]$log2FoldChange < (l2fc_thresh * -1) & contrast_list[[x]]$padj < p_thresh,
        set1pal[3],
        ifelse(
            contrast_list[[x]]$log2FoldChange > l2fc_thresh & contrast_list[[x]]$padj < p_thresh,
            set1pal[1], "grey45"
        )
    )
    keyvals[is.na(keyvals)] <- "grey45"
    names(keyvals)[keyvals == set1pal[1]] <- "up"
    names(keyvals)[keyvals == "grey45"] <- "NS"
    names(keyvals)[keyvals == set1pal[3]] <- "down"

    group1 <- (comb_list[[x]])[1]
    group2 <- (comb_list[[x]])[2]

    volcano_plot <- EnhancedVolcano::EnhancedVolcano(
        contrast_list[[x]],
        lab = rownames(contrast_list[[x]]),
        selectLab = keepLabs,
        x = "log2FoldChange", y = "padj",
        pointSize = 3.0,
        labSize = 4.0,
        pCutoff = p_thresh,
        FCcutoff = l2fc_thresh,
        colCustom = keyvals,
        drawConnectors = TRUE,
        max.overlaps = 20,
        arrowheads = FALSE,
        min.segment.length = 1.5,
        title = paste0(group1, " vs ", group2),
        subtitle = ""
    ) +
        theme(
            axis.text.x = element_text(size = 12 * 1.5),
            axis.title.x = element_text(size = 12 * 1.5),
            axis.text.y = element_text(size = 12 * 1.5),
            axis.title.y = element_text(size = 12 * 1.5),
            plot.subtitle = element_blank(),
            plot.caption = element_blank(),
            legend.position = "none",
            plot.title = element_text(size = 12 * 1.5, hjust = 0.5)
        ) + ylab(bquote(~ -Log[10] ~ adjusted ~ italic(P)))

    ggsave(
        volcano_plot,
        file = file.path(outdir, paste0("volcano_plot_", names(contrast_list)[x], ".png")),
        device = "png", units = "in",
        width = 8, height = 7, dpi = 300
    )
})
