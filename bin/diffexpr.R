#!/usr/bin/env Rscript

library(optparse)
library(DESeq2)
library(RColorBrewer)
library(RSQLite)
library(plyr)
library(tidyverse)

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


## -----------------------------------------------------------------------------
## Read and process data
## -----------------------------------------------------------------------------
counts_tab <- read_tsv(counts_f, na = c("", "NA"))
meta_tab <- read_tsv(meta_f)
contrast_tab <- read_tsv(cont_tab_f, col_types = cols(.default = "c"))

contrast_tab <- contrast_tab %>%
    mutate(across(everything(), ~ str_replace_all(str_remove_all(.x, "\\s+"), "-", ".")))

counts_tab <- counts_tab %>%
    mutate(across(-feature_id, as.numeric)) %>%
    column_to_rownames(var = "feature_id")

colnames(counts_tab) <- make.names(colnames(counts_tab))

meta_tab <- meta_tab %>%
    mutate(
        sample = make.names(str_replace_all(sample, "-", ".")),
        group  = factor(str_replace_all(group, "-", "."))
    )

common_samples <- intersect(colnames(counts_tab), meta_tab$sample)

if (length(common_samples) == 0) {
    stop("No matching samples found between counts_tab and meta_tab! Check your sample naming.")
}

counts_tab <- counts_tab[, common_samples, drop = FALSE]

meta_tab <- meta_tab %>%
    filter(sample %in% common_samples) %>%
    arrange(match(sample, common_samples)) %>%
    column_to_rownames(var = "sample")

stopifnot(identical(colnames(counts_tab), rownames(meta_tab)))


## -----------------------------------------------------------------------------
## Differential gene expression
## -----------------------------------------------------------------------------
comb_list <- contrast_tab %>%
    mutate(contrast_name = paste(.[[1]], .[[2]], sep = "_")) %>%
    {
        setNames(
            map2(.[[1]], .[[2]], c),
            .$contrast_name
        )
    }

dds <- DESeqDataSetFromMatrix(
    countData = round(counts_tab),
    colData   = meta_tab,
    design    = ~group
)
dds <- DESeq(dds)

contrast_list <- lapply(comb_list, function(x) {
    gp_1 <- x[1]
    gp_2 <- x[2]
    res <- lfcShrink(dds, contrast = c("group", gp_1, gp_2), type = "normal")
    res
})

lapply(seq_along(contrast_list), function(x) {
    contrast_name <- names(contrast_list)[x]

    res_df <- tibble::rownames_to_column(
        as.data.frame(contrast_list[[x]]), "feature_id"
    )

    write.table(
        res_df,
        file      = file.path(outdir, paste0("DGE_", contrast_name, ".tsv")),
        quote     = FALSE,
        sep       = "\t",
        row.names = FALSE,
        col.names = TRUE
    )
})


## -----------------------------------------------------------------------------
## Volcano plots
## -----------------------------------------------------------------------------
set1pal <- brewer.pal(9, "Set1")

lapply(seq_along(contrast_list), function(x) {
    res_obj <- contrast_list[[x]]

    # Get max -log10 padj safely (avoiding -Inf errors)
    valid_p <- na.omit(res_obj$padj)
    ymax <- ifelse(length(valid_p) > 0, max(-log10(valid_p[valid_p > 0])), 10)

    # Custom gene label selection
    keepLabs <- rownames(
        subset(res_obj, padj < p_thresh & abs(log2FoldChange) > l2fc_thresh)
    )

    # Custom colour scheme
    keyvals <- ifelse(
        res_obj$log2FoldChange < (-1 * l2fc_thresh) & res_obj$padj < p_thresh,
        set1pal[3],
        ifelse(
            res_obj$log2FoldChange > l2fc_thresh & res_obj$padj < p_thresh,
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
        res_obj,
        lab                = rownames(res_obj),
        selectLab          = keepLabs,
        x                  = "log2FoldChange",
        y                  = "padj",
        pointSize          = 3.0,
        labSize            = 4.0,
        pCutoff            = p_thresh,
        FCcutoff           = l2fc_thresh,
        colCustom          = keyvals,
        drawConnectors     = TRUE,
        max.overlaps       = 20,
        arrowheads         = FALSE,
        min.segment.length = 1.5,
        title              = paste0(group1, " vs ", group2),
        subtitle           = ""
    ) +
        theme(
            axis.text.x = element_text(size = 18),
            axis.title.x = element_text(size = 18),
            axis.text.y = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            plot.subtitle = element_blank(),
            plot.caption = element_blank(),
            legend.position = "none",
            plot.title = element_text(size = 18, hjust = 0.5)
        ) + ylab(bquote(~ -Log[10] ~ adjusted ~ italic(P)))

    ggsave(
        volcano_plot,
        file   = file.path(outdir, paste0("volcano_plot_", names(contrast_list)[x], ".png")),
        device = "png",
        units  = "in",
        width  = 8,
        height = 7,
        dpi    = 300
    )
})
