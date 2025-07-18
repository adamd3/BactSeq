#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(readr)
library(tibble)
library(dplyr)

option_list <- list(
    make_option(c("-o", "--outdir"),
        type = "character", default = NULL,
        help = "output directory for files", metavar = "character"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

counts_f <- "cpm_counts.tsv"
meta_f <- "sample_metadata.tsv"
outdir <- opt$outdir


large_disc_pal <- brewer.pal.info[brewer.pal.info$category == "qual", ]
colpal_large <- unlist(
    mapply(brewer.pal, large_disc_pal$maxcolors, rownames(large_disc_pal))
)
colpal_large[c(5:8)] <- colpal_large[c(70:73)] ## replace to avoid colour clashes
colpal_large <- c(brewer.pal(9, "Set1"), colpal_large)


## ------------------------------------------------------------------------------
## Read and process data
## ------------------------------------------------------------------------------
norm_counts <- read_tsv(counts_f)
meta_tab <- read_tsv(meta_f)


gene_names <- norm_counts[["feature_id"]]
norm_counts <- norm_counts %>%
    select(-feature_id) %>%
    as.data.frame()
rownames(norm_counts) <- gene_names

meta_tab <- meta_tab %>%
    mutate(group = as.factor(as.character(group))) %>%
    arrange(match(sample, colnames(norm_counts)))



## ------------------------------------------------------------------------------
## PCA
## ------------------------------------------------------------------------------
pca_counts <- prcomp(t(norm_counts), center = TRUE, scale = FALSE)

pca_coords <- pca_counts$x %>%
    as_tibble(rownames = "sample") %>%
    select(sample, everything())

## save pca object
saveRDS(pca_counts, file.path(outdir, "pca.rds"))

## save pca coordinates
write_tsv(pca_coords, file.path(outdir, "pca_coords.tsv"))

p1 <- ggplot(pca_coords, aes(x = PC1, y = PC2)) +
    geom_point(
        size = 4, shape = 21, colour = "black",
        aes(fill = meta_tab$group)
    ) +
    theme_bw(base_size = 12) +
    scale_colour_manual(values = colpal_large, guide = "none") +
    scale_fill_manual("Group", values = colpal_large) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    theme(
        legend.position = "right",
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 12)
    )

ggsave(
    p1,
    file = file.path(outdir, "pca_grouped.png"),
    device = "png", units = "in",
    width = 9, height = 7, dpi = 300
)
