#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(readr)
library(tibble)

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
# norm_counts <- read.table(
#     counts_f,
#     header = TRUE, sep = "\t", stringsAsFactors = FALSE
# )
# meta_tab <- read.table(
#     meta_f,
#     header = TRUE, sep = "\t", stringsAsFactors = FALSE
# )

norm_counts <- read_tsv(counts_f)
meta_tab <- read_tsv(meta_f)


gene_names <- norm_counts[["feature_id"]]
norm_counts[["feature_id"]] <- NULL
norm_counts <- as.data.frame(norm_counts)
rownames(norm_counts) <- gene_names

## factorise group
meta_tab$group <- as.factor(as.character(meta_tab$group))

## order rows to match counts columns
meta_tab <- meta_tab[match(colnames(norm_counts), meta_tab$sample), ]



## ------------------------------------------------------------------------------
## PCA
## ------------------------------------------------------------------------------
pca_counts <- prcomp(t(norm_counts), center = TRUE, scale = FALSE)

pca_coords <- data.frame(pca_counts$x)
pca_coords$sample <- rownames(pca_coords)

# move sample to first column
pca_coords <- pca_coords[c(
    "sample", colnames(pca_coords)[1:ncol(pca_coords) - 1]
)]
rownames(pca_coords) <- NULL

## save pca object
saveRDS(pca_counts, file.path(outdir, "pca.rds"))

## save pca coordinates
write.table(
    pca_coords,
    file = file.path(outdir, "pca_coords.tsv"),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
)

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
