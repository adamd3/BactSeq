#!/usr/bin/env Rscript

## load / install packages
if (!require("optparse")){
    install.packages("optparse",repos = "http://cran.us.r-project.org")
}
if (!require("readr")){
    install.packages("readr",repos = "http://cran.us.r-project.org")
}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos = "http://cran.us.r-project.org")
if (!require("fgsea")){
    BiocManager::install("fgsea")
}
if (!require("ggplot2")){
    install.packages("ggplot2",repos = "http://cran.us.r-project.org")
}

library(optparse)
library(readr)
library(fgsea)
library(ggplot2)

option_list <- list(
    make_option(c("-r", "--res_f"), type="character", default=NULL,
        help="tab-delimited results table from differential expression analysis",
        metavar="character"),
    make_option(c("-a", "--annot_f"), type="character", default=NULL,
        help="GMT-format file containing functional annotations",
        metavar="character"),
    make_option(c("-p", "--p_threshold"), type="double", default=0.05,
        help="adjusted p-value threshold (default = 0.05)",
        metavar="character"),
    make_option(c("-t", "--log2fc_threshold"), type="double", default=1,
        help="absolute log2FoldChange threshold (default = 1)",
        metavar="character"),
    make_option(c("-x", "--max_size"), type="double", default=200,
        help="maximum size of the functional categories to test (default = 200)",
        metavar="character"),
    make_option(c("-n", "--min_size"), type="double", default=10,
        help="minimum size of the functional categories to test (default = 10)",
        metavar="character"),
    make_option(c("-l", "--label"), type="character", default=10,
        help="prefix name for output files",
        metavar="character"),
    make_option(c("-o", "--out_dir"), type="character", default=NULL,
        help="output dir for results", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

res_f <- opt$res_f
annot_f <- opt$annot_f
p_threshold <- opt$p_threshold
log2fc_threshold <- opt$log2fc_threshold
max_size <- opt$max_size
min_size <- opt$min_size
label <- opt$label
out_dir <- opt$out_dir


# res_f <- "/home/adam/BactSeq/GO_annotations/Full_DEG_CRISPR_MAB_4099c.tsv"
# annot_f <- "/home/adam/BactSeq/GO_annotations/Mabs_GO_BP.gmt"
# p_threshold <- 0.05
# log2fc_threshold <- 1
# max_size <- 200
# min_size <- 10
# label <- "CRISPR_MAB_4099c"
# out_dir <- "/home/adam/BactSeq/GO_annotations"


# Rscript functional_enrichment.R \
#     -r /home/adam/BactSeq/GO_annotations/Full_DEG_CRISPR_MAB_4099c.tsv \
#     -a /home/adam/BactSeq/GO_annotations/Mabs_GO_BP.gmt \
#     -l "CRISPR_MAB_4099c" \
#     -o /home/adam/BactSeq/GO_annotations



##------------------------------------------------------------------------------
## Read data
##------------------------------------------------------------------------------
res_tab <- read.table(res_f)
annot_tab <- read_tsv(annot_f, col_names=FALSE, show_col_types=FALSE)[c(1,2)]
annot_list <- gmtPathways(file.path(annot_f))

if(!length(annot_list)==nrow(annot_tab)){
    stop("annotation data was not successfully parsed. please check files.")
}



##------------------------------------------------------------------------------
## Perform enrichment testing
##------------------------------------------------------------------------------

res_tab <- na.omit(res_tab)
gene_universe <- rownames(res_tab)

padj_col <- colnames(res_tab)[grepl("padj", colnames(res_tab))]
l2fc_col <- colnames(res_tab)[grepl("log2FoldChange", colnames(res_tab))]

genes_up <- rownames(res_tab[(res_tab[[padj_col]] < p_threshold &
    res_tab[[l2fc_col]] > log2fc_threshold),])
genes_down <- rownames(res_tab[(res_tab[[padj_col]] < p_threshold &
    res_tab[[l2fc_col]] < (log2fc_threshold*-1)),])


fora_up <- fora(annot_list, genes_up, gene_universe, minSize = min_size, maxSize=max_size)
fora_down <- fora(annot_list, genes_down, gene_universe, minSize = min_size, maxSize=max_size)


###################################
## annotate + export results tables
###################################
colnames(annot_tab) <- c("pathway", "description")

fora_up <- merge(annot_tab, fora_up, all.y=TRUE)
fora_down <- merge(annot_tab, fora_down, all.y=TRUE)

fora_up <- fora_up[order(fora_up$pval),]
fora_down <- fora_down[order(fora_down$pval),]


fora_up_export <- as.data.frame(apply(fora_up,2,as.character))
fora_down_export <- as.data.frame(apply(fora_down,2,as.character))

write.table(
    fora_up_export, file.path(out_dir,paste0(label,"_up_enrich.tsv")),
    col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE
)
write.table(
    fora_down_export, file.path(out_dir,paste0(label,"_down_enrich.tsv")),
    col.names = TRUE,row.names = FALSE, sep = "\t", quote = FALSE
)


#############################
## export barplots of pvalues
#############################
up_col <- "#E41A1C"
down_col <- "#7FC97F"

## select top 20 terms for plotting
fora_up_plot <- fora_up_export[1:20,]
fora_down_plot <- fora_down_export[1:20,]

fora_up_plot[["order"]] <- factor(1:nrow(fora_up_plot), levels=rev(1:nrow(fora_up_plot)))
fora_down_plot[["order"]] <- factor(1:nrow(fora_down_plot), levels=rev(1:nrow(fora_down_plot)))

fora_up_plot[["pval"]] <- as.numeric(fora_up_plot[["pval"]])
fora_down_plot[["pval"]] <- as.numeric(fora_down_plot[["pval"]])

p1 <- ggplot(fora_up_plot,
    aes(x = order, y = -log10(pval))) +
    geom_bar(
        position = "dodge", stat = "identity",
        width = 0.7,  colour = "black", fill = up_col
    ) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    theme_bw() +
    ylab(expression("−log"[10]*"("*italic(P)~"value"*")")) +
    scale_x_discrete(labels = rev(fora_up_plot$description)) +
    coord_flip() +
    theme(
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(colour = "black", size = 12*1.25),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black", size = 12*1.25),
        axis.text.y = element_text(colour = "black", size = 12*1.25),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.text = element_text(colour = "black", size = 12*1.25),
        legend.title = element_text(colour = "black", size = 12*1.25)
    )

ggsave(
    p1, file = file.path(out_dir, paste0(label,"_up_enrich.png")),
    device = "png",
    width = 8, height = 8,
    dpi = 400
)


p1 <- ggplot(fora_down_plot,
    aes(x = order, y = -log10(pval))) +
    geom_bar(
        position = "dodge", stat = "identity",
        width = 0.7,  colour = "black",  fill = down_col
    ) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    theme_bw() +
    ylab(expression("−log"[10]*"("*italic(P)~"value"*")")) +
    scale_x_discrete(labels = rev(fora_down_plot$description)) +
    coord_flip() +
    theme(
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 12*1.25),
        axis.title.x = element_text(colour = "black", size = 12*1.25),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black", size = 12*1.25),
        axis.text.y = element_text(colour = "black", size = 12*1.25),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.text = element_text(colour = "black", size = 12*1.25),
        legend.title = element_text(colour = "black", size = 12*1.25)
    )

ggsave(
    p1, file = file.path(out_dir, paste0(label,"_down_enrich.png")),
    device = "png",
    width = 8, height = 8,
    dpi = 400
)
