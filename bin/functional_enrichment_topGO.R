#!/usr/bin/env Rscript

library(optparse)
library(topGO)
library(GO.db)
library(readr)
library(ggplot2)
library(stringr)
library(dplyr)

option_list <- list(
    make_option(c("-a", "--annot_f"),
        type = "character", default = NULL,
        help = "CSV file containing functional annotations",
        metavar = "character"
    ),
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
    make_option(c("-n", "--min_size"),
        type = "double", default = 10,
        help = "minimum size of the functional categories to test (default = 10)",
        metavar = "character"
    ),
    make_option(c("-o", "--out_dir"),
        type = "character", default = NULL,
        help = "output dir for results", metavar = "character"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

annot_f <- opt$annot_f
p_threshold <- opt$p_threshold
log2fc_threshold <- opt$log2fc_threshold
min_size <- opt$min_size
out_dir <- opt$out_dir

## colours for barplots
up_col <- "#E41A1C"
down_col <- "#7FC97F"


## ------------------------------------------------------------------------------
## Read + process data
## ------------------------------------------------------------------------------
dge_files <- dir(path = "./", pattern = "^DGE_.*\\.tsv$")


in_dat <- read_csv(annot_f, show_col_types = FALSE)
colnames(in_dat) <- c("Gene", "GO_terms")


# ## Df of GO term annotations
# term_list <- lapply(in_dat[["GO_terms"]], function(terms){
#     str_split(terms,",")
# })
#
# go_df <- data.frame(GO.ID = unique(unlist(term_list)))
# go_df[["Term"]] <- as.character(sapply(go_df[["GO.ID"]], Term))


## make a list of gene -> go terms
gene_terms_list <- lapply(1:nrow(in_dat), function(idx) {
    gene <- in_dat[idx, "Gene"]
    terms <- in_dat[idx, "GO_terms"]
    terms <- str_split(terms, ",")
    terms[[1]]
})
names(gene_terms_list) <- in_dat[["Gene"]]

## make a list of go term -> genes
all_terms <- unique(unlist(gene_terms_list))

term_genes_list <- lapply(all_terms, function(term) {
    term_presence <- sapply(gene_terms_list, function(gene_terms) {
        term %in% gene_terms
    })
    names(gene_terms_list)[term_presence]
})
names(term_genes_list) <- all_terms




## ------------------------------------------------------------------------------
## Perform enrichment testing
## ------------------------------------------------------------------------------

lapply(dge_files, function(res_f) {
    # print(res_f)

    res_tab <- na.omit(read_tsv(res_f))

    all_genes <- res_tab[["feature_id"]]


    label <- sub(
        pattern = "(.*)\\..*$", replacement = "\\1", basename(file.path(res_f)))

    ## get up/down-regulated genes
    padj_col <- colnames(res_tab)[grepl("padj", colnames(res_tab))]
    l2fc_col <- colnames(res_tab)[grepl("log2FoldChange", colnames(res_tab))]

    genes_up <- (res_tab[(res_tab[[padj_col]] < p_threshold &
        res_tab[[l2fc_col]] > log2fc_threshold), ])[["feature_id"]]
    genes_down <- (res_tab[(res_tab[[padj_col]] < p_threshold &
        res_tab[[l2fc_col]] < (log2fc_threshold * -1)), ])[["feature_id"]]

    ## enrichment testing is only performed if >=5 DEGs identified
    if (length(genes_up) >= 5) {
        term_genes_list_up <- lapply(term_genes_list, function(x) {
            x[x %in% genes_up]
        })
        term_genes_list_up <- term_genes_list_up[lapply(
            term_genes_list_up, length
        ) > 0]

        ## make factor coding for upregulated genes
        upList <- factor(as.integer(all_genes %in% genes_up))
        names(upList) <- all_genes

        GOData_up <- new(
            "topGOdata",
            ontology = "BP", allGenes = upList, nodeSize = min_size,
            annot = annFUN.gene2GO, gene2GO = gene_terms_list
        )

        resultsFisher_up <- runTest(
            GOData_up,
            algorithm = "weight01", statistic = "fisher"
        )

        allRes_up <- GenTable(
            GOData_up,
            p_value = resultsFisher_up,
            topNodes = length(usedGO(GOData_up))
        )

        ## annotate which genes are in each Term
        allRes_up <- data.frame(do.call(rbind, lapply(
            1:nrow(allRes_up), function(idx) {
                row <- allRes_up[idx, ]
                term <- row[["GO.ID"]]
                no_genes <- length(term_genes_list_up[[term]])
                row[["no_genes"]] <- no_genes
                row[["genes"]] <- paste(term_genes_list_up[[term]], collapse = ";")
                row[["Significant"]] <- NULL
                row
            }
        )))

        write.table(
            allRes_up, file.path(out_dir, paste0(label, "_up_enrich.tsv")),
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE
        )

        ## barplot of top 20 terms
        go_up_plot <- allRes_up[1:20, ]
        go_up_plot[["order"]] <- factor(
            1:nrow(go_up_plot),
            levels = rev(1:nrow(go_up_plot))
        )
        go_up_plot[["p_value"]] <- as.numeric(go_up_plot[["p_value"]])

        p1 <- ggplot(
            go_up_plot,
            aes(x = order, y = -log10(p_value))
        ) +
            geom_bar(
                position = "dodge", stat = "identity",
                width = 0.7, colour = "black", fill = up_col
            ) +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
            theme_bw() +
            ylab(expression("−log"[10] * "(" * italic(P) ~ "value" * ")")) +
            scale_x_discrete(labels = rev(go_up_plot$Term)) +
            coord_flip() +
            theme(
                legend.position = "none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                axis.title.x = element_text(colour = "black", size = 12 * 1.25),
                axis.title.y = element_blank(),
                axis.text.x = element_text(colour = "black", size = 12 * 1.25),
                axis.text.y = element_text(colour = "black", size = 12 * 1.25),
                panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
                legend.text = element_text(colour = "black", size = 12 * 1.25),
                legend.title = element_text(colour = "black", size = 12 * 1.25)
            )

        ggsave(
            p1,
            file = file.path(out_dir, paste0(label, "_up_enrich.png")),
            device = "png",
            width = 8, height = 8,
            dpi = 400
        )
    }

    if (length(genes_down) >= 5) {
        term_genes_list_down <- lapply(term_genes_list, function(x) {
            x[x %in% genes_down]
        })
        term_genes_list_down <- term_genes_list_down[lapply(
            term_genes_list_down, length
        ) > 0]

        ## make factor coding for downregulated genes
        downList <- factor(as.integer(all_genes %in% genes_down))
        names(downList) <- all_genes

        GOData_down <- new(
            "topGOdata",
            ontology = "BP", allGenes = downList, nodeSize = min_size,
            annot = annFUN.gene2GO, gene2GO = gene_terms_list
        )

        resultsFisher_down <- runTest(
            GOData_down,
            algorithm = "weight01", statistic = "fisher"
        )
        allRes_down <- GenTable(
            GOData_down,
            p_value = resultsFisher_down,
            topNodes = length(usedGO(GOData_down))
        )

        ## annotate which genes are in each Term
        allRes_down <- data.frame(do.call(rbind, lapply(1:nrow(allRes_down), function(idx) {
            row <- allRes_down[idx, ]
            term <- row[["GO.ID"]]
            no_genes <- length(term_genes_list_down[[term]])
            row[["no_genes"]] <- no_genes
            row[["genes"]] <- paste(term_genes_list_down[[term]], collapse = ";")
            row[["Significant"]] <- NULL
            row
        })))

        write.table(
            allRes_down, file.path(out_dir, paste0(label, "_down_enrich.tsv")),
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE
        )


        ## barplots of top 20 terms
        go_down_plot <- allRes_down[1:20, ]
        go_down_plot[["order"]] <- factor(
            1:nrow(go_down_plot),
            levels = rev(1:nrow(go_down_plot))
        )
        go_down_plot[["p_value"]] <- as.numeric(go_down_plot[["p_value"]])

        p1 <- ggplot(
            go_down_plot,
            aes(x = order, y = -log10(p_value))
        ) +
            geom_bar(
                position = "dodge", stat = "identity",
                width = 0.7, colour = "black", fill = down_col
            ) +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
            theme_bw() +
            ylab(expression("−log"[10] * "(" * italic(P) ~ "value" * ")")) +
            scale_x_discrete(labels = rev(go_down_plot$Term)) +
            coord_flip() +
            theme(
                legend.position = "none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                plot.title = element_text(hjust = 0.5, size = 12 * 1.25),
                axis.title.x = element_text(colour = "black", size = 12 * 1.25),
                axis.title.y = element_blank(),
                axis.text.x = element_text(colour = "black", size = 12 * 1.25),
                axis.text.y = element_text(colour = "black", size = 12 * 1.25),
                panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
                legend.text = element_text(colour = "black", size = 12 * 1.25),
                legend.title = element_text(colour = "black", size = 12 * 1.25)
            )

        ggsave(
            p1,
            file = file.path(out_dir, paste0(label, "_down_enrich.png")),
            device = "png",
            width = 8, height = 8,
            dpi = 400
        )
    }
})
