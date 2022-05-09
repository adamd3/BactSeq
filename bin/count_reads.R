#!/usr/bin/env Rscript

if (!require("optparse")){
    install.packages("optparse",repos = "http://cran.us.r-project.org")
}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos = "http://cran.us.r-project.org")
if (!require("Rsubread")){
    BiocManager::install("Rsubread")
}
if (!require("ape")){
    install.packages("ape",repos = "http://cran.us.r-project.org")
}
if (!require("stringr")){
    install.packages("stringr",repos = "http://cran.us.r-project.org")
}
if (!require("ggplot2")){
    install.packages("ggplot2",repos = "http://cran.us.r-project.org")
}
if (!require("scales")){
    install.packages("scales",repos = "http://cran.us.r-project.org")
}
if (!require("RColorBrewer")){
    install.packages("RColorBrewer",repos = "http://cran.us.r-project.org")
}

library(optparse)
library(Rsubread)
library(ape)
library(stringr)
library(ggplot2)
library(scales)
library(RColorBrewer)


option_list <- list(
    make_option(c("-m", "--metadata"), type="character", default=NULL,
        help="sample metadata tsv file", metavar="character"),
    make_option(c("-g", "--gff"), type="character", default=NULL,
        help="GFF annotation file for the reference strain", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

meta_f <- opt$metadata
gff_f <- opt$gff
outf <- opt$outf



##------------------------------------------------------------------------------
## Read data
##------------------------------------------------------------------------------
meta_tab <- read.table(meta_f, header = TRUE, sep = "\t")
## columns: sample	filename	group	repeat    path_to_file


## cat the counts files
# system("cat *.counts > merged_counts.txt")
total_counts_list <- lapply(meta_tab$sample, function(x){
    total_counts <- read.table(paste0(x,".counts"), header = FALSE, sep = "\t")
    total_counts$sample <- x
    colnames(total_counts) <- c("chr","chr_size","mapped","blank","sample")
    total_counts
})
merged_total_counts <- as.data.frame(do.call(rbind, total_counts_list))


##------------------------------------------------------------------------------
## Read genome annotation
##------------------------------------------------------------------------------
ref_annot <- ape::read.gff(gff_f, na.strings = c(".", "?"), GFF3 = TRUE)

ref_annot <- subset(ref_annot, type=="gene")
dim(ref_annot)
# [1] 4970    9

gene_attr <- stringr::str_split(ref_annot$attributes,";")
locus_tags <- unlist(lapply(gene_attr, function(x){x[grepl("locus_tag", x)]}))
gene_biotypes <- unlist(lapply(gene_attr, function(x){x[grepl("gene_biotype", x)]}))
common_gene_names <- unlist(lapply(gene_attr, function(x){
    x <- x[grepl("gene=", x)]
    x[identical(x, character(0))] <- ""
    x
}))
gene_lengths <- (ref_annot$end - ref_annot$start) + 1

ref_gene_df <- data.frame(
    locus_tag = locus_tags,
    biotype = gene_biotypes,
    gene_name = common_gene_names,
    gene_length = gene_lengths
)
ref_gene_df$locus_tag <- gsub("locus_tag=","",ref_gene_df$locus_tag)
ref_gene_df$biotype <- gsub("gene_biotype=","",ref_gene_df$biotype)
ref_gene_df$gene_name <- gsub("gene=","",ref_gene_df$gene_name)

write.table(
    ref_gene_df, "ref_gene_df.tsv", col.names = TRUE, row.names = FALSE,
    sep = "\t", quote = FALSE
)


##------------------------------------------------------------------------------
## Count reads mapping to genes
##------------------------------------------------------------------------------
bamfilesCount <- paste0(meta_tab$sample, ".bam")

gene_counts <- Rsubread::featureCounts(
    bamfilesCount, annot.ext = gff_f,
    isGTFAnnotationFile = TRUE,
    GTF.featureType = "gene", GTF.attrType = "locus_tag",
    nthreads = 16, countMultiMappingReads = TRUE,
    isPairedEnd = FALSE, strandSpecific = 2
)
colnames(gene_counts$counts) <- gsub(".bam", "", colnames(gene_counts$counts))
colnames(gene_counts$counts) <- gsub("\\.","_",colnames(gene_counts$counts))

write.table(
    gene_counts$counts, "gene_counts.tsv", col.names = TRUE, row.names = TRUE,
    sep = "\t", quote = FALSE
)



##------------------------------------------------------------------------------
## Plot library composition
##------------------------------------------------------------------------------
## set up plots
brewer_pallette1 <- brewer.pal(9,"Set1")
brewer_pallette3 <- brewer.pal(8,"Dark2")

gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}
ggColsDefault <- (gg_color_hue(4))
ggCols <- c(brewer_pallette1[1],brewer_pallette1[3])

## summarise counts per sample
counts_summary <- data.frame(
    sample = meta_tab$"sample",
    group = meta_tab$"group",
    rep = meta_tab$"rep",
    # protein_coding = colSums(
    #     gene_counts$counts[ref_gene_df$biotype=="protein_coding",]),
    # tRNA = colSums(
    #     gene_counts$counts[ref_gene_df$biotype=="tRNA",]),
    other_genes = colSums(
        gene_counts$counts[!(ref_gene_df$biotype=="rRNA"),]),
    rRNA_genes = colSums(
        gene_counts$counts[ref_gene_df$biotype=="rRNA",])
)


## order rows
# counts_summary <- counts_summary[order(
#     counts_summary$CRISPR_sytem, counts_summary$target, counts_summary$sample),]

## add total mapped counts
counts_summary <- merge(counts_summary,merged_total_counts, by = "sample")
counts_summary$non_rRNA_genes <- counts_summary$mapped-counts_summary$rRNA_genes

counts_summary <- counts_summary[rev(order(counts_summary$sample)),]

counts_melt <- reshape2::melt(
    counts_summary, id.vars = c("sample"),
    measure.vars = c("non_rRNA_genes", "rRNA_genes")
)
counts_melt$sample <- factor(
    counts_melt$sample, levels = rev(unique(sort(counts_melt$sample)))
)
counts_melt$variable <- factor(counts_melt$variable, levels=c("rRNA_genes", "non_rRNA_genes"))
# minUsable <- min(mergedDf$q15_dedup_reads)

cc1 <- 12
p1 <- ggplot(counts_melt,
        aes(x = sample, colour = variable, fill = variable, y = value)
    ) + geom_bar(position = "stack", stat = "identity", width = 0.7) +
    coord_flip() +
    xlab("Sample") + ylab("Million reads") +
    scale_fill_manual(
        "",
        # labels = c("Unaligned", "Low quality", "Duplicate", "Usable"),
        values = ggCols,
        guide = guide_legend(reverse = TRUE)
        # values = ggColsDefault
    ) +
    scale_colour_manual(values = ggCols, guide = FALSE) +
    scale_y_continuous(labels = unit_format(unit = "", scale = 1e-6)) +
    ## add a dashed line at the min usable number of reads
    geom_hline(yintercept = 5e6, linetype="dashed") +
    theme_bw(base_size = cc1*1.5) +
    theme(
        legend.position = "top",
        legend.title = element_blank(),
        legend.text=element_text(size=cc1*1.6),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")
    )

nsamps <- ncol(gene_counts$counts)

ggsave(
    p1, file = paste0('library_composition.png'),
    device = "png",
    width = 8, height = (nsamps/2),
    dpi = 300
)
