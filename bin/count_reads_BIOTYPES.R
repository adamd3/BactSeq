#!/usr/bin/env Rscript

if (!require("optparse")){
    install.packages("optparse",repos = "http://cran.us.r-project.org")
}
if (!require("tidyverse")){
    install.packages("tidyverse",repos = "http://cran.us.r-project.org")
}
if (!require("scales")){
    install.packages("scales",repos = "http://cran.us.r-project.org")
}
if (!require("RColorBrewer")){
    install.packages("RColorBrewer",repos = "http://cran.us.r-project.org")
} 


library(optparse)
library(tidyverse)
library(scales)
library(RColorBrewer)
library(reshape2)


option_list <- list(
    make_option(c("-m", "--metadata"), type="character", default=NULL,
        help="sample metadata tsv file", metavar="character"),
    make_option(c("-r", "--ref_gene_f"), type="character", default=NULL,
        help="Gene annotations in reference strain", 
        metavar="character"),
    make_option(c("-r", "--gene_counts_f"), type="character", default=NULL,
        help="Gene annotations in reference strain", 
        metavar="character"),
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

meta_f <- opt$metadata
ref_gene_f <- opt$ref_gene_f
gene_counts_f <- opt$gene_counts_f
outf <- opt$outf



##------------------------------------------------------------------------------
## Read data
##------------------------------------------------------------------------------
meta_tab <- read.table(meta_f, header = TRUE, sep = "\t")
## columns: sample	 file1   file2	group	rep_no    paired


## cat the counts files
# system("cat *.counts > merged_counts.txt")

total_counts_list <- lapply(meta_tab$sample, function(x){
    # total_counts <- read.table(paste0(x,".counts"), header = FALSE, sep = "\t")
    # total_counts$sample <- x
    # colnames(total_counts) <- c("chr","chr_size","mapped","blank","sample")
    # total_counts
    mapped_count <- read.table(paste0(x,".counts"), header = FALSE)
    colnames(mapped_count) <- "mapped"
    mapped_count$sample <- x
    mapped_count
})
merged_total_counts <- as.data.frame(do.call(rbind, total_counts_list))



##------------------------------------------------------------------------------
## Read genome annotation
##------------------------------------------------------------------------------
# ref_annot <- ape::read.gff(gff_f, na.strings = c(".", "?"), GFF3 = TRUE)

# ref_annot <- subset(ref_annot, type=="gene")

# gene_attr <- stringr::str_split(ref_annot$attributes,";")
# locus_tags <- unlist(lapply(gene_attr, function(x){x[grepl("locus_tag", x)]}))
# gene_biotypes <- unlist(lapply(gene_attr, function(x){x[grepl("gene_biotype", x)]}))
# common_gene_names <- unlist(lapply(gene_attr, function(x){
#     x <- x[grepl("gene=", x)]
#     x[identical(x, character(0))] <- ""
#     x
# }))
# gene_lengths <- (ref_annot$end - ref_annot$start) + 1

# ref_gene_df <- data.frame(
#     locus_tag = locus_tags,
#     biotype = gene_biotypes,
#     gene_name = common_gene_names,
#     gene_length = gene_lengths
# )
# ref_gene_df$locus_tag <- gsub("locus_tag=","",ref_gene_df$locus_tag)
# ref_gene_df$biotype <- gsub("gene_biotype=","",ref_gene_df$biotype)
# ref_gene_df$gene_name <- gsub("gene=","",ref_gene_df$gene_name)

# write.table(
#     ref_gene_df, "ref_gene_df.tsv", col.names = TRUE, row.names = FALSE,
#     sep = "\t", quote = FALSE
# )

ref_gene_df <- read_tsv(ref_gene_f)


##------------------------------------------------------------------------------
## Count reads mapping to genes
##------------------------------------------------------------------------------
# bamfilesCount <- paste0(meta_tab$sample, ".bam")

# # 0 (unstranded), 1 (stranded) and 2 (reversely stranded)
# strand <- switch(as.character(strandedness),
#        "unstranded" = 0,
#        "forward" = 1,
#        "reverse" = 2,
#        stop("Invalid input")
# )

# gene_counts <- Rsubread::featureCounts(
#     bamfilesCount, annot.ext = gff_f,
#     isGTFAnnotationFile = TRUE,
#     GTF.featureType = "gene", 
#     GTF.attrType = "locus_tag",
#     nthreads = threads, 
#     countMultiMappingReads = TRUE, 
#     fraction = TRUE, ## assign fractional counts to multimappers
#     isPairedEnd = ispaired, 
#     strandSpecific = strand
# )
# colnames(gene_counts$counts) <- gsub(".bam", "", colnames(gene_counts$counts))
# colnames(gene_counts$counts) <- gsub("\\.","_",colnames(gene_counts$counts))


# counts_mat <- gene_counts$counts
# counts_mat <- tibble::rownames_to_column(as.data.frame(counts_mat), "feature_id")


# write.table(
#     counts_mat, "gene_counts.tsv", col.names = TRUE, row.names = FALSE,
#     sep = "\t", quote = FALSE
# )

# ## protein-coding genes only
# gene_counts_pc <- counts_mat[ref_gene_df$biotype=="protein_coding",]
# write.table(
#     gene_counts_pc, "gene_counts_pc.tsv", col.names = TRUE, row.names = FALSE,
#     sep = "\t", quote = FALSE
# )

gene_counts <- read_tsv(gene_counts_f)

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
ggCols <- brewer_pallette1[c(1,3,4,5,2,7,8)]

## summarise counts per sample

all_biotypes <- unique(ref_gene_df$biotype)

biotype_counts <- data.frame(do.call(cbind,
    lapply(all_biotypes, function(biotype){
    colSums(gene_counts$counts[ref_gene_df$biotype==biotype,,drop=FALSE])
})))
colnames(biotype_counts) <- all_biotypes



counts_summary <- data.frame(
    sample = meta_tab$"sample",
    group = meta_tab$"group",
    rep = meta_tab$"rep_no"#,
    # protein_coding = colSums(
    #     gene_counts$counts[ref_gene_df$biotype=="protein_coding",]),
    # tRNA = colSums(
        # gene_counts$counts[ref_gene_df$biotype=="tRNA",]),
    # rRNA = colSums(
    #     gene_counts$counts[ref_gene_df$biotype=="rRNA",])
)

counts_summary <- cbind(counts_summary,biotype_counts)

## add total mapped counts
counts_summary <- merge(counts_summary,merged_total_counts, by = "sample")
# counts_summary$other <- counts_summary$mapped-counts_summary$rRNA


counts_summary$antisense_other <- counts_summary$mapped - rowSums(biotype_counts)#(
    #counts_summary$rRNA + counts_summary$protein_coding)

write.table(
    counts_summary, "counts_summary.tsv", col.names = TRUE, row.names = FALSE,
    sep = "\t", quote = FALSE
)

write.table(
    biotype_counts, "biotype_counts.tsv", col.names = TRUE, row.names = FALSE,
    sep = "\t", quote = FALSE
)



counts_summary <- counts_summary[rev(order(counts_summary$sample)),]


cc1 <- 12

all_biotypes <- c(all_biotypes, "antisense_other")
non_rRNA_btypes <- all_biotypes[!all_biotypes=="rRNA"]


#############################
## raw counts plot
#############################
counts_melt <- reshape2::melt(
    counts_summary, id.vars = c("sample"),
    measure.vars = c(
        # "protein_coding", 
        # "antisense_other",
        # "rRNA"
        all_biotypes
        )
)
counts_melt$sample <- factor(
    counts_melt$sample, levels = rev(unique(sort(counts_melt$sample)))
)
counts_melt$variable <- factor(counts_melt$variable, levels=c(
    "rRNA", 
    non_rRNA_btypes
    # "antisense_other", 
    # "protein_coding"
    )
)
# minUsable <- min(mergedDf$q15_dedup_reads)


ylabel <- ifelse(isTRUE(ispaired), "Million read pairs", "Million reads")

p1 <- ggplot(counts_melt,
        aes(x = sample, colour = variable, fill = variable, y = value)
    ) + 
    geom_bar(position = "stack", stat = "identity", width = 0.7) +
    coord_flip() +
    xlab("Sample") + ylab(ylabel) +
    scale_fill_manual(
        "",
        values = ggCols,
        guide = guide_legend(reverse = TRUE)
    ) +
    scale_colour_manual(values = ggCols, guide = FALSE) +
    scale_y_continuous(labels = unit_format(unit = "", scale = 1e-6)) +
    ## add a dashed line at the min usable number of reads
    # geom_hline(yintercept = 5e6, linetype="dashed") +
    theme_bw(base_size = cc1*1.3) +
    theme(
        legend.position = "top",
        legend.title = element_blank(),
        legend.text=element_text(size=cc1),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")
    )

nsamps <- ncol(gene_counts$counts)

ggsave(
    p1, file = paste0('library_composition.png'),
    device = "png",
    width = 8, height = (nsamps/2.2),
    dpi = 300
)


#############################
## proportions plot
#############################
## get the proportions of reads per library
# propCols <- (mergedDf[,c(3,13,14,5)] / mergedDf[,2])

propCols <- counts_summary[c(
    # "protein_coding", 
    # "antisense_other", 
    # "rRNA"
    all_biotypes
    )] / counts_summary$mapped

propCols$sample <- counts_summary$sample

# rowSums(propCols) ## each row should sum to 1
prop_melt <- melt(
    propCols, id.vars = c("sample"),
    measure.vars = c(
        # "protein_coding", 
        # "antisense_other", 
        # "rRNA"
        all_biotypes
        )
)
prop_melt$sample <- factor(
    prop_melt$sample, levels = rev(unique(sort(prop_melt$sample)))
)
prop_melt$variable <- factor(prop_melt$variable, levels=c(
    "rRNA", 
    non_rRNA_btypes
    # "antisense_other", 
    # "protein_coding"
    )
)

p2 <- ggplot(prop_melt,
         aes(x = sample, colour = variable, fill = variable, y = value)
    ) + geom_bar(stat = "identity",  width = 0.7) +
    coord_flip()  +
    xlab("Sample") + ylab("Proportion of reads") +
    scale_fill_manual(
        "",
        values = ggCols,
        guide = guide_legend(reverse = TRUE)
    ) +
    scale_colour_manual(values = ggCols, guide = FALSE) +
    scale_y_continuous(labels = comma) +
    theme_bw(base_size = cc1*1.3) +
    theme(
        legend.position="top", 
        legend.text=element_text(size=cc1),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")
        )

ggsave(
    p2, file = 'library_composition_proportions.png',
    device = "png",
    width = 8, height = (nsamps/2.2),
    dpi = 300
)