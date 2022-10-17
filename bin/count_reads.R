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
if (!require("tibble")){
    install.packages("tibble",repos = "http://cran.us.r-project.org")
}


library(optparse)
library(Rsubread)
library(ape)
library(stringr)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(reshape2)
library(tibble)


option_list <- list(
    make_option(c("-m", "--metadata"), type="character", default=NULL,
        help="sample metadata tsv file", metavar="character"),
    make_option(c("-g", "--gff"), type="character", default=NULL,
        help="GFF annotation file for the reference strain", 
        metavar="character"),
    make_option(c("-p", "--is_paired"), type="character", default=NULL,
        help="are the reads paired-end? default = FALSE", 
        metavar="character"),
    make_option(c("-s", "--strandedness"), type="character", default=NULL,
        help="read strandedness. default = reverse", 
        metavar="character"),
    make_option(c("-t", "--threads"), type="numeric", default=1,
        help="number of threads to use. default = 1.", 
        metavar="numeric")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

meta_f <- opt$metadata
gff_f <- opt$gff
ispaired <- if(opt$is_paired == "TRUE") TRUE else FALSE
strandedness <- opt$strandedness
threads <- opt$threads
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
ref_annot <- ape::read.gff(gff_f, na.strings = c(".", "?"), GFF3 = TRUE)

ref_annot <- subset(ref_annot, type=="gene")

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

# 0 (unstranded), 1 (stranded) and 2 (reversely stranded)
strand <- switch(as.character(strandedness),
       "unstranded" = 0,
       "forward" = 1,
       "reverse" = 2,
       stop("Invalid strandedness")
)

gene_counts <- Rsubread::featureCounts(
    bamfilesCount, annot.ext = gff_f,
    isGTFAnnotationFile = TRUE,
    GTF.featureType = "gene", 
    GTF.attrType = "locus_tag",
    nthreads = threads, 
    countMultiMappingReads = TRUE, 
    fraction = TRUE, ## assign fractional counts to multimappers
    isPairedEnd = ispaired, 
    strandSpecific = strand
)
colnames(gene_counts$counts) <- gsub(".bam", "", colnames(gene_counts$counts))
colnames(gene_counts$counts) <- gsub("\\.","_",colnames(gene_counts$counts))


counts_mat <- gene_counts$counts
counts_mat <- tibble::rownames_to_column(as.data.frame(counts_mat), "feature_id")


write.table(
    counts_mat, "gene_counts.tsv", col.names = TRUE, row.names = FALSE,
    sep = "\t", quote = FALSE
)

## protein-coding genes only
gene_counts_pc <- counts_mat[ref_gene_df$biotype=="protein_coding",]
write.table(
    gene_counts_pc, "gene_counts_pc.tsv", col.names = TRUE, row.names = FALSE,
    sep = "\t", quote = FALSE
)


## count antisense (if libraries are stranded)
if (!strand==0){
    as_strand <- switch(as.character(strand),
       "1" = 2,
       "2" = 1,
       stop("Invalid strandedness")
    )
    as_gene_counts <- Rsubread::featureCounts(
        bamfilesCount, annot.ext = gff_f,
        isGTFAnnotationFile = TRUE,
        GTF.featureType = "gene", 
        GTF.attrType = "locus_tag",
        nthreads = threads, 
        countMultiMappingReads = TRUE, 
        fraction = TRUE, ## assign fractional counts to multimappers
        isPairedEnd = ispaired, 
        strandSpecific = as_strand
    )
    colnames(as_gene_counts$counts) <- gsub(".bam", "", colnames(as_gene_counts$counts))
    colnames(as_gene_counts$counts) <- gsub("\\.","_",colnames(as_gene_counts$counts))
    as_counts_mat <- as_gene_counts$counts
} else {
    ## set antisense counts to 0
    as_counts_mat <- gene_counts$counts
    as_counts_mat[!as_counts_mat==0] <- 0 
}

##------------------------------------------------------------------------------
## Plot library composition
##------------------------------------------------------------------------------
## set up plots
brewer_pallette1 <- brewer.pal(9,"Set1")
brewer_pallette3 <- brewer.pal(8,"Dark2")

gg_color_hue <- function(n){
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}
ggColsDefault <- (gg_color_hue(4))
ggCols <- c(brewer_pallette1[c(1,3,4,5,2,7,8)],"grey",brewer_pallette3)

## summarise counts per sample

all_biotypes <- unique(ref_gene_df$biotype)

biotype_counts <- data.frame(do.call(cbind,
    lapply(all_biotypes, function(biotype){
    colSums(gene_counts$counts[ref_gene_df$biotype==biotype,,drop=FALSE])
})))
colnames(biotype_counts) <- all_biotypes

## add antisense counts 
biotype_counts[["antisense"]] <- colSums(as_counts_mat)

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


counts_summary$other <- counts_summary$mapped - rowSums(biotype_counts)#(
    #counts_summary$rRNA + counts_summary$protein_coding)

counts_summary <- counts_summary[rev(order(counts_summary$sample)),]


cc1 <- 12

all_biotypes <- c(all_biotypes, "antisense", "other")
non_rRNA_btypes <- all_biotypes[!all_biotypes=="rRNA"]


#############################
## raw counts plot
#############################
counts_melt <- reshape2::melt(
    counts_summary, id.vars = c("sample"),
    measure.vars = c(
        # "protein_coding", 
        # "other",
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
    # "other", 
    # "protein_coding"
    )
)
# minUsable <- min(mergedDf$q15_dedup_reads)


p1 <- ggplot(counts_melt,
        aes(x = sample, colour = variable, fill = variable, y = value)
    ) + 
    geom_bar(position = "stack", stat = "identity", width = 0.7) +
    coord_flip() +
    xlab("Sample") + ylab("Million reads") +
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
    # "other", 
    # "rRNA"
    all_biotypes
    )] / counts_summary$mapped

propCols$sample <- counts_summary$sample

# rowSums(propCols) ## each row should sum to 1
prop_melt <- melt(
    propCols, id.vars = c("sample"),
    measure.vars = c(
        # "protein_coding", 
        # "other", 
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
    # "other", 
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