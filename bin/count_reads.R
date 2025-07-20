#!/usr/bin/env Rscript
library(optparse)
library(Rsubread)
library(ape)
library(tidyverse)
library(scales)
library(RColorBrewer)


option_list <- list(
    make_option(c("-m", "--metadata"),
        type = "character", default = NULL,
        help = "sample metadata tsv file", metavar = "character"
    ),
    make_option(c("-g", "--gff"),
        type = "character", default = NULL,
        help = "GFF annotation file for the reference strain",
        metavar = "character"
    ),
    make_option(c("-p", "--is_paired"),
        type = "character", default = NULL,
        help = "are the reads paired-end? default = FALSE",
        metavar = "character"
    ),
    make_option(c("-s", "--strandedness"),
        type = "character", default = NULL,
        help = "read strandedness. default = reverse",
        metavar = "character"
    ),
    make_option(c("-t", "--threads"),
        type = "numeric", default = 1,
        help = "number of threads to use. default = 1.",
        metavar = "numeric"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

meta_f <- opt$metadata
gff_f <- opt$gff
ispaired <- if (opt$is_paired == "TRUE") TRUE else FALSE
strandedness <- opt$strandedness
threads <- opt$threads



## -----------------------------------------------------------------------------
## Functions
## -----------------------------------------------------------------------------
parse_gff_attributes <- function(
    attributes_string, attribute_name, start_pos, end_pos, default_val = NULL) { 
  attr_pairs_raw <- stringr::str_split(attributes_string, ";")[[1]]
  attr_pairs <- stringr::str_trim(
    stringr::str_replace_all(attr_pairs_raw, "[\\s\\p{Zs}\\p{C}]", " "))
  search_pattern_regex <- paste0("^", attribute_name, "=")
  target_attr <- attr_pairs[grepl(search_pattern_regex, attr_pairs)] 
  if (length(target_attr) > 0) {
    value <- stringr::str_remove(target_attr[1], search_pattern_regex)
    return(value)
  } else {
    if (!is.null(default_val)) {
      return(default_val)
    } else {
      default_value_generated <- paste0("unknown_", start_pos, "_", end_pos)
      return(default_value_generated)
    }
  }
}


## -----------------------------------------------------------------------------
## Read data
## -----------------------------------------------------------------------------
meta_tab <- read_tsv(meta_f)

total_counts_list <- map(meta_tab$sample, function(x) {
    mapped_count <- read_tsv(paste0(x, ".counts"), col_names = "mapped")
    mapped_count$sample <- x
    mapped_count
})
merged_total_counts <- bind_rows(total_counts_list)



## -----------------------------------------------------------------------------
## Read genome annotation
## -----------------------------------------------------------------------------
ref_annot <- ape::read.gff(gff_f, na.strings = c(".", "?"), GFF3 = TRUE)

ref_annot <- subset(ref_annot, type == "gene")

gene_attr <- stringr::str_split(ref_annot$attributes, ";")

locus_tags <- unname(mapply(parse_gff_attributes,
                     ref_annot$attributes,
                     attribute_name = "locus_tag",
                     start_pos = ref_annot$start,
                     end_pos = ref_annot$end))

gene_biotypes <- unname(mapply(parse_gff_attributes,
                    ref_annot$attributes,
                    attribute_name = "gene_biotype",
                    start_pos = ref_annot$start,
                    end_pos = ref_annot$end,
                    "unknown"))

common_gene_names <- unname(mapply(parse_gff_attributes,
                    ref_annot$attributes,
                    attribute_name = "gene",
                    start_pos = ref_annot$start,
                    end_pos = ref_annot$end,
                    "unknown"))


gene_lengths <- (ref_annot$end - ref_annot$start) + 1

ref_gene_df <- data.frame(
    locus_tag = locus_tags,
    biotype = gene_biotypes,
    gene_name = common_gene_names,
    gene_length = gene_lengths
)
ref_gene_df$locus_tag <- sub("^.*=", "", ref_gene_df$locus_tag)
ref_gene_df$biotype <- sub("^.*=", "", ref_gene_df$biotype)
ref_gene_df$gene_name <- sub("^.*=", "", ref_gene_df$gene_name)

write_tsv(ref_gene_df, "ref_gene_df.tsv")


## -----------------------------------------------------------------------------
## Count reads mapping to genes
## -----------------------------------------------------------------------------
bamfilesCount <- paste0(meta_tab$sample, ".bam")

# 0 (unstranded), 1 (stranded) and 2 (reversely stranded)
strand <- switch(as.character(strandedness),
    "unstranded" = 0,
    "forward" = 1,
    "reverse" = 2,
    stop("Invalid input")
)

gene_counts <- Rsubread::featureCounts(
    bamfilesCount,
    annot.ext = gff_f,
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
colnames(gene_counts$counts) <- gsub("\\.", "_", colnames(gene_counts$counts))


counts_mat <- gene_counts$counts %>%
    as_tibble(rownames = "feature_id")

write_tsv(counts_mat, "gene_counts.tsv")

## protein-coding genes only
gene_counts_pc <- counts_mat[ref_gene_df$biotype == "protein_coding", ]
write_tsv(gene_counts_pc, "gene_counts_pc.tsv")


## -----------------------------------------------------------------------------
## Plot library composition
## -----------------------------------------------------------------------------
## set up plots
brewer_pallette1 <- brewer.pal(9, "Set1")
brewer_pallette3 <- brewer.pal(8, "Dark2")

gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}
ggColsDefault <- (gg_color_hue(4))
ggCols <- brewer_pallette1[c(1, 3, 4, 5, 2, 7, 8)]

## summarise counts per sample

all_biotypes <- unique(ref_gene_df$biotype)

biotype_counts <- data.frame(do.call(
    cbind,
    lapply(all_biotypes, function(biotype) {
        colSums(gene_counts$counts[ref_gene_df$biotype == biotype, , drop = FALSE])
    })
))
colnames(biotype_counts) <- all_biotypes



counts_summary <- data.frame(
    sample = meta_tab$"sample",
    group = meta_tab$"group",
    rep = meta_tab$"rep_no"
)

counts_summary$rRNA <- biotype_counts$rRNA

## add total mapped counts
counts_summary <- merge(counts_summary, merged_total_counts, by = "sample")

counts_summary$other <- counts_summary$mapped - counts_summary$rRNA

write_tsv(counts_summary, "counts_summary.tsv")


counts_summary <- counts_summary[rev(order(counts_summary$sample)), ]
all_biotypes <- c(all_biotypes, "other")
non_rRNA_btypes <- all_biotypes[!all_biotypes == "rRNA"]


#############################
## raw counts plot
#############################
counts_melt <- counts_summary %>%
    pivot_longer(
        cols = c("other", "rRNA"),
        names_to = "variable",
        values_to = "value"
    ) %>%
    mutate(
        sample = factor(sample, levels = rev(unique(sort(sample)))),
        variable = factor(variable, levels = c("rRNA", "other"))
    )
# minUsable <- min(mergedDf$q15_dedup_reads)


ylabel <- ifelse(isTRUE(ispaired), "Million read pairs", "Million reads")

p1 <- ggplot(
    counts_melt,
    aes(x = sample, colour = variable, fill = variable, y = value)
) +
    geom_bar(position = "stack", stat = "identity", width = 0.7) +
    coord_flip() +
    xlab("Sample") +
    ylab(ylabel) +
    scale_fill_manual(
        "",
        values = ggCols,
        guide = guide_legend(reverse = TRUE)
    ) +
    scale_colour_manual(values = ggCols, guide = FALSE) +
    scale_y_continuous(labels = unit_format(unit = "", scale = 1e-6)) +
    ## add a dashed line at the min usable number of reads
    # geom_hline(yintercept = 5e6, linetype="dashed") +
    theme_bw(base_size = 15) +
    theme(
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")
    )

nsamps <- ncol(gene_counts$counts)

ggsave(
    p1,
    file = paste0("library_composition.png"),
    device = "png",
    width = 8, height = (nsamps / 2.2),
    dpi = 300
)


#############################
## proportions plot
#############################
## get the proportions of reads per library
# propCols <- (mergedDf[,c(3,13,14,5)] / mergedDf[,2])

propCols <- counts_summary %>%
    mutate(
        other = other / mapped,
        rRNA = rRNA / mapped
    ) %>%
    select(sample, other, rRNA)

prop_melt <- propCols %>%
    pivot_longer(
        cols = c("other", "rRNA"),
        names_to = "variable",
        values_to = "value"
    ) %>%
    mutate(
        sample = factor(sample, levels = rev(unique(sort(sample)))),
        variable = factor(variable, levels = c("rRNA", "other"))
    )

p2 <- ggplot(
    prop_melt,
    aes(x = sample, colour = variable, fill = variable, y = value)
) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_flip() +
    xlab("Sample") +
    ylab("Proportion of reads") +
    scale_fill_manual(
        "",
        values = ggCols,
        guide = guide_legend(reverse = TRUE)
    ) +
    scale_colour_manual(values = ggCols, guide = FALSE) +
    scale_y_continuous(labels = comma) +
    theme_bw(base_size = 15) +
    theme(
        legend.position = "top",
        legend.text = element_text(size = 12),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")
    )

ggsave(
    p2,
    file = "library_composition_proportions.png",
    device = "png",
    width = 8, height = (nsamps / 2.2),
    dpi = 300
)
