process COUNT_READS {
    tag "$gff"
    label 'process_high'
    publishDir "${params.outdir}/read_counts", mode: 'copy'

    input:
    path bam
    path bai
    path counts
    path meta
    path gff

    output:
    path 'gene_counts.tsv', emit: counts_df
    path 'ref_gene_df.tsv', emit: ref_gene_df
    path 'library_composition.png', emit: libcomp_plot

    script:
    """
    if (meta.paired_end) {
        """
        count_reads.R -m $meta -g $gff -p TRUE
        """
    } else {
        """
        count_reads.R -m $meta -g $gff -p FALSE
        """
    }

    """
}
