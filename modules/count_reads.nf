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
    val paired
    val strandedness

    output:
    path 'gene_counts.tsv', emit: counts_df
    path 'gene_counts_pc.tsv', emit: counts_df_pc
    path 'ref_gene_df.tsv', emit: ref_gene_df
    path 'library_composition.png', emit: libcomp_plot

    script:

    if (paired) {
        """
        count_reads.R -p TRUE -s $strandedness -m $meta -g $gff -p TRUE
        """
    } else {
        """
        count_reads.R -p FALSE -s $strandedness -m $meta -g $gff -p FALSE
        """
    }
}
