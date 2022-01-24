process TMM_NORMALISE_COUNTS {
    tag "$gene_counts"
    label 'process_medium'
    publishDir "${params.outdir}/read_counts", mode: 'copy'

    input:
    path gene_counts
    path ref_gene_df

    output:
    path 'cpm_counts.tsv', emit: cpm_counts
    path 'rpkm_counts.tsv', emit: rpkm_counts

    script:
    """
    TMM_normalise_counts.R -t TRUE -o ./
    """
}
