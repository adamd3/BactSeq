process TMM_NORMALISE_COUNTS {
    tag "$gene_counts"
    label 'process_medium'
    publishDir "${params.outdir}/read_counts", mode: 'copy'

    input:
    path gene_counts

    output:
    path 'norm_counts.tsv', emit: tmm_counts

    script:
    """
    TMM_normalise_counts.R -t TRUE -o norm_counts.tsv
    """
}
