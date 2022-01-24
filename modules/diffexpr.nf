process DIFF_EXPRESSION {
    tag "$gene_counts"
    label 'process_high'
    publishDir "${params.outdir}/diff_expr", mode: 'copy'

    input:
    path gene_counts
    path ref_gene_df

    output:
    path '*.tsv', emit: deseq_out

    script:
    """
    diffexpr.R -t TRUE -o ./
    """
}
