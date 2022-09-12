process DIFF_EXPRESSION {
    tag "$gene_counts"
    label 'process_high'
    publishDir "${params.outdir}/diff_expr", mode: 'copy'

    input:
    path gene_counts
    path meta_merged
    path cont_tabl
    val p_thresh
    val l2fc_thresh

    output:
    path '*.tsv', emit: deseq_res
    path '*.png', emit: deseq_volcano

    script:
    """
    diffexpr.R -p $p_thresh -l $l2fc_thresh -o ./ 
    """
}
