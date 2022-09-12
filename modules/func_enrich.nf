process FUNC_ENRICHMENT {
    tag "$ch_func_file"
    label 'process_high'
    publishDir "${params.outdir}/func_enrich", mode: 'copy'

    input:
    path ch_func_file
    path ch_deseq_res
    val p_thresh
    val l2fc_thresh

    output:
    path '*.tsv', emit: func_res
    path '*.png', emit: func_plots

    script:
    """
    functional_enrichment_topGO.R -a $ch_func_file -p $p_thresh -l $l2fc_thresh -o ./
    """
}
