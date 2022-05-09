process FUNC_ENRICHMENT {
    tag "$gene_counts"
    label 'process_medium'
    publishDir "${params.outdir}/func_enrich", mode: 'copy'

    input:
    path ch_func_file
    path ch_deseq_res

    output:
    path '*.tsv', emit: func_res
    path '*.png', emit: func_plots

    script:
    """
    for resf in ./DGE*tsv; do
        functional_enrichment.R  \
            -r $resf -a $ch_func_file -o ./
    done
    """
}
