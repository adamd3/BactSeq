process FUNC_ENRICHMENT {
    tag "$gene_counts"
    label 'process_medium'
    publishDir "${params.outdir}/diff_expr", mode: 'copy'

    input:
    path ch_ann_file
    path ch_deseq_res

    output:
    path '*.tsv', emit: func_res
    path '*.png', emit: func_plots

    script:
    """
    for resf in ./DGE*tsv; do
        bname=$(basename $resf .tsv)
        functional_enrichment.R  \
            -r $resf -a $ch_ann_file  \
            -l $bname -o ./
    done
    """
}
