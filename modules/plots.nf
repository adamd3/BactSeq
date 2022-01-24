process PCA_SAMPLES {
    tag "$tmm_counts"
    label 'process_medium'
    publishDir "${params.outdir}/PCA_samples", mode: 'copy'

    input:
    path tmm_counts
    path meta_merged

    output:
    path '*.{rds,png}', emit: pca_out

    script:
    """
    pca.R -o ./
    """
}
