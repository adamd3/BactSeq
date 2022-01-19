process SUBSET_GENES {
    tag "$st_file"
    label 'process_medium'
    publishDir "${params.outdir}/gene_counts", mode: 'copy'

    input:
    path gpa_file
    path meta_merged
    path st_file
    val perc

    output:
    path 'gene_set_ST.tsv', emit: gene_subset

    script:
    """
    subset_genes.py \
        --gene_presence_absence=$gpa_file \
        --metadata_merged=$meta_merged \
        --ST_file=$st_file --perc=$perc \
        --outf=gene_set_ST.tsv
    """
}

process LENGTH_SCALE_COUNTS {
    tag "$name"
    label 'process_medium'
    publishDir "${params.outdir}/gene_counts", mode: 'copy'

    input:
    tuple val(name), path(merged_counts)
    tuple val(name), path(merged_lens)
    path gene_subset

    output:
    path 'kallisto_scaled_counts.tsv', emit: scaled_counts

    script:
    """
    length_scale_counts.R \
        -c $merged_counts \
        -l $merged_lens \
        -g $gene_subset \
        -i TRUE -p TRUE \
        -o kallisto_scaled_counts.tsv
    """
}

process TMM_NORMALISE_COUNTS {
    tag "$name"
    label 'process_medium'
    publishDir "${params.outdir}/gene_counts", mode: 'copy'

    input:
    tuple val(name), path(merged_counts)
    tuple val(name), path(merged_lens)
    path gene_subset

    output:
    path 'kallisto_scaled_counts.tsv', emit: tmm_counts

    script:
    """
    TMM_normalise_counts.R \
        -c $merged_counts \
        -l $merged_lens \
        -g $gene_subset \
        -r FALSE -p TRUE -t TRUE \
        -o kallisto_tmm_counts.tsv
    """
}
