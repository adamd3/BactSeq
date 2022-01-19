process MAKE_META_FILE {
    tag "$sample_file"
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path sample_file

    output:
    path 'sample_metadata.tsv', emit: sample_metadata

    script:
    """
    make_meta_file.py $sample_file ${params.data_dir} sample_metadata.tsv
    """
}
