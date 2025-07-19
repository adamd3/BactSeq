process MAKE_META_FILE {
    tag "$sample_file"
    label 'process_low'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path sample_file

    output:
    path 'sample_metadata.tsv', emit: sample_metadata

    script:
    """
    echo "Starting make_meta_file.py"
    echo "Current directory: \$(pwd)"
    echo "Files before:"
    ls -la
    echo "Running: make_meta_file.py $sample_file ${params.data_dir} sample_metadata.tsv"
    make_meta_file.py $sample_file ${params.data_dir} sample_metadata.tsv
    echo "Exit code: \$?"
    echo "Files after:"
    ls -la
    echo "Checking for target file:"
    ls -la sample_metadata.tsv || echo "sample_metadata.tsv not found"
    """
}
