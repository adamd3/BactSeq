process TRIMGALORE {
    tag "$name"
    label 'process_high'
    publishDir "${params.outdir}/trim_galore", mode: 'copy',
        saveAs: { filename ->
                      if (filename.endsWith('.html')) "fastqc/$filename"
                      else if (filename.endsWith('.zip')) "fastqc/zips/$filename"
                      else if (filename.endsWith('trimming_report.txt')) "logs/$filename"
                      else params.save_trimmed ? filename : null
                }

    input:
    tuple val(name), path(reads)

    output:
    tuple val(name), path('*.fq.gz'), emit: trimmed_reads
    path '*.txt', emit: trimgalore_results_mqc
    path '*.{zip,html}', emit: trimgalore_fastqc_reports_mqc

    script:
    // Calculate number of --cores for TrimGalore based on value of task.cpus
    // See: https://github.com/FelixKrueger/TrimGalore/blob/master/Changelog.md#version-060-release-on-1-mar-2019
    // See: https://github.com/nf-core/atacseq/pull/65
    // (Max no cores = 4, since there are diminishing returns beyond this)
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 3
        if (cores < 1) cores = 1
        if (cores > 4) cores = 4
    }

    // Add symlinks to original fastqs for consistent naming in MultiQC
    """
    [ ! -f  ${name}.fastq.gz ] && ln -s $reads ${name}.fastq.gz
    trim_galore --cores $cores --fastqc --gzip ${name}.fastq.gz
    """
}
