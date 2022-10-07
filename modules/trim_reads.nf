process TRIMGALORE {
    tag "$meta.sample_id"
    label 'process_high'
    publishDir "${params.outdir}/trim_galore", mode: 'copy',
        saveAs: { filename ->
                      if (filename.endsWith('.html')) "fastqc/$filename"
                      else if (filename.endsWith('.zip')) "fastqc/zips/$filename"
                      else if (filename.endsWith('trimming_report.txt')) "logs/$filename"
                      else params.save_trimmed ? filename : null
                }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*{trimmed,val}*.fq.gz"), emit: trimmed_reads
    tuple val(meta), path("*.txt")                , emit: trimgalore_results_mqc
    tuple val(meta), path("*.{zip,html}")         , emit: trimgalore_fastqc_reports_mqc
    tuple val(meta), path("*unpaired*.fq.gz")     , emit: unpaired, optional: true

    script:
    // Calculate number of --cores for TrimGalore based on value of task.cpus
    // See: https://github.com/FelixKrueger/TrimGalore/blob/master/Changelog.md#version-060-release-on-1-mar-2019
    // See: https://github.com/nf-core/atacseq/pull/65
    // (Max no cores = 4, since there are diminishing returns beyond this)
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 3
        if (meta.paired_end) cores = (task.cpus as int) - 4
        if (cores < 1) cores = 1
        if (cores > 4) cores = 4
    }

    def name = task.ext.prefix ?: "${meta.sample_id}"

    if (meta.paired_end) {
        """
        [ ! -f  ${name}_1.fastq.gz ] && ln -s ${reads[0]} ${name}_1.fastq.gz
        [ ! -f  ${name}_2.fastq.gz ] && ln -s ${reads[1]} ${name}_2.fastq.gz
        trim_galore --cores $cores --fastqc --paired --gzip \\
            ${name}_1.fastq.gz ${name}_2.fastq.gz
        """
    } else {
        """
        [ ! -f  ${name}.fastq.gz ] && ln -s $reads ${name}.fastq.gz
        trim_galore --cores $cores --fastqc --gzip ${name}.fastq.gz
        """
    }
}
