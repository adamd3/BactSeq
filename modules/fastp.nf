process FASTP {
    tag "$meta.sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/fastp", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.html') || filename.endsWith('.json')) filename
            else params.save_trimmed ? filename : null
        }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastp.fq.gz"), emit: reads
    tuple val(meta), path("*.json")       , emit: json
    tuple val(meta), path("*.html")       , emit: html

    script:
    def name = task.ext.prefix ?: "${meta.sample_id}"
    def cores = task.cpus ? Math.min(task.cpus as int, 16) : 2
    def extra_args = params.fastp_args ?: ""

    if (meta.paired_end) {
        """
        fastp \\
            --in1 ${reads[0]} \\
            --in2 ${reads[1]} \\
            --out1 ${name}_1.fastp.fq.gz \\
            --out2 ${name}_2.fastp.fq.gz \\
            --json ${name}.fastp.json \\
            --html ${name}.fastp.html \\
            --thread $cores \\
            --length_required ${params.min_read_len} \\
            --detect_adapter_for_pe \\
            ${extra_args}
        """
    } else {
        """
        fastp \\
            --in1 ${reads} \\
            --out1 ${name}.fastp.fq.gz \\
            --json ${name}.fastp.json \\
            --html ${name}.fastp.html \\
            --thread $cores \\
            --length_required ${params.min_read_len} \\
            ${extra_args}
        """
    }
}
