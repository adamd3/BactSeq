process MAKE_BWA_INDEX {
    tag "$ref_genome"
    label 'process_high'
    publishDir "${params.outdir}/bwa_idx", mode: 'copy'

    input:
    path ref_genome

    output:
    path '*', emit: bwa_idx

    script:
    """
    bwa index $ref_genome -p ref_idx
    """
}

process BWA_ALIGN {
    tag "$name"
    label 'process_high'
    publishDir "${params.outdir}/bwa_aln", mode: 'copy'

    input:
    tuple val(name), path(reads)
    path idx

    output:
    path '*.bam', emit: bam_files
    path '*.bai', emit: bai_files
    path '*.counts', emit: count_files

    script:
    """
    bwa mem -t ${task.cpus} ref_idx ${reads} | samtools sort -@ ${task.cpus - 1} -O bam - > ${name}.bam
    samtools index -@ ${task.cpus} ${name}.bam
    samtools idxstats ${name}.bam | head -n 1 > ${name}.counts
    """
}
