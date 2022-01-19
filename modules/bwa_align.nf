process MAKE_BWA_INDEX {
    tag "$name"
    label 'process_high'
    publishDir "${params.outdir}/bwa_idx", mode: 'copy'

    input:
    path ref_genome

    output:
    path '*', emit: bwa_idx

    script:
    """
    bwa index $ref_genome
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
    tuple val(name), path('*.bam'), path('*.bai'), emit: bwa_out

    script:
    """
    bwa mem -t ${task.cpus} ${idx} ${reads} | samtools sort -@ ${task.cpus - 1} -O bam - > ${name}.bam
    samtools index -@ ${task.cpus} ${name}.bam
    """
}

process MERGE_COUNTS {
    tag "$name"
    label 'process_high'
    publishDir "${params.outdir}/bwa_aln", mode: 'copy'

    input:
    path gpa_file
    tuple val(name), path(bwa_dir)
    path meta_merged
    path st_file

    output:
    tuple val(name), path('*.tsv'), emit: bwa_merged_counts

    script:
    """
    merge_bwa_counts.py \
        --gene_presence_absence=$gpa_file \
        --quant_dir=$bwa_dir \
        --metadata_merged=$meta_merged \
        --ST_file=$st_file \
        --outf=bwa_merged_counts.tsv
    """
}

process MERGE_LENS {
    tag "$name"
    label 'process_high'
    publishDir "${params.outdir}/bwa_aln", mode: 'copy'

    input:
    path gpa_file
    tuple val(name), path(bwa_dir)
    path meta_merged
    path st_file

    output:
    tuple val(name), path('*.tsv'), emit: bwa_merged_lens

    script:
    """
    merge_bwa_lens.py \
        --gene_presence_absence=$gpa_file \
        --quant_dir=$bwa_dir \
        --metadata_merged=$meta_merged \
        --ST_file=$st_file \
        --outf=bwa_merged_lens.tsv
    """
}
