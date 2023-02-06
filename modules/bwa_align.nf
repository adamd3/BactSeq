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
    tag "$meta.sample_id"
    label 'process_high'
    publishDir "${params.outdir}/bwa_aln", mode: 'copy'

    input:
    tuple val(meta), path(trimmed_reads)
    path idx

    output:
    path '*.bam', emit: bam_files
    path '*.bai', emit: bai_files
    path '*.counts', emit: count_files

    script:

    def name = task.ext.prefix ?: "${meta.sample_id}"

    if (meta.paired_end) {
        """
        bwa mem -t ${task.cpus} ref_idx \\
            ${name}_1_val_1.fq.gz ${name}_2_val_2.fq.gz | \\
            samtools sort -@ ${task.cpus - 1} -O bam - > ${name}.bam
        samtools index -@ ${task.cpus} ${name}.bam
        # samtools idxstats ${name}.bam | head -n 1 > ${name}.counts
        samtools view -F 0x4 ${name}.bam | cut -f 1 | sort | uniq | wc -l > ${name}.counts
        """
    } else {
        """
        bwa mem -t ${task.cpus} ref_idx ${name}_trimmed.fq.gz \\
            | samtools sort -@ ${task.cpus - 1} -O bam - > ${name}.bam
        samtools index -@ ${task.cpus} ${name}.bam
        # samtools idxstats ${name}.bam | head -n 1 > ${name}.counts
        samtools view -F 0x4 ${name}.bam | cut -f 1 | sort | uniq | wc -l > ${name}.counts
        """
    }
}

process COUNT_READS {
    tag "$gff"
    label 'process_high'
    publishDir "${params.outdir}/read_counts", mode: 'copy'

    input:
    path bam
    path bai
    path counts
    path meta
    path gff
    val paired
    val strandedness

    output:
    path 'gene_counts.tsv', emit: counts_df
    path 'gene_counts_pc.tsv', emit: counts_df_pc
    path 'counts_summary.tsv', emit: counts_summary
    path 'ref_gene_df.tsv', emit: ref_gene_df
    path 'library_composition.png', emit: libcomp_plot
    path 'library_composition_proportions.png', emit: libcomp_plot_prop

    script:

    if (paired) {
        """
        count_reads.R -p TRUE -s $strandedness -m $meta -g $gff -t ${task.cpus}
        """
    } else {
        """
        count_reads.R -p FALSE -s $strandedness -m $meta -g $gff -t ${task.cpus}
        """
    }
}
