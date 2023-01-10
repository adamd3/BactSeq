
process MAKE_KALLISTO_IDX {
    tag "$ref_genome"
    label 'process_high'
    publishDir "${params.outdir}/kallisto_idx", mode: 'copy'

    input:
    path ref_genome

    output:
    path '*.kidx', emit: kallisto_idx

    """
    kallisto index -i ref_genome.kidx $ref_genome
    """
}

process KALLISTO_QUANT {
    tag "$meta.sample_id"
    label 'process_high'
    publishDir "${params.outdir}/kallisto_quant", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path idx
    val strandedness

    output:
    path "kallisto_${meta.sample_id}", emit: kallisto_out_dirs

    script:
    def name = task.ext.prefix ?: "${meta.sample_id}"

    if (strandedness == 'unstranded'){
        strand_arg = ""
    } else if (strandedness == 'forward'){
        strand_arg = "--fr-stranded"
    } else if (strandedness == 'reverse'){
        strand_arg = "--rf-stranded"
    }

    if (meta.paired_end) {
        // if trimming has not been performed, must symlink to match expected file names
        // kallisto params -l, -s are estimated from paired end data, but are required when using --single
        """
        [ ! -f  ${name}_1_val_1.fq.gz ] && ln -s ${reads[0]} ${name}_1_val_1.fq.gz
        [ ! -f  ${name}_2_val_2.fq.gz ] && ln -s ${reads[1]} ${name}_2_val_2.fq.gz
        kallisto quant -t $task.cpus -i ref_genome.kidx \
            ${strand_arg} -o kallisto_${name} ${name}_1_val_1.fq.gz ${name}_2_val_2.fq.gz
        """
    } else {
        """
        [ ! -f  ${name}_trimmed.fq.gz ] && ln -s $reads ${name}_trimmed.fq.gz
        kallisto quant -t $task.cpus -i ref_genome.kidx \
            ${strand_arg} --single -l ${params.fragment_len} -s ${params.fragment_sd} \
            -o kallisto_${name} ${name}_trimmed.fq.gz
        """
    }
}

process MERGE_COUNTS {
    tag "merge_counts"
    label 'process_high'
    publishDir "${params.outdir}/kallisto_quant", mode: 'copy'

    input:
    path kallisto_dirs
    path meta

    output:
    path 'kallisto_merged_counts.tsv', emit: kallisto_merged_out
    path 'ref_gene_df.tsv', emit: ref_gene_df

    script:
    """
    merge_kallisto_counts.py \
        --metadata_f=$meta \
        --out_dir="./"
    """
}

