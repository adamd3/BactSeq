process MAKE_KALLISTO_INDEX {
    tag "$name"
    label 'process_high'
    publishDir "${params.outdir}/kallisto_idx", mode: 'copy'

    input:
    tuple val(name), path(clone_fasta)

    output:
    tuple val(name), path('*.kidx'), emit: kallisto_idx

    script:
    """
    kallisto index -i ${name}.kidx ${clone_fasta}
    """
}

process KALLISTO_QUANT {
    tag "$name"
    label 'process_high'
    publishDir "${params.outdir}/kallisto_quant", mode: 'copy'

    input:
    tuple val(name), path(reads)
    tuple val(name2), path(idx)

    output:
    tuple val(name), path(name), emit: kallisto_out

    script:
    """
    kallisto quant -t $task.cpus --single -i $idx \
        --fr-stranded --single -l 150 -s 20 -o $name $reads
    """
}

process MERGE_COUNTS {
    tag "$name"
    label 'process_high'
    publishDir "${params.outdir}/kallisto_quant", mode: 'copy'

    input:
    path gpa_file
    tuple val(name), path(kallisto_dir)
    path meta_merged
    path st_file

    output:
    tuple val(name), path('*.tsv'), emit: kallisto_merged_counts

    script:
    """
    merge_kallisto_counts.py \
        --gene_presence_absence=$gpa_file \
        --quant_dir=$kallisto_dir \
        --metadata_merged=$meta_merged \
        --ST_file=$st_file \
        --outf=kallisto_merged_counts.tsv
    """
}

process MERGE_LENS {
    tag "$name"
    label 'process_high'
    publishDir "${params.outdir}/kallisto_quant", mode: 'copy'

    input:
    path gpa_file
    tuple val(name), path(kallisto_dir)
    path meta_merged
    path st_file

    output:
    tuple val(name), path('*.tsv'), emit: kallisto_merged_lens

    script:
    """
    merge_kallisto_lens.py \
        --gene_presence_absence=$gpa_file \
        --quant_dir=$kallisto_dir \
        --metadata_merged=$meta_merged \
        --ST_file=$st_file \
        --outf=kallisto_merged_lens.tsv
    """
}
