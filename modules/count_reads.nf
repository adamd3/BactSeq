process COUNT_READS {
    tag "$bam"
    label 'process_high'
    publishDir "${params.outdir}/read_counts", mode: 'copy'

    input:
    path bam
    path bai
    path counts
    path meta
    path gff

    output:
    path '*.{tsv,png}', emit: counts_out

    script:
    """
    count_reads.R -m $meta -g $gff
    """
}
