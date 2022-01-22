
process COUNT_READS {
    tag "$bam"
    label 'process_high'
    publishDir "${params.outdir}/bwa_aln", mode: 'copy'

    input:
    path bam
    path bai
    path counts
    path meta
    path gff

    output:
    path "gene_counts.tsv", emit: counts_out

    script:
    """
    count_reads.R -m $meta -g $gff
    """
}
