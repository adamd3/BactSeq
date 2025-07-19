#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
================================================================================
    BactSeq: bacterial RNA-Seq data analysis
================================================================================
    Github : [github.com/adamd3/BactSeq]
*/


// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}


/*
================================================================================
    Validate inputs and create channels for files
================================================================================
*/

// required inputs
if (params.sample_file) {
    ch_samples = file(params.sample_file, checkIfExists: true)
} else { exit 1, 'Sample file not specified!' }

if (params.ref_genome) {
    ch_fasta_file = file(params.ref_genome, checkIfExists: true)
} else { exit 1, 'Reference genome FASTA file not specified!' }


// optional inputs
if (params.ref_ann) {
    ch_gff_file = file(params.ref_ann, checkIfExists: true)
} else {
    ch_gff_file = Channel.empty()
}

if (params.contrast_file) {
    ch_cont_file = file(params.contrast_file, checkIfExists: true)
} else {
    ch_cont_file = Channel.empty()
}
if (params.func_file) {
    ch_func_file = file(params.func_file, checkIfExists: true)
} else {
    ch_func_file = Channel.empty()
}






/*
================================================================================
    Modules
================================================================================
*/
include {MAKE_META_FILE} from './modules/metadata'
include {TRIMGALORE} from './modules/trim_reads'
include {MAKE_BWA_INDEX; BWA_ALIGN; COUNT_READS} from './modules/bwa_align'
include {MAKE_KALLISTO_IDX; KALLISTO_QUANT; MERGE_COUNTS} from './modules/kallisto'
include {NORMALISE_COUNTS} from './modules/normalisation'
include {PCA_SAMPLES} from './modules/plots'
include {DIFF_EXPRESSION} from './modules/diffexpr'
include {FUNC_ENRICHMENT} from './modules/func_enrich'



/*
================================================================================
    Functions
================================================================================
*/
// Function to get list of [ meta, [ fastq_1, file2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create sample metadata
    def meta = [:]
    meta.sample_id    = row.sample
    meta.paired_end   = row.paired.toBoolean()
   
    // add path(s) of the fastq file(s) to the metadata
    def fastq_meta = []
    if (!file(row.file1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.file1}"
    }
    if (meta.paired_end) {
        if (!file(row.file2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.file2}"
        }
        fastq_meta = [ meta, [ file(row.file1), file(row.file2) ] ]
    } else {
        fastq_meta = [ meta, [ file(row.file1) ] ]
    }
    return fastq_meta
}



/*
================================================================================
    Main workflow
================================================================================
*/
workflow {

    /*
     * Parameter validation
     */
    if (params.p_thresh <= 0 || params.p_thresh > 1) {
        exit 1, "ERROR: p_thresh must be between 0 and 1, got: ${params.p_thresh}"
    }
    if (params.l2fc_thresh < 0) {
        exit 1, "ERROR: l2fc_thresh must be >= 0, got: ${params.l2fc_thresh}"
    }
    if (!(params.aligner in ['bwa', 'kallisto'])) {
        exit 1, "ERROR: aligner must be 'bwa' or 'kallisto', got: ${params.aligner}"
    }
    if (!(params.strandedness in ['unstranded', 'forward', 'reverse'])) {
        exit 1, "ERROR: strandedness must be 'unstranded', 'forward', or 'reverse', got: ${params.strandedness}"
    }
    if (params.aligner == 'bwa' && !params.ref_ann) {
        exit 1, "ERROR: BWA aligner requires --ref_ann parameter"
    }

    /*
     * Make metadata file linking samples with FastQ files
     */
    MAKE_META_FILE (
        ch_samples
    )
    ch_metadata = MAKE_META_FILE.out.sample_metadata


    /*
     *  Create channels for input files
     */

    ch_metadata
        .splitCsv(header: true, sep:'\t')
        .map { create_fastq_channel(it) }
        .set { ch_raw_reads_trimgalore }


    /*
     *  Trim reads
     */
    if (params.skip_trimming) {
        ch_trimmed_reads = ch_raw_reads_trimgalore
        ch_trimgalore_results_mqc = Channel.empty()
        ch_trimgalore_fastqc_reports_mqc = Channel.empty()
    } else {
        TRIMGALORE (
            ch_raw_reads_trimgalore
        )
        ch_trimmed_reads = TRIMGALORE.out.trimmed_reads
        ch_trimgalore_results_mqc = TRIMGALORE.out.trimgalore_results_mqc
        ch_trimgalore_fastqc_reports_mqc = TRIMGALORE.out.trimgalore_fastqc_reports_mqc
    }



    /*
     *  Align / pseudo-align reads
     */
    if (params.aligner == "bwa") {

        MAKE_BWA_INDEX (
            ch_fasta_file
        )
        ch_bwa_idx = MAKE_BWA_INDEX.out.bwa_idx
        ch_ref_fasta = MAKE_BWA_INDEX.out.ref_fasta

        BWA_ALIGN (
            ch_trimmed_reads,
            ch_bwa_idx,
            ch_ref_fasta
        )
        ch_bwa_out_bam = BWA_ALIGN.out.bam_files.collect()
        ch_bwa_out_bai = BWA_ALIGN.out.bai_files.collect()
        ch_bwa_out_count = BWA_ALIGN.out.count_files.collect()

        /*
         *  Count reads mapped per gene; summarise library composition
         */
        COUNT_READS (
            ch_bwa_out_bam,
            ch_bwa_out_bai,
            ch_bwa_out_count,
            ch_metadata,
            ch_gff_file,
            params.paired,
            params.strandedness
        )
        ch_readcounts_df = COUNT_READS.out.counts_df
        ch_readcounts_df_pc = COUNT_READS.out.counts_df_pc
        ch_refgene_df = COUNT_READS.out.ref_gene_df

 
    } else if (params.aligner == "kallisto") {

        MAKE_KALLISTO_IDX (
            ch_fasta_file
        )
        ch_kallisto_idx = MAKE_KALLISTO_IDX.out.kallisto_idx

        KALLISTO_QUANT (
            ch_trimmed_reads,
            ch_kallisto_idx,
            params.strandedness,
        )
        ch_kallisto_out_dirs = KALLISTO_QUANT.out.kallisto_out_dirs.collect()

        /*
         *  Merge counts
         */
        MERGE_COUNTS (
            ch_kallisto_out_dirs,
            // ch_gff_file,
            ch_metadata
        )
        ch_readcounts_df = MERGE_COUNTS.out.counts_df
        ch_readcounts_df_pc = MERGE_COUNTS.out.counts_df_pc
        ch_refgene_df = MERGE_COUNTS.out.ref_gene_df


    } else { exit 1, 'aligner not valid: please choose one of `bwa` or `kallisto`' }

    /*
     *  Get normalised read counts per gene
     */
    NORMALISE_COUNTS (
        ch_readcounts_df_pc,
        ch_refgene_df
    )
    ch_deseq_counts = NORMALISE_COUNTS.out.deseq_counts
    ch_cpm_counts = NORMALISE_COUNTS.out.cpm_counts
    ch_rpkm_counts = NORMALISE_COUNTS.out.rpkm_counts
    // NB the resulting counts are log-transformed by default


    /*
     *  Principal component analysis (PCA) of samples
     */
    PCA_SAMPLES (
        ch_cpm_counts,
        ch_metadata
    )
    ch_pca_out = PCA_SAMPLES.out.pca_out

    /*
     *  Differential gene expression (DESeq2)
     */
    if (params.contrast_file) {
        DIFF_EXPRESSION (
            ch_readcounts_df_pc,
            ch_metadata,
            ch_cont_file,
            params.p_thresh,
            params.l2fc_thresh
        )
        ch_deseq_res = DIFF_EXPRESSION.out.deseq_res.collect()
    }

    /*
     *  Functional enrichment of DEGs (optional)
     */
    if (params.func_file && params.contrast_file) {
        FUNC_ENRICHMENT (
            ch_func_file,
            ch_deseq_res,
            params.p_thresh,
            params.l2fc_thresh
        )
        ch_func_enrich = FUNC_ENRICHMENT.out.func_res
    }

}

/*
================================================================================
    Completion summary
================================================================================
*/

c_green = "\033[0;32m";
c_reset = "\033[0m"

workflow.onComplete {
    log.info"""
    Execution status: ${ workflow.success ? 'OK' : 'failed' }
    ${c_green}Results are reported here: $params.outdir${c_reset}
    ￼""".stripIndent()
}


def helpMessage() {
    log.info"""
    Usage:
    nextflow run BactSeq --data_dir [dir] --sample_file [file] --ref_genome [file] --ref_ann [file] -profile conda [other_options]

    Mandatory arguments:
    --data_dir [file]               Path to directory containing FastQ files.
    --ref_genome [file]             Path to FASTA file containing reference genome sequence (bwa) or multi-FASTA file containing coding gene sequences (kallisto).
    --ref_ann [file]                Path to GFF file containing reference genome annotation.
    --sample_file [file]            Path to file containing sample information.
    -profile [str]                  Configuration profile to use.
                                    Available: conda, docker, singularity.

    Other options:
    --aligner [str]                 (Pseudo-)aligner to be used. Options: `bwa`, `kallisto`. Default = bwa.
    --contrast_file [file]          Path to tsv file containing contrasts to be performed for differential expression.
    --fragment_len [str]            Estimated average fragment length for kallisto transcript quantification (only required for single-end reads). Default = 150.
    --fragment_sd [str]             Estimated standard deviation of fragment length for kallisto transcript quantification (only required for single-end reads). Default = 20.
    --func_file [file]              Path to GFF3-format file containing functional annotations.
    --l2fc_thresh [str]             Absolute log2(FoldChange) threshold for identifying differentially expressed genes. Default = 1.
    --outdir [file]                 The output directory where the results will be saved (Default: './results').
    --paired [str]                  Data are paired-end.
    --p_thresh [str]                Adjusted p-value threshold for identifying differentially expressed genes. Default = 0.05.
    --skip_trimming [bool]          Do not trim adaptors from FastQ files.
    --strandedness [str]            Is data stranded? Options: `unstranded`, `forward`, `reverse`. Default = reverse.
    -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    -resume                         Re-start the pipeline if it has been previously run.

    """.stripIndent()
}
