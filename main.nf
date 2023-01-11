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

if (params.ref_ann) {
    ch_gff_file = file(params.ref_ann, checkIfExists: true)
} else { exit 1, 'Reference genome GFF file not specified!' }

// optional inputs
if (params.cont_tabl) {
    ch_cont_file = file(params.cont_tabl, checkIfExists: true)
}
if (params.func_file) {
    ch_func_file = file(params.func_file, checkIfExists: true)
}



// optional functional enrichment step
//ch_func_file = ( params.func_file
//            ? Channel.empty()
//            : file(params.func_file, checkIfExists: true) )



/*
================================================================================
    Modules
================================================================================
*/
include {MAKE_META_FILE} from './modules/metadata'
include {TRIMGALORE} from './modules/trim_reads'
include {MAKE_BWA_INDEX; BWA_ALIGN; COUNT_READS} from './modules/bwa_align'
include {MAKE_KALLISTO_IDX; KALLISTO_QUANT; MERGE_COUNTS} from './modules/kallisto'
include {TMM_NORMALISE_COUNTS} from './modules/normalisation'
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
     * Make metadata file linking samples with FastQ files
     */
    MAKE_META_FILE (
        ch_samples
    )
    ch_metadata = MAKE_META_FILE.out.sample_metadata


    /*
     *  Create channels for input files
     */
    // ch_metadata
    //     .splitCsv(header:true, sep:'\t')
    //     .map {
    //         row -> [ row.sample, [ file(row.path_to_file, checkIfExists: true) ] ]
    //     }
    //     .set { ch_raw_reads_trimgalore }

    // ch_metadata
    //     .splitCsv(header: true, sep:'\t')
    //     .map { row -> row.sample }
    //     .set { ch_sample_ids }

    ch_metadata
        .splitCsv(header: true, sep:'\t')
        .map { create_fastq_channel(it) }
        .set { ch_raw_reads_trimgalore }


    /*
     *  Trim reads
     */
    if (params.skip_trimming) {
        ch_trimmed_reads = ch_raw_reads_trimgalore.collect()
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

        BWA_ALIGN (
            ch_trimmed_reads,
            ch_bwa_idx
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
            ch_gff_file,
            ch_metadata
        )
        ch_readcounts_df = MERGE_COUNTS.out.counts_df
        ch_readcounts_df_pc = MERGE_COUNTS.out.counts_df_pc
        ch_refgene_df = MERGE_COUNTS.out.ref_gene_df


    } else { exit 1, 'aligner not valid: please choose one of `bwa` or `kallisto`' }

    /*
     *  Get normalised read counts per gene
     */
    TMM_NORMALISE_COUNTS (
        ch_readcounts_df_pc,
        ch_refgene_df
    )
    ch_cpm_counts = TMM_NORMALISE_COUNTS.out.cpm_counts
    ch_rpkm_counts = TMM_NORMALISE_COUNTS.out.rpkm_counts
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
    if (params.cont_tabl) {
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
    if (params.func_file) {
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
    The typical command for running the pipeline is as follows:
      nextflow run BactSeq --data_dir [dir] --sample_file [file] --ref_genome [file] --ref_ann [file] -profile docker

    Mandatory arguments:
      --data_dir [file]               Path to directory containing FastQ files.
      --sample_file [file]            Path to file containing sample information.
      --ref_genome [file]             Path to FASTA file containing reference genome sequence (bwa) or multi-FASTA file containing coding gene sequences (kallisto).
      --ref_ann [file]                Path to GFF file containing reference genome annotation.
      --aligner [str]                 (Pseudo-)aligner to be used. Options: `bwa`, `kallisto`. Default = bwa.
      --paired [str]                  Is data paired-end? Default = FALSE.
      --strandedness [str]            Is data stranded? Options: `unstranded`, `forward`, `reverse`. Default = reverse.
      -profile [str]                  Configuration profile to use.
                                      Available: conda, docker, singularity.

    Other options:
      --fragment_len [str]            Estimated average fragment length for kallisto transcript quantification (only required for single-end reads). Default = 150.
      --fragment_sd [str]             Estimated standard deviation of fragment length for kallisto transcript quantification (only required for single-end reads). Default = 20.
      --cont_tabl [file]              Path to tsv file containing contrasts to be performed for differential expression.
      --func_file [file]              Path to GMT-format file containing functional annotation.
      --p_thresh [str]                Adjusted p-value threshold for identifying differentially expressed genes. Default = 0.05.
      --l2fc_thresh [str]             Absolute log2(FoldChange) threshold for identifying differentially expressed genes. Default = 1.
      --skip_trimming [bool]          Do not trim adaptors from FastQ files.
      --outdir [file]                 The output directory where the results will be saved (Default: './results').
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    """.stripIndent()
}
