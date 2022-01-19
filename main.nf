#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
================================================================================
    BactSeq: bacterial RNA-Seq data analysis
================================================================================
    Github : [github.com/adamd3/BactSeq]
*/

// /home/adam/nextflow run /home/adam/strain_seq \
//     --data_dir /projects/pseudomonas_transcriptomics/storage/fastq_files \
//     --sample_file /projects/pseudomonas_transcriptomics/storage/hzi_meta.txt \
//     -profile docker -resume


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

// ## to get the sample sheet from an ENA file list (for fetchngs):
// awk -F '\t' 'BEGIN{print "sample\tfilename\tgroup\trepeat"} FNR > 1 {printf "%s\t%s\t%s\t%s\n", $3,$3"_T1.fastq.gz",$10,"1"}' ENA_file_list.tsv > sample_info.tsv

if (params.sample_file) {
    ch_samples = file(params.sample_file, checkIfExists: true)
} else { exit 1, 'Sample file not specified!' }

// if (params.multifasta_file) {
//     ch_multifasta_file = file(params.multifasta_file, checkIfExists: true)
// } else { exit 1, 'Multi-fasta file not specified!' }

// if (params.faidx_file) {
//     ch_faidx_file = file(params.faidx_file, checkIfExists: true)
// } else { exit 1, 'Index for multi-fasta file not specified!' }

// if (params.gpa_file) {
//     ch_gpa_file = file(params.gpa_file, checkIfExists: true)
// } else { exit 1, 'Gene presence/absence file not specified!' }




/*
================================================================================
    Modules
================================================================================
*/
include {MAKE_META_FILE} from './modules/metadata'
include {TRIMGALORE} from './modules/trim_reads'
// include {MAKE_KALLISTO_INDEX; KALLISTO_QUANT; MERGE_COUNTS; MERGE_LENS} from './modules/kallisto'
// include {SUBSET_GENES; LENGTH_SCALE_COUNTS; TMM_NORMALISE_COUNTS} from './modules/normalisation'




/*
================================================================================
    Main workflow
================================================================================
*/
workflow {

    // /*
    //  * Make metadata file linking samples with FastQ files
    //  */
    MAKE_META_FILE (
        ch_samples
    )
    ch_metadata = MAKE_META_FILE.out.sample_metadata


    /*
     *  Create channels for input files
     */
    ch_metadata
        .splitCsv(header:true, sep:'\t')
        .map { row -> [ row.sample_id, [ file(row.path_to_file, checkIfExists: true) ] ] }
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


    // /*
    //  *  Create a Kallisto index for each strain
    //  */
    // MAKE_KALLISTO_INDEX (
    //     ch_clone_fasta
    // )
    // ch_kallisto_idx = MAKE_KALLISTO_INDEX.out.kallisto_idx
    //
    //
    // /*
    //  *  Quantify gene expression using Kallisto
    //  */
    // KALLISTO_QUANT (
    //     ch_trimmed_reads,
    //     ch_kallisto_idx
    // )
    // // NOTE: the output is a _directory_ containing the kallisto results
    // ch_kallisto_out = KALLISTO_QUANT.out.kallisto_out
    //
    // /*
    //  *  Merge counts
    //  */
    // MERGE_COUNTS (
    //     ch_gpa_file,
    //     ch_kallisto_out,
    //     ch_metadata,
    //     params.st_file
    // )
    // ch_kallisto_counts = MERGE_COUNTS.out.kallisto_merged_counts
    //
    // /*
    //  *  Merge effective gene lengths
    //  */
    // MERGE_LENS (
    //     ch_gpa_file,
    //     ch_kallisto_out,
    //     ch_metadata,
    //     params.st_file
    // )
    // ch_kallisto_lens = MERGE_LENS.out.kallisto_merged_lens
    //
    // /*
    //  *  Scale counts to median gene length across strains
    //  */
    // LENGTH_SCALE_COUNTS (
    //     ch_kallisto_counts,
    //     ch_kallisto_lens,
    //     ch_gene_subset
    // )
    // ch_scaled_counts = LENGTH_SCALE_COUNTS.out.scaled_counts
    //
    // /*
    //  *  Get size-factor-scaled, TMM-normalised counts
    //  */
    // TMM_NORMALISE_COUNTS (
    //     ch_kallisto_counts,
    //     ch_kallisto_lens,
    //     ch_gene_subset
    // )
    // ch_tmm_counts = TMM_NORMALISE_COUNTS.out.tmm_counts
    // // NB the resulting counts are log-transformed by default

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
      nextflow run strain_seq --data_dir [dir] --sample_file [file] --gpa_file [gene_presence_absence.csv] -profile docker

    Mandatory arguments:
      --data_dir [file]               Path to directory containing FastQ files.
      --sample_file [file]            Path to file containing sample information.
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated).
                                      Available: conda, docker

    Other options:
      --outdir [file]                 The output directory where the results will be saved (Default: './results').
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    """.stripIndent()
}
