process MAKE_CLONE_FASTA {
     tag "$name"
     maxForks 20 // maximum number of files to process in parallel (TODO: make this a parameter)
     publishDir "${params.outdir}/clone_fasta", mode: 'copy'

     input:
     path gpa 
     tuple val(name), path(genes)

     output:
     tuple val(name), path('*.fna'), emit: clone_fasta

     script:
     """
     make_single_clone_fasta.py $genes $gpa $name
     """
 }
