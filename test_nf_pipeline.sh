/home/adam/nextflow run /home/adam/strain_seq \
    --data_dir /projects/pseudomonas_transcriptomics/storage/fastq_files \
    --meta_file /projects/pseudomonas_transcriptomics/storage/hzi_meta.txt \
    --gpa_file /projects/pseudomonas_transcriptomics/storage/gene_presence_absence.csv \
    --st_file /home/adam/pseudomonas_transcriptomics/CF_STs.txt \
    --perc 99 \
    -profile docker
