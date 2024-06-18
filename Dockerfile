FROM continuumio/miniconda3

RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda \
    && conda config --add channels defaults

RUN conda install -y \
    numpy>=1.21 \
    python>=3.7,<3.9 \
    gawk>=5.1 \
    pigz>=2.3 \
    pandas>=1.4 \
    r-base>=3.6 \
    r-optparse>=1.7 \
    r-rcolorbrewer>=1 \
    r-reshape2>=1.4 \
    r-pheatmap>=1.0 \
    r-matrixstats>=0.6 \
    r-ape>=5.0 \
    r-xtable>=1.8 \
    r-rsqlite>=2.2 \
    r-plyr>=1.8 \
    r-fastmap>=1.1 \
    r-devtools>=2.4 \
    r-scales>=1.1 \
    r-tidyverse>=1.9 \
    fastqc>=0.11 \
    trim-galore>=0.6 \
    bwa>=0.7 \
    samtools>=1.15 \
    bioconductor-edger>=3.36 \
    bioconductor-deseq2>=1.34 \
    bioconductor-topgo>=2.00 \
    bioconductor-go.db>=3.00 \
    bioconductor-rsubread>=2.00 \
    bioconductor-enhancedvolcano>=1.00

WORKDIR /data/
