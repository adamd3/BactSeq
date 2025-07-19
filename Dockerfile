FROM continuumio/miniconda3

WORKDIR /app

COPY environment.yml /app

RUN conda env create -f environment.yml

# Set conda environment as default
ENV CONDA_DEFAULT_ENV=bact_seq-1.0.0
ENV PATH=/opt/conda/envs/bact_seq-1.0.0/bin:$PATH

SHELL ["/bin/bash", "-c"]

COPY . /app
