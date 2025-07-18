FROM continuumio/miniconda3

WORKDIR /app

COPY environment.yml /app

RUN conda env create -f environment.yml

RUN echo "conda activate bact_seq-1.0.0" >> ~/.bashrc
ENV PATH=/opt/conda/envs/bact_seq-1.0.0/bin:$PATH

SHELL ["/bin/bash", "-c"]

COPY . /app

ENTRYPOINT ["bash", "-c", "source activate bact_seq-1.0.0 && exec bash"]
