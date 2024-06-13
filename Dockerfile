FROM continuumio/miniconda3
LABEL maintainer="Adam Dinan <ad866@cam.ac.uk>"

WORKDIR /app

COPY environment.yml .

RUN conda env create -f environment.yml

RUN ENV_NAME=$(grep -m 1 -E '^name: *' environment.yml | cut -d' ' -f2) && \
    echo "conda activate $ENV_NAME" > ~/.bashrc

SHELL ["/bin/bash", "-c"]

ENTRYPOINT ["bash", "-l", "-c", "source ~/.bashrc && bash"]
