FROM continuumio/miniconda3
LABEL maintainer="Adam Dinan <ad866@cam.ac.uk>"

WORKDIR /app

COPY environment.yml .

RUN conda env create -f environment.yml

RUN echo "source activate $(grep -m 1 -E 'name: *' environment.yml | cut -d' ' -f2)" \
    > ~/.bashrc

SHELL ["conda", "run", "-n", \
    "$(grep -m 1 -E 'name: *' environment.yml | cut -d' ' -f2)", "/bin/bash", "-c"]

# Open a shell by default when the container is run interactively
ENTRYPOINT ["bash"]
