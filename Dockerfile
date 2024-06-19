FROM continuumio/miniconda3

ENV PATH=/opt/conda/envs/bact_seq-1.0.0/bin:$PATH
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN useradd -m -s /bin/bash nextflow_user
USER nextflow_user
WORKDIR /home/nextflow_user

COPY environment.yml .

RUN conda env create -f environment.yml && conda clean -a

SHELL ["conda", "run", "-n", "bact_seq-1.0.0", "/bin/bash", "-c"]

ENTRYPOINT ["conda", "run", "-n", "bact_seq-1.0.0", "/bin/bash", "-c"]

CMD ["bash"]
