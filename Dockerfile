FROM continuumio/miniconda3

ENV PATH=/opt/conda/envs/bact_seq-1.0.0/bin:$PATH
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "source activate bact_seq-1.0.0" >> ~/.bashrc && \
    git clone https://github.com/adamd3/BactSeq.git && \
    /opt/conda/bin/conda env create -f BactSeq/environment.yml && \
    rm -rf BactSeq

ENTRYPOINT ["bash", "-c"]

CMD ["exec \"$@\""]
