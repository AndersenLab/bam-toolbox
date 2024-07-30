FROM continuumio/miniconda3

RUN apt-get update \
    && apt-get install -y \
        build-essential \
        bzip2 \
        ca-certificates \
        git \
        libglib2.0-0 \
        libsm6 \
        libxext6 \
        libxrender1 \
        mariadb-server \
        mariadb-client \
        wget \
        zlib1g-dev \
        procps \
    && rm -rf /var/lib/apt/lists/*

RUN conda install -c bioconda -c conda-forge -c anaconda \
    colorama \
    termcolor \
    docopt 

RUN conda install -c bioconda -c conda-forge -c anaconda \
    pandas \
    pybedtools

RUN conda install -c bioconda -c conda-forge -c anaconda \
    samtools
ENV DISPLAY=:0
ENV LANG=C.UTF-8
WORKDIR /opt/pybedtools

RUN pip install https://github.com/AndersenLab/bam-toolbox/archive/1.0.0.tar.gz