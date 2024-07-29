FROM continuumio/miniconda

RUN apt-get  --allow-releaseinfo-change update -t oldoldstable && apt-get install -y \
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
    zlib1g-dev

RUN conda install -c daler \
    pip \
    cython \
    matplotlib \
    nose \
    numpydoc \
    pip \
    pandas \
    pyyaml \
    sphinx \
    pysam
RUN conda install -c daler \
    tabix \
    bedtools=2.25.0
ENV DISPLAY=:0
ENV LANG C.UTF-8
WORKDIR /opt/pybedtools

RUN pip install https://github.com/AndersenLab/bam-toolbox/archive/0.0.3.tar.gz
RUN apt-get install -y procps \
    && rm -rf /var/lib/apt/lists/*