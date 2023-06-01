FROM ubuntu:20.04

# Avoid timezone interactive configuration
ENV DEBIAN_FRONTEND=noninteractive

# Update the system and install necessary tools
RUN apt-get update -y && apt-get install -y \
    build-essential \
    wget \
    curl \
    python3 \
    python3-pip \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libnss3-dev \
    libblas-dev \
    liblapack-dev \
    gfortran \
    libreadline-dev \
    libx11-dev \
    libxt-dev \
    x11proto-core-dev \
    libcairo2-dev \
    libpng-dev \
    unzip \
    libxml2-dev \
    gawk

# Define environment variable with R version
ENV R_VERSION 4.3.0

# Download and install R from source
RUN wget https://cran.r-project.org/src/base/R-4/R-${R_VERSION}.tar.gz \
  && tar xzf R-${R_VERSION}.tar.gz \
  && cd R-${R_VERSION} \
  && ./configure --enable-R-shlib --with-blas --with-lapack \
  && make \
  && make install

# Installing Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2 && \
    tar -xvjf samtools-1.14.tar.bz2 && \
    cd samtools-1.14 && \
    ./configure && \
    make && \
    make install

#Installing Bowtie2
RUN apt-get update && apt-get install -y \
    bowtie2

# Installing Bioconductor and necessary packages
RUN R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('GenomicFeatures', 'TCGAbiolinks', 'SummarizedExperiment'))"

# Installing BEDTools
RUN apt-get update && apt-get install -y bedtools

WORKDIR /data
