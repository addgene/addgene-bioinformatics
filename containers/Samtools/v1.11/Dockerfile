##### Base image ###############################################################

FROM ralatsdio/biocontainers:v1.2.0

##### Metadata #################################################################

LABEL base.image="ralatsdio/biocontainers"
LABEL container.version="1.2.0"
LABEL software="Samtools"
LABEL software.version="1.11"
LABEL about.summary="Tools (written in C using htslib) for manipulating next-generation sequencing data"
LABEL about.home="https://github.com/samtools/samtools"
LABEL about.provides="samtools v1.11"
LABEL about.license_file="https://github.com/samtools/samtools/blob/develop/LICENSE"
LABEL about.license="MIT/Expat"
LABEL about.tags="Genomics"
LABEL extra.identifiers.biotools="Samtools"
LABEL extra.binaries="samtools"

##### Installation #############################################################

USER root
ENV INSTALL_PFX=/home/biodocker

# Install dependencies
RUN apt-get update && \
    apt-get install -y \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install Samtools 
ENV SAMTOOLS_VER=1.11
RUN cd $INSTALL_PFX && \
    mkdir samtools && \
    cd samtools && \
    wget https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VER/samtools-$SAMTOOLS_VER.tar.bz2 && \
    tar -jxvf samtools-$SAMTOOLS_VER.tar.bz2 && \
    rm samtools-$SAMTOOLS_VER.tar.bz2 && \
    cd samtools-$SAMTOOLS_VER && \
    ./configure --prefix=$INSTALL_PFX && \
    make && \
    make install

##### Configuration ############################################################

USER biodocker
WORKDIR /data
# CMD ['samtools']

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
