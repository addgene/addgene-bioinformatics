##### Base image ###############################################################

FROM ralatsdio/biocontainers:v1.2.0

##### Metadata #################################################################

LABEL base.image="ralatsdio/biocontainers"
LABEL container.version="1.2.0"
LABEL software="Unicycler"
LABEL software.version="0.4.8"
LABEL about.summary="Unicycler is an assembly pipeline for bacterial genomes"
LABEL about.home="https://github.com/rrwick/Unicycler"
LABEL about.provides="Unicycler v0.4.8"
LABEL about.license_file="https://github.com/rrwick/Unicycler#license"
LABEL about.license="?"
LABEL about.tags="Genomics"
LABEL extra.identifiers.biotools="Unicycler"
LABEL extra.binaries="unicyler.py"

##### Installation #############################################################

USER root
ENV INSTALL_PFX=/home/biodocker

# Install dependencies
RUN apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository -y ppa:openjdk-r/ppa && \
    apt-get update && \
    apt-get install -y \
        bowtie2 \
        cmake \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
        libncursesw5-dev \
        libz-dev \
        ncbi-blast+ \
        openjdk-11-jre \
        python3 \
        python3-pip && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install Pilon
ENV PILON_VER=1.22
RUN cd $INSTALL_PFX && \
    wget https://github.com/broadinstitute/pilon/releases/download/v$PILON_VER/pilon-$PILON_VER.jar -P pilon && \
    echo '#!/bin/bash\njava -jar $INSTALL_PFX/pilon/pilon-$PILON_VER.jar $@' > pilon/pilon && \
    chmod +x pilon/pilon

# Install SPAdes
ENV SPADES_VER=3.13.0
RUN cd $INSTALL_PFX && \
    wget http://cab.spbu.ru/files/release$SPADES_VER/SPAdes-$SPADES_VER-Linux.tar.gz && \
    mkdir spades && \
    tar -zxvf SPAdes-$SPADES_VER-Linux.tar.gz -C spades --strip-components 1 && \
    rm SPAdes-$SPADES_VER-Linux.tar.gz

# Install Samtools and BCFtools
ENV SAMTOOLS_VER=1.9
ENV BCFTOOLS_VER=1.9
RUN cd $INSTALL_PFX && \
    wget https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VER/samtools-$SAMTOOLS_VER.tar.bz2 && \
    mkdir samtools && \
    tar -jxvf samtools-$SAMTOOLS_VER.tar.bz2 -C samtools --strip-components 1 && \
    rm samtools-$SAMTOOLS_VER.tar.bz2 && \
    cd samtools && \
    make && \
    make install && \
    cd $INSTALL_PFX && \
    wget https://github.com/samtools/bcftools/releases/download/$BCFTOOLS_VER/bcftools-$BCFTOOLS_VER.tar.bz2 && \
    mkdir bcftools && \
    tar -jxvf bcftools-$BCFTOOLS_VER.tar.bz2 -C bcftools --strip-components 1 && \
    rm bcftools-$BCFTOOLS_VER.tar.bz2 && \
    cd bcftools && \
    make && \
    make install

# Install Racon
ENV RACON_VER=1.3.1
RUN cd $INSTALL_PFX && \
    wget https://github.com/isovic/racon/releases/download/$RACON_VER/racon-v$RACON_VER.tar.gz && \
    mkdir racon && \
    tar -zxvf racon-v$RACON_VER.tar.gz -C racon --strip-components 1 && \
    rm racon-v$RACON_VER.tar.gz && \
    cd racon && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release ../ && \
    make && \
    make install

# Install Miniasm
ENV MINIASM_VER=0.3
RUN cd $INSTALL_PFX && \
    wget https://github.com/lh3/miniasm/archive/v$MINIASM_VER.tar.gz && \
    mkdir miniasm && \
    tar -zxvf v$MINIASM_VER.tar.gz -C miniasm --strip-components 1 && \
    rm v$MINIASM_VER.tar.gz && \
    cd miniasm && \
    make

# Install Unicycler
ENV UNICYCLER_VER=0.4.8
RUN cd $INSTALL_PFX && \
    wget https://github.com/rrwick/Unicycler/archive/v$UNICYCLER_VER.tar.gz && \
    mkdir unicycler && \
    tar -zxvf v$UNICYCLER_VER.tar.gz -C unicycler --strip-components 1 && \
    rm v$UNICYCLER_VER.tar.gz && \
    python3 $INSTALL_PFX/unicycler/setup.py install
    
##### Configuration ############################################################

USER biodocker
WORKDIR /data
# CMD ['unicylcer']

# Set path to executables
ENV PATH="$INSTALL_PFX/pilon:${PATH}"
ENV PATH="$INSTALL_PFX/spades/bin:${PATH}"
ENV PATH="$INSTALL_PFX/miniasm:${PATH}"

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
