##### Base image ###############################################################

FROM ralatsdio/biocontainers:v1.2.0

##### Metadata #################################################################

LABEL base_image="ralatsdio/biocontainers"
LABEL version="1.2.0"
LABEL software="REPdenovo"
LABEL software.version="v0.1.0"
LABEL about.summary="Designed for constructing repeats directly from sequence reads"
LABEL about.home="https://git@github.com:ralatsdc/REPdenovo.git"
LABEL about.provides="main.py v0.1.0"
LABEL about.license_file="MIT License"
LABEL about.license="https://github.com/ralatsdc/REPdenovo/blob/master/LICENSE"
LABEL about.tags="Genomics"
LABEL extra.identifiers.biotools="REPdenovo"
LABEL extra.binaries="main.py"

##### Installation #############################################################

USER root
ENV INSTALL_PFX=/home/biodocker

# Install Jellyfish 2.2.10
ENV JELLYFISH_VER=cf201901
RUN conda install -c bioconda/label/$JELLYFISH_VER jellyfish

# Build and install Velvet 1.2.10
ENV VELVET_VER=1.2.10
RUN cd $INSTALL_PFX && \
    wget https://github.com/BioDocker/software-archive/releases/download/velvet/velvet_$VELVET_VER.tgz && \
    tar -zxvf velvet_$VELVET_VER.tgz && \
    rm velvet_$VELVET_VER.tgz && \
    cd velvet_$VELVET_VER && \
    make MAXKMERLENGTH=60 && \
    cp velvet* $INSTALL_PFX/bin && \
    cd .. && \
    rm -rf velvet_$VELVET_VER

# Install BWA-MEM 0.7.17
ENV BWA_VER=cf201901
RUN conda install -c bioconda/label/$BWA_VER bwa

# Install Samtools 1.10
ENV SAMTOOLS_VER=cf201901
RUN conda install -c bioconda/label/$SAMTOOLS_VER samtools

# Install REPdenovo
ENV REPDENOVO_VER=0.1.0
RUN cd $INSTALL_PFX && \
    wget https://github.com/ralatsdc/REPdenovo/archive/v$REPDENOVO_VER.tar.gz -O REPdenovo-$REPDENOVO_VER.tar.gz && \
    tar -zxvf REPdenovo-$REPDENOVO_VER.tar.gz && \
    rm REPdenovo-$REPDENOVO_VER.tar.gz && \
    mv REPdenovo-$REPDENOVO_VER/* $INSTALL_PFX/bin && \
    rm -rf REPdenovo-$REPDENOVO_VER && \
    chown -R biodocker:biodocker $INSTALL_PFX/bin && \
    chmod +x $INSTALL_PFX/bin/TERefiner_1 && \
    chmod +x $INSTALL_PFX/bin/ContigsMerger && \
    chmod +x $INSTALL_PFX/bin/main.py

##### Configuration ############################################################

USER biodocker
WORKDIR /data
# CMD ['main.py']

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
