##### Base image ###############################################################

FROM ralatsdio/biocontainers:v1.1.0

##### Metadata #################################################################

LABEL base.image="ralatsdio/biocontainers"
LABEL version="1.1.0"
LABEL software="MaSuRCA"
LABEL software.version="v3.3.1"
LABEL about.summary="The Maryland Super Read Cabog Assembler (MaSuRCA)"
LABEL about.home="https://github.com/alekseyzimin/masurca"
LABEL about.provides="masurca v3.3.1"
LABEL about.license_file="https://github.com/alekseyzimin/masurca/blob/master/LICENSE.txt"
LABEL about.license="GNU General Public License v3.0"
LABEL about.tags="Genomics"
LABEL extra.identifiers.biotools="MaSuRCA"
LABEL extra.binaries="masurca"

##### Installation #############################################################

USER root
ENV INSTALL_PFX=/home/biodocker

# Install dependencies
RUN apt-get update && \
    apt-get install -y \
        libboost-all-dev \
	libbz2-dev && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install MaSuRCA
ENV MASURCA_VERSION=3.3.1
WORKDIR /software
RUN cd $INSTALL_PFX && \
    wget https://github.com/alekseyzimin/masurca/releases/download/v$MASURCA_VERSION/MaSuRCA-$MASURCA_VERSION.tar.gz
    tar -xf *.tar.gz && \
    rm *.tar.gz && \
    mv * masurca && \
    cd masurca && \
    chmod +x install.sh && \
    ./install.sh && \
    cp bin/masurca $INSTALL_PFX/bin && \
    chmod +x $INSTALL_PFX/bin

##### Configuration ############################################################

USER biodocker
WORKDIR /data
# CMD ['masurca']

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
