##### Base image ###############################################################

FROM ralatsdio/biocontainers:v1.1.0

##### Metadata #################################################################

LABEL base.image="ralatsdio/biocontainers"
LABEL version="1.1.0"
LABEL software="MaSuRCA"
LABEL software.version="v4.0.1"
LABEL about.summary="The Maryland Super Read Cabog Assembler (MaSuRCA)"
LABEL about.home="https://github.com/alekseyzimin/masurca"
LABEL about.provides="masurca v4.0.1"
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
ENV MASURCA_VER=4.0.1
ENV DEST=$INSTALL_PFX
RUN cd $INSTALL_PFX && \
    wget https://github.com/alekseyzimin/masurca/releases/download/v$MASURCA_VER/MaSuRCA-$MASURCA_VER.tar.gz && \
    tar -zxvf MaSuRCA-$MASURCA_VER.tar.gz && \
    cd MaSuRCA-$MASURCA_VER && \
    ./install.sh && \
    cd $INSTALL_PFX/bin && \
    echo '#!/bin/bash' > masurca.sh && \
    echo 'masurca $1 && ./assemble.sh' >> masurca.sh && \
    chmod +x masurca.sh

##### Configuration ############################################################

USER biodocker
WORKDIR /data
# CMD ['masurca']

# Set path to executables
ENV PATH="$INSTALL_PFX/MaSuRCA-$MASURCA_VER/bin:${PATH}"

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
