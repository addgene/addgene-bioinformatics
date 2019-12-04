##### Base image ###############################################################

FROM ralatsdio/biocontainers:v1.1.0

##### Metadata #################################################################

LABEL base.image="ralatsdio/biocontainers"
LABEL version="1.1.0"
LABEL software="BBTools"
LABEL software.version="v38.72"
LABEL about.summary="BBMap, a short read aligner, as well as various other bioinformatic tools"
LABEL about.home="https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/"
LABEL about.provides="BBTools v38.72"
LABEL about.license_file="NA"
LABEL about.license="NA"
LABEL about.tags="Genomics"
LABEL extra.identifiers.biotools="BBTools"
LABEL extra.binaries="Various"

##### Installation #############################################################

USER root
ENV INSTALL_PFX=/home/biodocker

# Install dependencies
RUN apt-get update && \
    apt-get install -y \
        openjdk-9-jre && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install BBTools
RUN mkdir $INSTALL_PFX/bbtools && \
  cd $INSTALL_PFX/bbtools && \
  wget https://sourceforge.net/projects/bbmap/files/BBMap_38.73.tar.gz && \
  tar -zxvf BBMap_38.73.tar.gz \
  cp bbmap/* $INSTALL_PFX/bin

##### Configuration ############################################################

USER biodocker
WORKDIR /data
# CMD ['Various']

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
