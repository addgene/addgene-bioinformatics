##### Base image ###############################################################

FROM ralatsdio/biocontainers:v1.1.0

##### Metadata #################################################################

LABEL base.image="ralatsdio/biocontainers"
LABEL version="1.1.0"
LABEL software="SKESA"
LABEL about.summary="Strategic Kmer Extension for Scrupulous Assemblies"
LABEL about.home="https://github.com/ncbi/SKESA"
LABEL about.license_file="https://github.com/ncbi/SKESA/blob/master/LICENSE"
LABEL about.license="?"
LABEL about.tags="Genomics"
LABEL extra.identifiers.biotools="SKESA"
LABEL extra.binaries="skesa"

##### Installation #############################################################

USER root
ENV INSTALL_PFX=/home/biodocker

# Install SKESA
RUN cd $INSTALL_PFX && \
    git clone https://github.com/ncbi/SKESA && \
    cd SKESA && \
    apt-get update && apt-get -y -qq install git build-essential bison libxml2 libboost-all-dev libz-dev && \
    make 
    # cp SKESA/skesa $INSTALL_PFX/bin/ && \
    # chmod +x $INSTALL_PFX/bin/skesa

##### Configuration ############################################################

USER biodocker
WORKDIR /data
# CMD ['skesa']

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
