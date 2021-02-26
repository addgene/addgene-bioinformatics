##### Base image ###############################################################

FROM ralatsdio/biocontainers:v1.2.0

##### Metadata #################################################################

LABEL base.image="ralatsdio/biocontainers"
LABEL version="1.2.0"
LABEL software="SKESA"
LABEL software.version="2.4.0"
LABEL about.summary="Strategic Kmer Extension for Scrupulous Assemblies"
LABEL about.home="https://github.com/ncbi/SKESA"
LABEL about.provides="SKESA v2.4.0"
LABEL about.license_file="https://github.com/ncbi/SKESA/blob/master/LICENSE"
LABEL about.license="?"
LABEL about.tags="Genomics"
LABEL extra.identifiers.biotools="SKESA"
LABEL extra.binaries="skesa"

##### Installation #############################################################

USER root
ENV INSTALL_PFX=/home/biodocker

# Install SKESA
ENV SKESA_VER=2.4.0
ENV CENTOS_VER=7.7
RUN cd $INSTALL_PFX && \
    mkdir skesa && \
    cd skesa && \  
    wget https://github.com/ncbi/SKESA/releases/download/$SKESA_VER/skesa.centos.$CENTOS_VER && \
    mv skesa.centos.$CENTOS_VER $INSTALL_PFX/bin/skesa && \
    chmod +x $INSTALL_PFX/bin/skesa

##### Configuration ############################################################

USER biodocker
WORKDIR /data
# CMD ['skesa']

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
