##### Base image ###############################################################

FROM ralatsdio/biocontainers:v1.2.0

##### Metadata #################################################################

LABEL base.image="ralatsdio/biocontainers"
LABEL version="1.2.0"
LABEL software="SSAKE"
LABEL software.version="v4.0.1"
LABEL about.summary="Genome assembly with short sequence reads"
LABEL about.home="https://github.com/ralatsdc/SSAKE.git"
LABEL about.provides="ssake v4.0.1"
LABEL about.license_file="https://github.com/ralatsdc/SSAKE/blob/master/LICENSE"
LABEL about.license="GNU General Public License v2.0"
LABEL about.tags="Genomics"
LABEL extra.identifiers.biotools="SSAKE"
LABEL extra.binaries=""

##### Installation #############################################################

USER root
ENV INSTALL_PFX=/home/biodocker

# Install dependencies
RUN apt-get update && \
    apt-get install -y \
	time && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install Perl modules
RUN cpan install Statistics::Descriptive

# Install SSAKE
ENV SSAKE_VER=v4.0.1
RUN cd $INSTALL_PFX && \
    git clone https://github.com/ralatsdc/SSAKE.git && \
    cd SSAKE && \
    git checkout $SSAKE_VER && \
    cp ssake/SSAKE $INSTALL_PFX/bin/ && \
    cp ssake/tools/* $INSTALL_PFX/bin/ && \
    cd .. && \
    rm -rf ssake

##### Configuration ############################################################

USER biodocker
WORKDIR /data
# CMD ['runSSAKE.sh']

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
