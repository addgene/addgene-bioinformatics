##### Base image ###############################################################

FROM ralatsdio/biocontainers:v1.2.0

##### Metadata #################################################################

LABEL base_image="ralatsdio/biocontainers"
LABEL version="1.2.0"
LABEL software="IDBA"
LABEL software.version="v1.1.4"
LABEL about.summary="the basic iterative de Bruijn graph assembler for second-generation sequencing reads"
LABEL about.home="https://github.com/ralatsdc/idba.git"
LABEL about.provides="idba v1.1.4"
LABEL about.license_file="https://github.com/ralatsdc/idba/blob/master/README.md"
LABEL about.license="GNU General Public License v2.0, or any later version"
LABEL about.tags="Genomics"
LABEL extra.identifiers.biotools="idba"
LABEL extra.binaries="idba"

##### Installation #############################################################

USER root
ENV INSTALL_PFX=/home/biodocker

# Install IDBA (Use ralatsdc fork)
ENV IDBA_VER=1.1.4
RUN cd $INSTALL_PFX && \
    git clone https://github.com/ralatsdc/idba.git && \
    cd idba && \
    git checkout $CSC_VER && \
    ./build.sh && \
    cd $INSTALL_PFX && \
    find $INSTALL_PFX/idba/bin -maxdepth 1 -type f -perm -111 -exec ln -s {} $INSTALL_PFX/bin \;

##### Configuration ############################################################

USER biodocker
WORKDIR /data
# CMD ['idba']

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
