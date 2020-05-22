##### Base image ###############################################################

FROM ralatsdio/biocontainers:v1.1.0

##### Metadata #################################################################

LABEL base_image="ralatsdio/biocontainers"
LABEL version="1.1.0"
LABEL software="SPAdes"
LABEL software.version="3.14.1"
LABEL about.summary="SPAdes – St. Petersburg genome assembler – is an assembly toolkit containing various assembly pipelines."
LABEL about.home="https://github.com/ablab/spades"
LABEL about.provides="spades.pl v3.14.1"
LABEL about.license_file="https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html"
LABEL about.license="GPLv2"
LABEL about.tags="Genomics"
LABEL extra.identifiers.biotools="SPAdes"
LABEL extra.binaries="spades.py"

##### Installation #############################################################

USER root
ENV INSTALL_PFX=/home/biodocker

# Install SPAdes
ENV SPADES_VER=3.14.1
RUN cd $INSTALL_PFX && \
    wget http://cab.spbu.ru/files/release$SPADES_VER/SPAdes-$SPADES_VER-Linux.tar.gz && \
    tar -zxvf SPAdes-$SPADES_VER-Linux.tar.gz && \
    rm SPAdes-$SPADES_VER-Linux.tar.gz && \
    cd SPAdes-$SPADES_VER-Linux && \
    mv bin/* $INSTALL_PFX/bin && \
    mv share $INSTALL_PFX && \
    cd .. && \
    rm -rf SPAdes-$SPADES_VER-Linux

##### Configuration ############################################################

USER biodocker
WORKDIR /data
# CMD ['spades.py']

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
