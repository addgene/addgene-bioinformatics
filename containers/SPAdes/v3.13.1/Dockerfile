##### Base image ###############################################################

FROM ralatsdio/biocontainers:v1.1.0

##### Metadata #################################################################

LABEL base_image="ralatsdio/biocontainers"
LABEL version="1.1.0"
LABEL software="SPAdes"
LABEL software.version="3.13.1"
LABEL about.summary="SPAdes – St. Petersburg genome assembler – is an assembly toolkit containing various assembly pipelines."
LABEL about.home="https://github.com/ablab/spades"
LABEL about.provides="spades.pl v3.13.1"
LABEL about.license_file="https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html"
LABEL about.license="GPLv2"
LABEL about.tags="Genomics"
LABEL extra.identifiers.biotools="SPAdes"
LABEL extra.binaries="spades.py"

##### Installation #############################################################

USER root
ENV INSTALL_PFX=/home/biodocker

# Install SPAdes
ENV ZIP=SPAdes-3.13.1-Linux.tar.gz
ENV URL=http://cab.spbu.ru/files/release3.13.1
ENV FOLDER=SPAdes-3.13.1-Linux
RUN wget $URL/$ZIP -O $INSTALL_PFX/$ZIP && \
    tar xvf $INSTALL_PFX/$ZIP -C $INSTALL_PFX && \
    rm $INSTALL_PFX/$ZIP && \
    cd $INSTALL_PFX/$FOLDER && \
    mv $INSTALL_PFX/$FOLDER/bin/* $INSTALL_PFX/bin && \
    mv $INSTALL_PFX/$FOLDER/share $INSTALL_PFX && \
    rm -rf $INSTALL_PFX/$FOLDER

##### Configuration ############################################################

USER biodocker
WORKDIR /data
# CMD ['NOVOPlasty.pl']

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
