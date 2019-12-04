##### Base image ###############################################################

FROM biocontainers/biocontainers:v1.0.0_cv4

##### Metadata #################################################################

LABEL base_image="biocontainers:v1.0.0_cv4"
LABEL version="3"
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

# Install BBTools
ENV BBTOOLS_URL=https://sourceforge.net/projects/bbmap/files/BBMap_38.72.tar.gz/download
ENV NOVOPLASTY_SRC=NOVOPlasty
ENV NOVOPLASTY_TAG=NOVOPlasty3.7
RUN git clone $NOVOPLASTY_URL && \
    cd $NOVOPLASTY_SRC && \
    git checkout $NOVOPLASTY_TAG && \
    cp NOVOPlasty3.7.pl $INSTALL_PFX/bin/NOVOPlasty.pl && \
    cd .. && \
    rm -rf $NOVOPLASTY_SRC

##### Configuration ############################################################

USER biodocker
WORKDIR /data
# CMD ['Various']

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
