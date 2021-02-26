##### Base image ###############################################################

FROM ralatsdio/biocontainers:v1.2.0

##### Metadata #################################################################

LABEL base.image="ralatsdio/biocontainers"
LABEL version="1.2.0"
LABEL software="KMC"
LABEL software.version="v3.1.1"
LABEL about.summary="A disk-based programm for counting k-mers from (possibly gzipped) FASTQ/FASTA files"
LABEL about.home="https://github.com/refresh-bio/KMC"
LABEL about.provides="kmc v3.1.1"
LABEL about.license_file="https://www.gnu.org/licenses/gpl-3.0.en.html"
LABEL about.license="GNU General Public License v3.0"
LABEL about.tags="Genomics"
LABEL extra.identifiers.biotools="KMC"
LABEL extra.binaries="kmc, kmc_dump, kmc_tools"

##### Installation #############################################################

USER root
ENV INSTALL_PFX=/home/biodocker

# Install KMC
ENV KMC_VER=3.1.1
RUN cd $INSTALL_PFX && \
    wget https://github.com/refresh-bio/KMC/releases/download/v$KMC_VER/KMC$KMC_VER.linux.tar.gz && \
    tar -zxvf KMC$KMC_VER.linux.tar.gz && \
    rm KMC$KMC_VER.linux.tar.gz && \
    mv kmc* $INSTALL_PFX/bin/

##### Configuration ############################################################

USER biodocker
WORKDIR /data
# CMD ['kmc']

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
