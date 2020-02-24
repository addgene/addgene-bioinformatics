##### Base image ###############################################################

FROM ralatsdio/biocontainers:v1.1.0

##### Metadata #################################################################

LABEL base_image="ralatsdio/biocontainers"
LABEL version="1.1.0"
LABEL software="SPADE"
LABEL software.version="v0.1.0"
LABEL about.summary="Explore various periodic repeat regions comprehensively from large genomic and protein data resources"
LABEL about.home="git@github.com:ralatsdc/SPADE.git"
LABEL about.provides="SPADE.py v0.1.0"
LABEL about.license_file="GNU General Public License v3.0"
LABEL about.license="https://github.com/ralatsdc/SPADE/blob/master/LICENSE"
LABEL about.tags="Genomics"
LABEL extra.identifiers.biotools="SPADE"
LABEL extra.binaries="SPADE.py"

##### Installation #############################################################

USER root
ENV INSTALL_PFX=/home/biodocker

# Install Python requirements
RUN pip install --upgrade pip && \
    pip install matplotlib==2.2.3 && \
    pip install seaborn==0.8.1 && \
    pip install weblogo==3.6.0 && \
    pip install biopython

# Install MAFFT 7.455
ENV MAFFT_VER=cf201901
RUN conda install -c bioconda/label/$MAFFT_VER mafft

# Install BLAST 2.9.0
ENV BLAST_VER=cf201901
RUN conda install -c bioconda/label/$BLAST_VER blast

# Install SPADE
ENV SPADE_VER=0.1.0
RUN cd $INSTALL_PFX && \
    wget https://github.com/ralatsdc/SPADE/archive/v$SPADE_VER.tar.gz -O SPADE-$SPADE_VER.tar.gz && \
    tar -zxvf SPADE-$SPADE_VER.tar.gz && \
    rm SPADE-$SPADE_VER.tar.gz && \
    mv SPADE-$SPADE_VER/* $INSTALL_PFX/bin && \
    rm -rf SPADE-$SPADE_VER

##### Configuration ############################################################

USER biodocker
WORKDIR /data
# CMD ['SPADE.py']

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
