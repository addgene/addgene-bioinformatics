##### Base image ###############################################################

FROM ralatsdio/biocontainers:v1.2.0

##### Metadata #################################################################

LABEL base_image="ralatsdio/biocontainers"
LABEL version="1.2.0"
LABEL software="apc"
LABEL software.version="v0.1.0"
LABEL about.summary="Perl script to test each sequence in a multi-fasta file for overlapping ends"
LABEL about.home="https://github.com/ralatsdc/apc.git"
LABEL about.provides="apc.pl v0.1.0"
LABEL about.license_file="NA"
LABEL about.license="NA"
LABEL about.tags="Genomics"
LABEL extra.identifiers.biotools="apc"
LABEL extra.binaries="apc.pl"

##### Installation #############################################################

USER root
ENV INSTALL_PFX=/home/biodocker

# Install LAST
ENV LAST_VER=992
RUN cd $INSTALL_PFX && \
    hg clone http://last.cbrc.jp/last/ && \
    cd last && \
    hg update $LAST_VER && \
    make && \
    make install prefix=$INSTALL_PFX && \
    cd .. && \
    rm -rf last

# Install BWA
ENV BWA_VER=v0.7.17
RUN cd $INSTALL_PFX && \
    git clone https://github.com/lh3/bwa.git && \
    cd bwa && \
    git checkout $BWA_VER && \
    make && \
    cp bwa $INSTALL_PFX/bin && \
    cd .. && \
    rm -rf bwa

# Install HTSlib (Use BioContainers install)
ENV HTSLIB_VER=1.2.1
RUN cd $INSTALL_PFX && \
    wget https://github.com/BioDocker/software-archive/releases/download/htslib/htslib-$HTSLIB_VER.tar.bz2 && \
    tar -jxvf htslib-$HTSLIB_VER.tar.bz2 && \
    rm htslib-$HTSLIB_VER.tar.bz2 && \
    cd htslib-$HTSLIB_VER && \
    make && \
    make install && \
    cd .. && \
    rm -rf htslib-$HTSLIB_VER

# Install Samtools (Use BioContainers install)
ENV SAMTOOLS_VER=1.7.0
RUN conda install ncurses samtools=$SAMTOOLS_VER

# Install apc (Use ralatsdc fork)
ENV APC_VER=v0.1.0
RUN cd $INSTALL_PFX && \
    git clone https://github.com/ralatsdc/apc.git && \
    cd apc && \
    git checkout $APC_VER && \
    cp apc.pl $INSTALL_PFX/bin && \
    cd .. && \
    rm -rf apc

##### Configuration ############################################################

USER biodocker
WORKDIR /data
# CMD ['apc.pl']

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
