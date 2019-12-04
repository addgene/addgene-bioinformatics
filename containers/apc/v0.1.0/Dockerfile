##### Base image ###############################################################

FROM ralatsdio/biocontainers:v1.1.0

##### Metadata #################################################################

LABEL base_image="ralatsdio/biocontainers"
LABEL version="1.1.0"
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
ENV LAST_URL=http://last.cbrc.jp/last/
ENV LAST_SRC=last
ENV LAST_TAG=992
RUN hg clone $LAST_URL && \
    cd $LAST_SRC && \
    hg update $LAST_TAG && \
    make && \
    make install prefix=$INSTALL_PFX && \
    cd .. && \
    rm -rf $LAST_SRC

# Install BWA
ENV BWA_URL=https://github.com/lh3/bwa.git
ENV BWA_SRC=bwa
ENV BWA_TAG=v0.7.17
RUN git clone $BWA_URL && \
    cd $BWA_SRC && \
    git checkout $BWA_TAG && \
    make && \
    cp bwa $INSTALL_PFX/bin && \
    cd .. && \
    rm -rf $BWA_SRC

# Install HTSlib (Use BioContainers install)
ENV HTSLIB_ZIP=htslib-1.2.1.tar.bz2
ENV HTSLIB_URL=https://github.com/BioDocker/software-archive/releases/download/htslib
ENV HTSLIB_SRC=htslib-1.2.1
RUN wget $HTSLIB_URL/$HTSLIB_ZIP -O $HTSLIB_ZIP && \
    tar xvf $HTSLIB_ZIP && \
    rm $HTSLIB_ZIP && \
    cd $HTSLIB_SRC && \
    make && \
    make install && \
    cd .. && \
    rm -rf $HTSLIB_SRC

# Install Samtools (Use BioContainers install)
ENV SAMTOOLS_VER=1.7.0
RUN conda install ncurses samtools=$SAMTOOLS_VER

# Install apc (Use ralatsdc fork)
ENV APC_URL=https://github.com/ralatsdc/apc.git
ENV APC_SRC=apc
ENV APC_TAG=v0.1.0
RUN git clone $APC_URL && \
    cd $APC_SRC && \
    git checkout $APC_TAG && \
    cp apc.pl $INSTALL_PFX/bin && \
    cd .. && \
    rm -rf $APC_SRC

##### Configuration ############################################################

USER biodocker
WORKDIR /data
# CMD ['apc.pl']

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
