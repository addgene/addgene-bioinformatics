##### Base image ###############################################################

FROM ralatsdio/biocontainers:v1.2.0

##### Metadata #################################################################

LABEL base.image="ralatsdio/biocontainers"
LABEL container.version="1.2.0"
LABEL software="Shovill"
LABEL software.version="1.1.0"
LABEL about.summary="Faster than SPAdes de novo DBG genome assembler (with assembler options!)"
LABEL about.home="https://github.com/tseemann/shovill"
LABEL about.provides="shovill v1.1.0"
LABEL about.license_file="https://github.com/tseemann/shovill/blob/master/LICENSE"
LABEL about.license="?"
LABEL about.tags="Genomics"
LABEL extra.identifiers.biotools="Shovill"
LABEL extra.binaries="shovill"

##### Installation #############################################################

USER root
ENV INSTALL_PFX=/home/biodocker

# Install dependencies
RUN apt-get update && \
    apt-get install -y python \
        libbz2-dev \
        libcurl4-gnutls-dev \
        libfindbin-libs-perl \
        liblzma-dev \
        libncurses5-dev \
        libpthread-stubs0-dev \
        libssl-dev \
        openjdk-11-jre \
        pigz && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install SPAdes
ENV SPADES_VER=3.14.1
RUN cd $INSTALL_PFX && \
    wget http://cab.spbu.ru/files/release$SPADES_VER/SPAdes-$SPADES_VER-Linux.tar.gz && \
    tar -zxvf SPAdes-$SPADES_VER-Linux.tar.gz && \
    rm SPAdes-$SPADES_VER-Linux.tar.gz

# Install Seqtk
ENV SEQTK_VER=1.3
RUN cd $INSTALL_PFX && \
    mkdir seqtk && \ 
    cd seqtk && \
    wget https://github.com/lh3/seqtk/archive/v$SEQTK_VER.tar.gz && \
    tar -zxvf v$SEQTK_VER.tar.gz && \
    rm v$SEQTK_VER.tar.gz && \
    cd seqtk-$SEQTK_VER/ && \
    make && \
    make install

# Install Mash
ENV MASH_VER=2.1
RUN cd $INSTALL_PFX && \
    wget https://github.com/marbl/Mash/releases/download/v$MASH_VER/mash-Linux64-v$MASH_VER.tar && \
    tar -xvf mash-Linux64-v$MASH_VER.tar && \
    rm -rf mash-Linux64-v$MASH_VER.tar

# Install Lighter
ENV LIGHTER_VER=1.1.1
RUN cd $INSTALL_PFX && \
    wget https://github.com/mourisl/Lighter/archive/v$LIGHTER_VER.tar.gz && \
    tar -zxvf v$LIGHTER_VER.tar.gz && \
    rm -rf v$LIGHTER_VER.tar.gz && \
    cd Lighter-$LIGHTER_VER && \
    make

# Install Trimmomatic
ENV TRIMMOMATIC_VER=0.38
RUN cd $INSTALL_PFX && \
    mkdir trimmomatic && \
    cd trimmomatic && \
    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-$TRIMMOMATIC_VER.zip && \
    unzip Trimmomatic-$TRIMMOMATIC_VER.zip && \
    rm -rf Trimmomatic-$TRIMMOMATIC_VER.zip && \
    chmod +x Trimmomatic-$TRIMMOMATIC_VER/trimmomatic-$TRIMMOMATIC_VER.jar && \
    echo "#!/bin/bash" >> trimmomatic && \
    echo "exec java -jar $INSTALL_PFX/trimmomatic/Trimmomatic-$TRIMMOMATIC_VER/trimmomatic-$TRIMMOMATIC_VER.jar """"$""@"""" " >> trimmomatic && \
    chmod +x trimmomatic

# Install BWA-MEM
ENV BWA_VER=0.7.17
RUN cd $INSTALL_PFX && \
    mkdir bwa && \
    cd bwa && \
    wget https://github.com/lh3/bwa/releases/download/v$BWA_VER/bwa-$BWA_VER.tar.bz2 && \
    tar -jxvf bwa-$BWA_VER.tar.bz2 && \
    rm bwa-$BWA_VER.tar.bz2 && \
    cd bwa-$BWA_VER && \
    make

# Install Samtools 
ENV SAMTOOLS_VER=1.10
RUN cd $INSTALL_PFX && \
    mkdir samtools && \
    cd samtools && \
    wget https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VER/samtools-$SAMTOOLS_VER.tar.bz2 && \
    tar -jxvf samtools-$SAMTOOLS_VER.tar.bz2 && \
    rm samtools-$SAMTOOLS_VER.tar.bz2 && \
    cd samtools-$SAMTOOLS_VER && \
    ./configure && \
    make && \
    make install

# Install SKESA
ENV SKESA_VER=2.3.0
ENV CENTOS_VER=6.10
RUN cd $INSTALL_PFX && \
    mkdir skesa && \
    cd skesa && \  
    wget https://github.com/ncbi/SKESA/releases/download/v$SKESA_VER/skesa.centos$CENTOS_VER && \
    mv skesa.centos$CENTOS_VER $INSTALL_PFX/bin/skesa && \
    chmod +x $INSTALL_PFX/bin/skesa

# Install MEGAHIT
ENV MEGAHIT_VER=1.1.4
RUN cd $INSTALL_PFX && \
    mkdir megahit && \
    cd megahit && \
    wget https://github.com/voutcn/megahit/releases/download/v$MEGAHIT_VER/megahit_v${MEGAHIT_VER}_LINUX_CPUONLY_x86_64-bin.tar.gz && \
    tar -zxvf megahit_v${MEGAHIT_VER}_LINUX_CPUONLY_x86_64-bin.tar.gz && \
    rm megahit_v${MEGAHIT_VER}_LINUX_CPUONLY_x86_64-bin.tar.gz

# Install Velvet
ENV VELVET_VER=1.2.10
RUN cd $INSTALL_PFX && \
    mkdir velvet && \
    cd velvet && \
    wget https://github.com/dzerbino/velvet/archive/v$VELVET_VER.tar.gz && \
    tar -zxvf v$VELVET_VER.tar.gz && \
    rm -rf v$VELVET_VER.tar.gz && \
    cd velvet-$VELVET_VER && \
    make

# Install Flash
ENV FLASH_VER=1.2.11
RUN cd $INSTALL_PFX && \
    mkdir flash && \
    cd flash && \
    wget https://sourceforge.net/projects/flashpage/files/FLASH-$FLASH_VER.tar.gz && \
    tar -zxvf FLASH-$FLASH_VER.tar.gz && \
    rm -rf FLASH-$FLASH_VER.tar.gz && \
    cd FLASH-$FLASH_VER && \
    make

# Install Pilon
ENV PILON_VER=1.22
RUN cd $INSTALL_PFX && \
    mkdir pilon && \
    cd pilon && \
    wget https://github.com/broadinstitute/pilon/releases/download/v$PILON_VER/pilon-$PILON_VER.jar && \
    chmod +x pilon-$PILON_VER.jar && \
    echo "#!/bin/bash" >> pilon && \
    echo "exec java -jar $INSTALL_PFX/pilon/pilon-$PILON_VER.jar """"$""@"""" " >> pilon && \
    chmod +x pilon

# Install Samclip
RUN cd $INSTALL_PFX && \
    mkdir samclip && \
    cd samclip && \
    wget https://raw.githubusercontent.com/tseemann/samclip/master/samclip && \
    chmod +x samclip

# Install KMC
ENV KMC_VER=3.1.1
RUN cd $INSTALL_PFX && \
    wget https://github.com/refresh-bio/KMC/releases/download/v$KMC_VER/KMC$KMC_VER.linux.tar.gz && \
    tar -zxvf KMC$KMC_VER.linux.tar.gz && \
    rm KMC$KMC_VER.linux.tar.gz && \
    mv kmc* $INSTALL_PFX/bin/

# Install shovill
ENV SHOVILL_ENV=1.1.0
RUN cd $INSTALL_PFX && \
    mkdir shovill && \
    cd shovill && \
    wget https://github.com/tseemann/shovill/archive/v$SHOVILL_ENV.tar.gz && \
    tar -zxvf v$SHOVILL_ENV.tar.gz && \
    rm v$SHOVILL_ENV.tar.gz

##### Configuration ############################################################

USER biodocker
WORKDIR /data
# CMD ['shovill.pl']

# Set path to executables
# TODO: Sort which executables are required, then create link in $INSTALL_PFX/bin
ENV PATH="${PATH}:\
$INSTALL_PFX/SPAdes-$SPADES_VER-Linux/bin:\
$INSTALL_PFX/mash-Linux64-v$MASH_VER:\
$INSTALL_PFX/Lighter-$LIGHTER_VER:\
$INSTALL_PFX/trimmomatic:\
$INSTALL_PFX/bwa/bwa-$BWA_VER:\
$INSTALL_PFX/skesa:\
$INSTALL_PFX/megahit/megahit_v${MEGAHIT_VER}_LINUX_CPUONLY_x86_64-bin:\
$INSTALL_PFX/velvet/velvet-$VELVET_VER:\
$INSTALL_PFX/flash/FLASH-$FLASH_VER:\
$INSTALL_PFX/shovill/shovill-$SHOVILL_ENV/bin:\
$INSTALL_PFX/pilon:\
$INSTALL_PFX/samclip"

# Set Perl locale settings
ENV LC_ALL=C

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
