##### Base image ###############################################################

FROM ubuntu:xenial

##### Metadata #################################################################

LABEL base.image="ubuntu:xenial"
LABEL version="1"
LABEL software="Unicycler"
LABEL software.version="0.4.7"
LABEL description="Unicycler is an assembly pipeline for bacterial genomes"
LABEL website="https://github.com/rrwick/Unicycler"
LABEL license="https://github.com/rrwick/Unicycler#license"
LABEL maintainer="Kelsey Florek"
LABEL maintainer.email="kelsey.florek@slh.wisc.edu"

##### Installation #############################################################

RUN apt-get update && apt-get install -y software-properties-common && \
 add-apt-repository -y ppa:openjdk-r/ppa && apt-get update && apt-get install -y \
 python3 \
 python3-pip \
 libz-dev \
 bowtie2 \
 ncbi-blast+ \
 openjdk-11-jre \
 wget \
 bzip2 \
 libncurses5-dev \
 libncursesw5-dev \
 libbz2-dev \
 liblzma-dev \
 cmake

# Install Unicycler
RUN wget https://github.com/rrwick/Unicycler/archive/v0.4.7.tar.gz && \
 mkdir unicycler && \
 tar -xzvf v0.4.7.tar.gz -C unicycler --strip-components 1 && \
 rm v0.4.7.tar.gz

# Install Pilon
RUN wget https://github.com/broadinstitute/pilon/releases/download/v1.22/pilon-1.22.jar -P /pilon

# Install SPAdes
RUN wget http://cab.spbu.ru/files/release3.13.0/SPAdes-3.13.0-Linux.tar.gz && \
 mkdir spades && \
 tar -xzvf SPAdes-3.13.0-Linux.tar.gz -C spades --strip-components 1 && \
 rm SPAdes-3.13.0-Linux.tar.gz

# Install Samtools and BCFtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
 wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
 mkdir samtools && \
 tar -vxjf samtools-1.9.tar.bz2 -C samtools --strip-components 1 && \
 rm samtools-1.9.tar.bz2 && \
 mkdir bcftools && \
 tar -vxjf bcftools-1.9.tar.bz2 -C bcftools --strip-components 1 && \
 rm bcftools-1.9.tar.bz2 && \
 cd samtools && \
 make && \
 make install && \
 cd ../bcftools && \
 make && \
 make install

# Install Racon
RUN wget https://github.com/isovic/racon/releases/download/1.3.1/racon-v1.3.1.tar.gz && \
 mkdir racon && \
 tar -xzvf racon-v1.3.1.tar.gz -C racon --strip-components 1 && \
 rm racon-v1.3.1.tar.gz && \
 cd racon && \
 mkdir build && \
 cd build && \
 cmake -DCMAKE_BUILD_TYPE=Release ../ && \
 make && \
 make install

# Install Miniasm
RUN wget https://github.com/lh3/miniasm/archive/v0.3.tar.gz && \
 mkdir miniasm && \
 tar -xzvf v0.3.tar.gz -C miniasm --strip-components 1 && \
 rm v0.3.tar.gz && \
 cd miniasm && \
 make

# Install Unicycler
RUN python3 /unicycler/setup.py install

##### Configuration ############################################################

RUN mkdir /data
WORKDIR /data
ENV PATH="/pilon:${PATH}"
ENV PATH="/spades/bin:${PATH}"
ENV PATH="/miniasm:${PATH}"
# CMD ['unicycler.py']

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
