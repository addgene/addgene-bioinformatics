##### Base image ###############################################################

FROM ubuntu:18.04

##### Metadata #################################################################

LABEL base_image="ubuntu:18.04"
LABEL software="MaSuRCA"
LABEL software.version="v3.3.1"
LABEL about.summary="The Maryland Super Read Cabog Assembler (MaSuRCA)"
LABEL about.home="https://github.com/alekseyzimin/masurca"
LABEL about.provides="masurca v3.3.1"
LABEL about.license_file="https://github.com/alekseyzimin/masurca/blob/master/LICENSE.txt"
LABEL about.license="GNU General Public License v3.0"
LABEL about.tags="Genomics"

##### Installation #############################################################

ARG MASURCA_VERSION=3.3.1

RUN apt-get update && \
    apt-get install -y
        wget \
        libboost-all-dev \
	libbz2-dev && \
    apt-get clean
		       
# Install MaSuRCA
WORKDIR /software
RUN wget https://github.com/alekseyzimin/masurca/releases/download/v${MASURCA_VERSION}/MaSuRCA-${MASURCA_VERSION}.tar.gz
RUN tar -xf *.tar.gz && \
    rm *.tar.gz && \
    mv * masurca
WORKDIR /software/masurca
RUN chmod +x install.sh && \
    ./install.sh && \
    chmod +x bin/masurca

##### Configuration ############################################################

RUN mkdir /data
WORKDIR /data
ENV PATH="/software/masurca/bin:${PATH}"

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
