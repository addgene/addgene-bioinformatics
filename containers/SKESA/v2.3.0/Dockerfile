##### Base image ###############################################################

FROM ubuntu:xenial

##### Metadata #################################################################

LABEL base.image="ubuntu:xenial"
LABEL version="1.0.0"
LABEL software="SKESA"
LABEL software.version="2.3.0"
LABEL description="Strategic Kmer Extension for Scrupulous Assemblies"
LABEL website="https://github.com/ncbi/SKESA"
LABEL lisence="https://github.com/ncbi/SKESA/blob/master/LICENSE"
LABEL maintainer="Curtis Kapsak"
LABEL maintainer.email="pjx8@cdc.gov"

##### Installation #############################################################

RUN apt-get update && apt-get install -y wget

# Install SKESA
RUN mkdir skesa && \
  cd skesa && \
  wget https://github.com/ncbi/SKESA/releases/download/v2.3.0/skesa.centos6.10 && \
  mv skesa.centos6.10 skesa && \
  chmod +x skesa

##### Configuration ############################################################

RUN mkdir /data
WORKDIR /data
ENV PATH=${PATH}:/skesa

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
