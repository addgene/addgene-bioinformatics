##### Base image ###############################################################

FROM ubuntu:xenial

##### Metadata #################################################################

LABEL base.image="ubuntu:xenial"
LABEL version="1.0.0"
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

RUN apt-get update && \
    apt-get install -y \
        openjdk-9-jre \
        wget

# Install BBTools
RUN mkdir bbtools && \
  cd bbtools && \
  wget https://sourceforge.net/projects/bbmap/files/BBMap_38.73.tar.gz && \
  tar -zxvf BBMap_38.73.tar.gz

##### Configuration ############################################################

RUN mkdir /data
WORKDIR /data
ENV PATH=${PATH}:/bbtools/bbmap
# CMD ['Various']

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
