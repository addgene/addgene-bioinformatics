##### Base image ###############################################################

FROM ralatsdio/biocontainers:v1.1.0

##### Metadata #################################################################

LABEL base_image="ralatsdio/biocontainers"
LABEL version="1.1.0"
LABEL software="CSC"
LABEL software.version="v0.1.0"
LABEL about.summary="Finds the cyclic rotation of x that minimises the blockwise q-gram distance from y"
LABEL about.home="https://github.com/ralatsdc/csc.git"
LABEL about.provides="csc v0.1.0"
LABEL about.license_file="https://github.com/ralatsdc/csc/blob/master/LICENSE"
LABEL about.license="GNU General Public License v3.0"
LABEL about.tags="Genomics"
LABEL extra.identifiers.biotools="csc"
LABEL extra.binaries="csc"

##### Installation #############################################################

USER root
ENV INSTALL_PFX=/home/biodocker

# Install apc (Use ralatsdc fork)
ENV CSC_VER=v0.1.0
RUN cd $INSTALL_PFX && \
    git clone https://github.com/ralatsdc/csc.git && \
    cd csc && \
    git checkout $CSC_VER && \
    ./pre-install.sh && \
    make -f Makefile.64-bit.gcc && \
    ln -s $INSTALL_PFX/csc/csc $INSTALL_PFX/bin

##### Configuration ############################################################

USER biodocker
WORKDIR /data
# CMD ['csc']

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
