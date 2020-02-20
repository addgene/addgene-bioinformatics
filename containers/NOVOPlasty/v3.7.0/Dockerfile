##### Base image ###############################################################

FROM ralatsdio/biocontainers:v1.1.0

##### Metadata #################################################################

LABEL base_image="ralatsdio/biocontainers"
LABEL version="1.1.1"
LABEL software="NOVOPlasty"
LABEL software.version="v3.7.0"
LABEL about.summary="The organelle assembler and heteroplasmy caller"
LABEL about.home="https://github.com/ndierckx/NOVOPlasty"
LABEL about.provides="NOVOPlasty.pl v3.7.0"
LABEL about.license_file="https://github.com/ndierckx/NOVOPlasty/blob/master/LICENSE"
LABEL about.license="Custom"
LABEL about.tags="Genomics"
LABEL extra.identifiers.biotools="NOVOPlasty"
LABEL extra.binaries="NOVOPlasty.pl"

##### Installation #############################################################

USER root
ENV INSTALL_PFX=/home/biodocker

# Install NOVOPlasty
ENV NOVOPLASTY_URL=https://github.com/ndierckx/NOVOPlasty
ENV NOVOPLASTY_SRC=NOVOPlasty
ENV NOVOPLASTY_TAG=NOVOPlasty3.7
RUN cd $INSTALL_PFX && \
    git clone $NOVOPLASTY_URL && \
    cd $NOVOPLASTY_SRC && \
    git checkout $NOVOPLASTY_TAG && \
    cp NOVOPlasty3.7.pl $INSTALL_PFX/bin/NOVOPlasty.pl && \
    cd .. && \
    rm -rf $NOVOPLASTY_SRC

##### Configuration ############################################################

USER biodocker
WORKDIR /data
# CMD ['NOVOPlasty.pl']

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
