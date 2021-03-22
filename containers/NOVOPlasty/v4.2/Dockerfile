##### Base image ###############################################################

FROM ralatsdio/biocontainers:v1.2.0

##### Metadata #################################################################

LABEL base_image="ralatsdio/biocontainers"
LABEL version="1.2.0"
LABEL software="NOVOPlasty"
LABEL software.version="v4.2"
LABEL about.summary="The organelle assembler and heteroplasmy caller"
LABEL about.home="https://github.com/ndierckx/NOVOPlasty"
LABEL about.provides="NOVOPlasty.pl v4.2"
LABEL about.license_file="https://github.com/ndierckx/NOVOPlasty/blob/master/LICENSE"
LABEL about.license="Custom"
LABEL about.tags="Genomics"
LABEL extra.identifiers.biotools="NOVOPlasty"
LABEL extra.binaries="NOVOPlasty.pl"

##### Installation #############################################################

USER root
ENV INSTALL_PFX=/home/biodocker

# Install NOVOPlasty
ENV NOVOPLASTY_VER=4.2
RUN cd $INSTALL_PFX && \
    git clone https://github.com/ndierckx/NOVOPlasty && \
    cd NOVOPlasty && \
    git checkout NOVOPlasty$NOVOPLASTY_VER && \
    cp NOVOPlasty$NOVOPLASTY_VER.pl $INSTALL_PFX/bin/NOVOPlasty.pl && \
    cd .. && \
    rm -rf NOVOPlasty

##### Configuration ############################################################

USER biodocker
WORKDIR /data
# CMD ['NOVOPlasty.pl']

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
