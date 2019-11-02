############################### Dockerfile ##################################

FROM biocontainers/biocontainers:v1.0.0_cv4

################################# METADATA ###################################

LABEL base_image="biocontainers:v1.0.0_cv4"
LABEL version="3"
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

################## BEGIN INSTALLATION ######################

# Change user to root
USER root

# Shared environment
ENV INSTALL_PFX=/home/biodocker

# Install NOVOPlasty
ENV NOVOPlasty_URL=https://github.com/ndierckx/NOVOPlasty
ENV NOVOPlasty_SRC=apc
ENV NOVOPlasty_TAG=NOVOPlasty3.7
RUN git clone $NOVOPlasty_URL && \
    cd $NOVOPlasty_SRC && \
    git checkout $NOVOPlasty_TAG && \
    cp NOVOPlasty.pl $INSTALL_PFX/bin && \
    cd .. && \
    rm -rf $NOVOPlasty_SRC

# Change user to back to biodocker
USER biodocker

# Change WORKDIR to /data
WORKDIR /data

# Define default command
# CMD ['NOVOPlasty.pl']

##################### INSTALLATION END #####################

# File Author / Maintainer
MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
