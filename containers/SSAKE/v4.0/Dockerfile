##### Base image ###############################################################

FROM ralatsdio/biocontainers:v1.1.0

##### Metadata #################################################################

LABEL base.image="ralatsdio/biocontainers"
LABEL version="1.1.0"
LABEL software="SSAKE"
LABEL software.version="v4.0"
LABEL about.summary="Genome assembly with short sequence reads"
LABEL about.home="https://github.com/bcgsc/SSAKE"
LABEL about.provides="ssake v4.0"
LABEL about.license_file="https://github.com/bcgsc/SSAKE/blob/master/LICENSE"
LABEL about.license="GNU General Public License v2.0"
LABEL about.tags="Genomics"
LABEL extra.identifiers.biotools="SSAKE"
LABEL extra.binaries=""

##### Installation #############################################################

USER root
ENV INSTALL_PFX=/home/biodocker

# Install Perl modules
RUN cpan install Statistics::Descriptive

# Install SSAKE
ENV SSAKE_VER=4.0
ENV SSAKE_TGZ=ssake_v4-0.tar.gz
RUN cd $INSTALL_PFX && \
    wget https://github.com/bcgsc/SSAKE/releases/download/v$SSAKE_VER/$SSAKE_TGZ && \
    tar -zxvf $SSAKE_TGZ && \
    rm $SSAKE_TGZ && \
    cp ssake/SSAKE $INSTALL_PFX/bin/ && \
    cp ssake/tools/* $INSTALL_PFX/bin/ && \
    cat ssake/tools/runSSAKE.sh \
        | sed "s%/usr/bin/bash%/bin/bash%" \
        | sed "s%\.\./tools/%%" \
        | sed "s%/usr/bin/time -v -o \$4\.time \.\./%%" \
        > $INSTALL_PFX/bin/runSSAKE.sh && \
    cat ssake/tools/getStats.pl \
        | sed "s%/home/martink/bin%/usr/bin%" \
        > $INSTALL_PFX/bin/getStats.pl && \
    rm -rf ssake

##### Configuration ############################################################

USER biodocker
WORKDIR /data
# CMD ['ssake']

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
