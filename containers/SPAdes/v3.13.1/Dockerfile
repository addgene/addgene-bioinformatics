##### Base image ###############################################################

FROM biocontainers/biocontainers:v1.0.0_cv4

##### Metadata #################################################################

LABEL base_image="biocontainers:v1.0.0_cv4"
LABEL version="1"
LABEL software="SPAdes"
LABEL software.version="3.13.1"
LABEL about.summary="SPAdes – St. Petersburg genome assembler – is an assembly toolkit containing various assembly pipelines."
LABEL about.home="https://github.com/ablab/spades"
LABEL about.provides="spades.pl v3.13.1"
LABEL about.license_file="https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html"
LABEL about.license="GPLv2"
LABEL about.tags="Genomics"
LABEL extra.identifiers.biotools="SPAdes"

##### Installation #############################################################

USER root
ENV INSTALL_PFX=/home/biodocker

# SPAdes
ENV ZIP=SPAdes-3.13.1-Linux.tar.gz
ENV URL=http://cab.spbu.ru/files/release3.13.1
ENV FOLDER=SPAdes-3.13.1-Linux
ENV DST=/home/biodocker/bin
RUN wget $URL/$ZIP -O $DST/$ZIP && \
    tar xvf $DST/$ZIP -C $DST && \
    rm $DST/$ZIP && \
    cd $DST/$FOLDER && \
    mv $DST/$FOLDER/bin/* /home/biodocker/bin && \
    mv $DST/$FOLDER/share /home/biodocker && \
    rm -rf $DST/$FOLDER

##### Configuration ############################################################

USER biodocker
WORKDIR /data
# CMD ['NOVOPlasty.pl']

##### Maintainer ###############################################################

MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
