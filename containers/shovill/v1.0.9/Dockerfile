############################### Dockerfile ##################################

FROM biocontainers/biocontainers:v1.0.0_cv4

################################# METADATA ###################################

LABEL base_image="biocontainers:v1.0.0_cv4"
LABEL version="3"
LABEL software="shovill"
LABEL software.version="v1.0.9"
LABEL about.summary="Assemble bacterial isolate genomes from Illumina paired-end reads"
LABEL about.home="https://github.com/tseemann/shovill"
LABEL about.provides="shovill v1.0.9"
LABEL about.license_file="https://github.com/tseemann/shovill/blob/master/LICENSE"
LABEL about.license="GNU General Public License v3.0"
LABEL about.tags="Genomics"
LABEL extra.identifiers.biotools="shovill"
LABEL extra.binaries="shovill"

################## BEGIN INSTALLATION ######################

# Change user to root
USER root

# Shared environment
ENV INSTALL_PFX=/home/biodocker

# Install shovill
RUN conda install -c conda-forge -c bioconda -c defaults shovill

# Change user to back to biodocker
USER biodocker

# CHANGE WORKDIR TO /DATA
WORKDIR /data

# DEFINE DEFAULT COMMAND
# CMD ['shovill']

##################### INSTALLATION END #####################

# File Author / Maintainer
MAINTAINER Raymond LeClair <raymond.leclair@springbok.io>
