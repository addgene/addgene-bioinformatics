set -x

# Build base images
docker build --tag ralatsdio/biocontainers:v1.2.0 BioContainers/v1.2.0 \
     2>&1 | tee BioContainers/v1.2.0/build.log

# Build assembler images
docker build --tag ralatsdio/idba:v1.1.4 IDBA/v1.1.4 \
     2>&1 | tee IDBA/v1.1.4/build.log
docker build --tag ralatsdio/masurca:v4.0.1 MaSuRCA/v4.0.1 \
     2>&1 | tee MaSuRCA/v4.0.1/build.log
docker build --tag ralatsdio/novoplasty:v4.2 NOVOPlasty/v4.2 \
     2>&1 | tee NOVOPlasty/v4.2/build.log
docker build --tag ralatsdio/skesa:v2.4.0 SKESA/v2.4.0 \
     2>&1 | tee SKESA/v2.4.0/build.log
docker build --tag ralatsdio/spades:v3.15.1 SPAdes/v3.15.1 \
     2>&1 | tee SPAdes/v3.15.1/build.log
docker build --tag ralatsdio/shovill:v1.1.0 Shovill/v1.1.0 \
     2>&1 | tee Shovill/v1.1.0/build.log
docker build --tag ralatsdio/unicycler:v0.4.8 Unicycler/v0.4.8 \
     2>&1 | tee Unicycler/v0.4.8/build.log

# Build tool images
docker build --tag ralatsdio/bbtools:v38.90 BBTools/v38.90 \
     2>&1 | tee BBTools/v38.90/build.log
docker build --tag ralatsdio/csc:v0.1.0 CSC/v0.1.0 \
     2>&1 | tee CSC/v0.1.0/build.log
docker build --tag ralatsdio/kmc:v3.1.1 KMC/v3.1.1 \
     2>&1 | tee KMC/v3.1.1/build.log
docker build --tag ralatsdio/repdenovo:v0.1.0 REPdenovo/v0.1.0 \
     2>&1 | tee REPdenovo/v0.1.0/build.log
docker build --tag ralatsdio/spade:v0.1.0 SPADE/v0.1.0 \
     2>&1 | tee SPADE/v0.1.0/build.log
docker build --tag ralatsdio/ssake:v4.0.1 SSAKE/v4.0.1 \
     2>&1 | tee SSAKE/v4.0.1/build.log
docker build --tag ralatsdio/samtools:v1.11 Samtools/v1.11 \
     2>&1 | tee Samtools/v1.11/build.log
docker build --tag ralatsdio/apc:v0.1.0 apc/v0.1.0 \
     2>&1 | tee apc/v0.1.0/build.log

# Summarize build logs
find . -name 'build.log' -exec echo ">{}" \; -exec tail -2 {} \;
