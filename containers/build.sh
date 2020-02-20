
# Build base image
docker build --tag ralatsdio/biocontainers:v1.1.0 BioContainers/v1.1.0 \
       | tee BioContainers/v1.1.0/build.log 2>&1

# Build assembler images
docker build --tag ralatsdio/masurca:v3.3.1 MaSuRCA/v3.3.1 \
       | tee MaSuRCA/v3.3.1/build.log 2>&1
docker build --tag ralatsdio/novoplasty:v3.7.0 NOVOPlasty/v3.7.0 \
       | tee NOVOPlasty/v3.7.0/build.log 2>&1
docker build --tag ralatsdio/skesa:v2.3.0 SKESA/v2.3.0 \
       | tee SKESA/v2.3.0/build.log 2>&1
docker build --tag ralatsdio/spades:v3.13.1 SPAdes/v3.13.1 \
       | tee SPAdes/v3.13.1/build.log 2>&1
docker build --tag ralatsdio/shovill:v1.0.9 Shovill/v1.0.9 \
       | tee Shovill/v1.0.9/build.log 2>&1
docker build --tag ralatsdio/unicycler:v0.4.7 Unicycler/v0.4.7 \
       | tee Unicycler/v0.4.7/build.log 2>&1

# Build tool images
docker build --tag ralatsdio/bbtools:v38.73 BBTools/v38.73 \
       | tee BBTools/v38.73/build.log 2>&1
docker build --tag ralatsdio/repdenovo:v0.1.0 REPdenovo/v0.1.0 \
       | tee REPdenovo/v0.1.0/build.log 2>&1
docker build --tag ralatsdio/apc:v0.1.0 apc/v0.1.0 \
       | tee apc/v0.1.0/build.log 2>&1
docker build --tag ralatsdio/csc:v0.1.0 csc/v0.1.0 \
       | tee csc/v0.1.0/build.log 2>&1

# Summarize build logs
find . -name 'build.log' -exec echo ">{}" \; -exec tail -2 {} \;
