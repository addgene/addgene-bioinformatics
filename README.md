# Introduction

The addgene-bioinformatics repository demonstrates a two step workflow for assembling plasmid next generation sequencing data using Toil.

# Set-up

## Install Docker Desktop, and start the Docker daemon

Docker will need to mount the directory Toil creates for a job. When using Docker Desktop on macOS, this mount is accomplished in Preferences > File Sharing by adding /var/folders to the list.

## Clone the repository:

```
$ git clone git@github.com:addgene/addgene-bioinformatics.git
```

## Create a virtual environment (using a Python 3 version), and install the requirements:

```
$ mkvirtualenv addgene-bioinformatics
$ pip install -r requirements.txt
```

## Pull required images from Docker Hub:

```
$ docker pull biocontainers/spades:v3.13.1_cv1
$ docker pull ralatsdio/apc:latest
```

# Demonstration

## Run one of the Toil jobs:

```
$ python src/python/SpadesJob.py spadesJobStore
$ python src/python/ApcJob.py apcJobStore
$ python src/python/WellAssemblyJob.py wellAssemblyJobStore
$ python src/python/PlateAssemblyJob.py plateAssemblyJobStore
```
