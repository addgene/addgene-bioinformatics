# Introduction

The addgene-bioinformatics repository demonstrates a two step workflow
for assembling plasmid next generation sequencing data using Toil.

# Set-up

## Install Docker Desktop, and start the Docker daemon

Download the Docker Desktop package, then install as usual. The Docker
daemon should be started by default.

Docker will need to mount the directory Toil creates for a job. When
using Docker Desktop on macOS, this mount is accomplished in
"Preferences > File Sharing" by adding "/var/folders" to the list.

## Clone the repository:

```
$ git clone git@github.com:addgene/addgene-bioinformatics.git
```

## Create a virtual environment (using a Python 2 version), and install the requirements:

```
$ mkvirtualenv addgene-bioinformatics
$ pip install -r requirements.txt
```

Note that Toil runs fine using a python 3 version, however, the Toil
appliance (a docker image) still uses a python 2 version.

## Pull required images from Docker Hub:

```
$ docker pull biocontainers/spades:v3.13.1_cv1
$ docker pull ralatsdio/apc:latest
```

# Demonstration

## Run the tests:

```
$ python src/python/JobsTest.py
```

## Run one of the jobs with default inputs locally:

```
$ python src/python/SpadesJob.py sjs
$ python src/python/ApcJob.py ajs
$ python src/python/WellAssemblyJob.py wajs
$ python src/python/PlateAssemblyJob.py pajs
```

## Run a well, or sample plate assembly job locally with data imported from S3:

```
$ python src/python/WellAssemblyJob.py  -s s3 -d addgene-sequencing-data/2018/FASTQ -p A11935X_sW0148 -w A01 wajs
$ python src/python/PlateAssemblyJob.py -s s3 -d addgene-sequencing-data/2018/FASTQ -p A11935X_sW0148 pajs
```

## Run one of the jobs on EC2

The following assumes the instructions for [preparing your AWS
environment](https://toil.readthedocs.io/en/latest/running/cloud/amazon.html#preparing-your-aws-environment)
have been completed.

### Launch the cluster leader:

```
$ toil launch-cluster --zone us-east-1a --keyPairName id_rsa --leaderNodeType t2.medium assembly-cluster
```

### Synchronize code and data:

```
$ src/sh/make-archives-for-leader.sh
$ toil rsync-cluster --zone us-east-1a assembly-cluster python.tar.gz :/root
$ toil rsync-cluster --zone us-east-1a assembly-cluster miscellaneous.tar.gz :/root
```

### Login to the cluster leader, and extract the archives:

```
$ toil ssh-cluster --zone us-east-1a assembly-cluster
# cd
# tar -zxvf python.tar.gz
# tar -zxvf miscellaneous.tar.gz
```

### Run the Toil sort example:

```
# python sort.py --provisioner aws --nodeTypes c3.large --maxNodes 2 --batchSystem mesos aws:us-east-1:sort
```

### Run the default plate assembly job on the cluster leader only with a local or S3 file store:

```
# python PlateAssemblyJob.py --data-directory miscellaneous --plate-spec A11967B_sW0154 pajs
# python PlateAssemblyJob.py --data-directory miscellaneous --plate-spec A11967B_sW0154 aws:us-east-1:pajs
```

## Run a well, or sample plate assembly job on the cluster leader only with data imported from S3:

```
# python WellAssemblyJob.py  -s s3 -d addgene-sequencing-data/2018/FASTQ -p A11935X_sW0148 -w A01 wajs
# python PlateAssemblyJob.py -s s3 -d addgene-sequencing-data/2018/FASTQ -p A11935X_sW0148 pajs
```

### Run the default or a larger plate assembly job using auto-scaling with an S3 file store:

```
# python PlateAssemblyJob.py --data-directory miscellaneous --plate-spec A11967B_sW0154 --provisioner aws --nodeTypes c3.large --maxNodes 2 --batchSystem mesos aws:us-east-1:pajs
# python PlateAssemblyJob.py --data-directory miscellaneous --plate-spec A11967A_sW0154 --provisioner aws --nodeTypes c3.large --maxNodes 2 --batchSystem mesos aws:us-east-1:pajs
```

### Destroy the cluster leader:

```
$ toil destroy-cluster --zone us-east-1a assembly-cluster
```
