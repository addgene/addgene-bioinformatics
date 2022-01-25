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

```shell
$ git clone git@github.com:addgene/addgene-bioinformatics.git
```

## Create a virtual environment (using a Python 3 version), and install the requirements:

```shell
$ mkvirtualenv addgene-bioinformatics
$ pip install -r requirements.txt
```

## Pull required images from Docker Hub:

```shell
$ docker pull ralatsdio/masurca:v3.3.1
$ docker pull ralatsdio/novoplasty:v3.7.0
$ docker pull ralatsdio/shovill:v1.0.9
$ docker pull ralatsdio/skesa:v2.3.0
$ docker pull ralatsdio/spades:v3.13.1
$ docker pull ralatsdio/unicycler:v0.4.7
$ docker pull ralatsdio/apc:v0.1.0
```

# Demonstration

## Run the tests:

```shell
$ python src/python/toil/JobsTest.py
```

## Run one of the jobs with default inputs locally:

```shell
$ python src/python/toil/SpadesJob.py sjfs
$ python src/python/toil/ApcJob.py ajfs
$ python src/python/toil/WellAssemblyJob.py wajfs
$ python src/python/toil/PlateAssemblyJob.py pajfs
```

## Run a well or sample plate assembly job locally with data imported from S3:

```shell
$ python src/python/toil/WellAssemblyJob.py  -s s3 -d addgene-sequencing-data/2018/FASTQ -p A11935X_sW0148 -w A01 wajfs
$ python src/python/toil/PlateAssemblyJob.py -s s3 -d addgene-sequencing-data/2018/FASTQ -p A11935X_sW0148 pajfs
```

## Run one of the jobs on EC2

The following assumes the instructions for [preparing your AWS
environment](https://toil.readthedocs.io/en/latest/running/cloud/amazon.html#preparing-your-aws-environment)
have been completed.

### Launch the cluster leader:

```shell
$ toil launch-cluster --zone us-east-1a --keyPairName id_rsa --leaderNodeType t2.medium assembly-cluster
```

### Synchronize code and data:

```shell
$ sh src/sh/make-archives-for-leader.sh
$ toil rsync-cluster --zone us-east-1a assembly-cluster python.tar.gz :/root
$ toil rsync-cluster --zone us-east-1a assembly-cluster miscellaneous.tar.gz :/root
```

### Login to the cluster leader, and extract the archives:

```shell
$ toil ssh-cluster --zone us-east-1a assembly-cluster
$ tar -zxvf python.tar.gz
$ tar -zxvf miscellaneous.tar.gz
```

### Run the default plate assembly job on the cluster leader only with a local or S3 file store:

```shell
$ cd python/src/toil
$ python PlateAssemblyJob.py --data-directory miscellaneous --plate-spec A11967B_sW0154 pajfs
$ python PlateAssemblyJob.py --data-directory miscellaneous --plate-spec A11967B_sW0154 aws:us-east-1:pajfs
```

## Run a well, or sample plate assembly job on the cluster leader only with data imported from S3:

```bash
$ cd python/src/toil
$ python WellAssemblyJob.py  -s s3 -d addgene-sequencing-data/2018/FASTQ -p A11935X_sW0148 -w A01 wajfs
$ python PlateAssemblyJob.py -s s3 -d addgene-sequencing-data/2018/FASTQ -p A11935X_sW0148 pajfs
```

### Run the default or a larger plate assembly job using auto-scaling with an S3 file store:

```shell
$ cd python/src/toil
$ python PlateAssemblyJob.py --data-directory miscellaneous --plate-spec A11967B_sW0154 --provisioner aws --nodeTypes c3.large --maxNodes 2 --batchSystem mesos aws:us-east-1:pajfs
$ python PlateAssemblyJob.py --data-directory miscellaneous --plate-spec A11967A_sW0154 --provisioner aws --nodeTypes c3.large --maxNodes 2 --batchSystem mesos aws:us-east-1:pajfs
```

### Destroy the cluster leader:

```shell
$ toil destroy-cluster --zone us-east-1a assembly-cluster
```

### Update python dependencies

```shell
pip install pip-tools
pip-compile -U
```
