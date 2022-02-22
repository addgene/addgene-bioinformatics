# Introduction

The addgene-bioinformatics repository demonstrates a two step workflow
for assembling plasmid next generation sequencing data using Toil.

# Set-up

## Install Docker Desktop, and start the Docker daemon

Download the Docker Desktop package, then install as usual. The Docker
daemon should be started by default.

Docker will need to mount the directory Toil creates for a job. When
using Docker Desktop on macOS, this mount is accomplished in
"Preferences > Resources > File Sharing" by adding "/var/folders" to the list.

Additionally, Docker will need to be allocated enough resources to run a job.
This can be done in "Preferences > Resources > Advanced" by setting (minimum)
CPUs to 2 and memory to 4.0GB

## Clone the repository:

```shell
$ git clone git@github.com:addgene/addgene-bioinformatics.git
```

## Create a virtual environment (using a Python 3 version), and install the requirements:

```shell
$ mkvirtualenv addgene-bioinformatics
$ pip install -r requirements.txt
```

## Build Required Docker Images:

```shell
$ cd containers
$ ./build.sh
```

# Demonstration

## Run the tests:

```shell
$ cd src/python/jobs
$ python JobsTest.py
```

## Run one of the jobs with default inputs locally:

```shell
$ python src/python/jobs/SpadesJob.py sjfs
$ python src/python/jobs/ApcJob.py ajfs
$ python src/python/jobs/WellAssemblyJob.py wajfs
$ python src/python/jobs/PlateAssemblyJob.py pajfs
```

## Run a well or sample plate assembly job locally with data imported from S3:

```shell
$ python src/python/jobs/WellAssemblyJob.py  -s s3 -d addgene-sequencing-data/2018/FASTQ -l A11935X_sW0148 -w A01 wajfs
$ python src/python/jobs/PlateAssemblyJob.py -s s3 -d addgene-sequencing-data/2018/FASTQ -l A11935X_sW0148 pajfs
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
$ python WellAssemblyJob.py  -s s3 -d addgene-sequencing-data/2018/FASTQ -l A11935X_sW0148 -w A01 wajfs
$ python PlateAssemblyJob.py -s s3 -d addgene-sequencing-data/2018/FASTQ -l A11935X_sW0148 pajfs
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
