#!/bin/bash

HOME="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd | sed s%/src/sh%% )"

pushd $HOME/src/python/jobs
tar -czvf python.tar.gz *.py Assembler.ini adapters.fa oriSeed.fasta output
mv python.tar.gz $HOME
popd

pushd $HOME/dat
tar -czvf miscellaneous.tar.gz \
    miscellaneous/A11967A_sW0154_FASTQ/A11967A_sW0154_A* \
    miscellaneous/A11967A_sW0154_FASTQ/A11967A_sW0154_B01_* \
    miscellaneous/A11967A_sW0154_FASTQ/A11967A_sW0154_G05_* \
    miscellaneous/A11967A_sW0154_FASTQ/A11967A_sW0154_G06_* \
    miscellaneous/A11967B_sW0154_FASTQ
mv miscellaneous.tar.gz $HOME
popd
