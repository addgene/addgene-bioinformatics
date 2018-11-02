#!/bin/bash

export PATH="../../local/bin:$PATH"

fastq_home="../../dat"
fastq_dir="representative-assortment"

fastq_fnm="A12180A_B11"

spades_start=$(date +%s)
spades.py \
    -1 ${fastq_home}/${fastq_dir}/${fastq_fnm}_R1_001.fastq.gz \
    -2 ${fastq_home}/${fastq_dir}/${fastq_fnm}_R2_001.fastq.gz \
    -o ${fastq_fnm} --cov-cutoff 100 \
    | tee spades_${fastq_fnm}.out \
    | grep "^=="
spades_end=$(date +%s)
spades_real=$(expr $spades_end - $spades_start)

echo "FASTQ filename : ${fastq_fnm} : SPAdes required (s) : ${spades_real}"

{ time ../apc/apc.pl \
       -b ${fastq_fnm} ${fastq_fnm}/scaffolds.fasta \
       > apc_${fastq_fnm}.out; } > apc_${fastq_fnm}_time.out 2>&1
apc_real=$(grep "real" apc_${fastq_fnm}_time.out | cut -f 2)

echo "FASTQ filename : ${fastq_fnm} : apc required : ${apc_real}"
