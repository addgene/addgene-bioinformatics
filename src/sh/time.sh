#!/bin/bash

export PATH="../../local/bin:$PATH"

fastq_home="../../dat"
fastq_dirs=("A12417_sW0200_FASTQ"
            "representative-assortment"
            "simple-data"
            "challenge-data"
           )

for fastq_dir in ${fastq_dirs}; do

    date > time_${fastq_dir}.out

    dir_size=$(du -m ${fastq_home}/${fastq_dir} | cut -f 1)
    echo
    echo "FASTQ directory : ${fastq_dir} : size (MB) : ${dir_size}" \
        | tee -a time_${fastq_dir}.out

    fastq_fnms=$(ls -1 ${fastq_home}/${fastq_dir}/*_R1_001.fastq.gz)
    for fastq_pth in ${fastq_fnms}; do

        fastq_fnm=$(basename ${fastq_pth} _R1_001.fastq.gz)

        echo
        spades_start=$(date +%s)
        spades.py \
            -1 ${fastq_home}/${fastq_dir}/${fastq_fnm}_R1_001.fastq.gz \
            -2 ${fastq_home}/${fastq_dir}/${fastq_fnm}_R2_001.fastq.gz \
            -o ${fastq_fnm} --cov-cutoff 100 \
            | tee spades_${fastq_fnm}.out \
            | grep "^=="
        spades_end=$(date +%s)
        spades_real=$(expr $spades_end - $spades_start)
        echo
        echo "FASTQ filename : ${fastq_fnm} : SPAdes required (s) : ${spades_real}" \
            | tee -a time_${fastq_dir}.out
        
        { time ../apc/apc.pl \
               -b ${fastq_fnm} ${fastq_fnm}/scaffolds.fasta \
               > apc_${fastq_fnm}.out; } > apc_${fastq_fnm}_time.out 2>&1
        apc_real=$(grep "real" apc_${fastq_fnm}_time.out | cut -f 2)
        echo
        echo "FASTQ filename : ${fastq_fnm} : apc required : ${apc_real}" \
            | tee -a time_${fastq_dir}.out

        break

    done

done
