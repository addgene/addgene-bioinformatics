#!/bin/bash

export PATH="../../local/bin:$PATH"

fastq_home="../../dat"
fastq_dirs=("A12417_sW0200_FASTQ"
            "representative-assortment"
            "simple-data"
            "challenge-data"
           )

fastq_home="../../dat/miscellaneous"
fastq_dirs=("A11967A_sW0154_FASTQ")

for fastq_dir in ${fastq_dirs}; do

    date > time_${fastq_dir}.out

    inp_size=$(du -m ${fastq_home}/${fastq_dir} | cut -f 1)
    echo
    echo "FASTQ directory : ${fastq_dir} : input size (MB) : ${inp_size}" \
        | tee -a time_${fastq_dir}.out

    fastq_fnms=$(ls -1 ${fastq_home}/${fastq_dir}/*_R1_001.fastq.gz)
    for fastq_pth in ${fastq_fnms}; do

        fastq_fnm=$(basename ${fastq_pth} _R1_001.fastq.gz)

        if [ -d ${fastq_fnm} ]; then
            continue

        fi
        
        echo
        spades_start=$(date +%s)

        spades.py \
            -1 ${fastq_home}/${fastq_dir}/${fastq_fnm}_R1_001.fastq.gz \
            -2 ${fastq_home}/${fastq_dir}/${fastq_fnm}_R2_001.fastq.gz \
            -o ${fastq_fnm} \
            --cov-cutoff 100 \
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

        apc_alns=$(ls apc_aln*)
        for apc_aln_old in $apc_alns; do
            apc_aln_new=$(echo $apc_aln_old | sed s/_aln/_${fastq_fnm}_aln/)
            mv $apc_aln_old $apc_aln_new

        done
 
     done

done

out_size=$(du -cm *A[0-9]*_[A-Z][0-9]* | grep 'total' | cut -f 1)
echo
echo "FASTQ directory : ${fastq_dir} : output size (MB) : ${out_size}" \
    | tee -a time_${fastq_dir}.out

grep "input size" time_${fastq_dir}.out \
    | cut -d : -f 4 \
          > time_${fastq_dir}_input_size.out

grep "SPAdes required" time_${fastq_dir}.out \
    | cut -d : -f 4 \
          > time_${fastq_dir}_spades.out

grep "apc required" time_${fastq_dir}.out \
    | cut -d : -f 4 \
    | cut -d m -f 2 \
    | sed s/s// \
          > time_${fastq_dir}_apc.out

grep "output size" time_${fastq_dir}.out \
    | cut -d : -f 4 \
          > time_${fastq_dir}_output_size.out
