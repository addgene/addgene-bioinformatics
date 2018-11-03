#! /bin/sh

# Exercise LAST programs, and compare the output to a reference
# output.  More tests should be added!

try () {
    echo TEST "$@"
    eval "$@"
    echo
}

cd $(dirname $0)

# Make sure we use this version of LAST:
PATH=../src:../scripts:$PATH

dnaSeq=galGal3-M-32.fa
protSeq=Q2LCP8.fa
fastq=SRR001981-1k.fastq
gc=../examples/vertebrateMito.gc
db=/tmp/last-test

trap 'rm -f $db*' EXIT

{
    # spaced seeds, soft-masking, centroid alignment, matrix file
    lastdb -c -m110 $db $dnaSeq
    try lastal -fMAF -u1 -j5 -p ../data/HOXD70.mat -z3400 -e2500 $db $dnaSeq

    # multiple volumes & query batches
    lastdb -m1 -s1 $db $dnaSeq
    try lastal -fTAB -i1 -w0 -e40 $db $dnaSeq

    # match-counting, with multiple query batches
    try lastal -j0 -i1 -s0 $db $dnaSeq

    # FASTQ quality scores
    try lastal -Q1 -e90 -a9 $db $fastq

    # gapless translated alignment & genetic code file
    lastdb -p $db $protSeq
    try lastal -F12 -pBL62 -e40 -G $gc -j1 $db $dnaSeq

    # subset seed file, soft-masking
    lastdb -c -u ../data/YASS.seed $db $dnaSeq
    try lastal -s0 -f0 -e18 $db $dnaSeq

    # asymmetric scoring matrix
    try lastal -s0 -f0 -p asymmetric.mat -e2000 $db $dnaSeq

    # FASTQ-Illumina quality scores
    lastdb -m1111110 $db $dnaSeq
    try lastal -Q3 -e110 $db illumina100.txt

    # PRB-format quality data
    try lastal -Q4 -e90 $db mouse_tss_prb.txt

    # probabilistic alignment with quality scores
    try lastal -Q1 -j6 -e90 -a9 $db $fastq

    # sparse index, generalized affine gap costs
    lastdb -w2 -c $db $dnaSeq
    try lastal -r3 -q3 -a21 -c2 -e60 -f0 $db $dnaSeq

    # generalized affine gaps, frameshifts, tabular output
    lastdb -p -c $db $protSeq
    try lastal -F12 -pBL62 -c2 -e40 -f0 $db $dnaSeq

    # gapless alignment, protein-protein alignment, seed freq
    try lastal -j1 -f0 -e37 -m100 $db $protSeq

    # fastq-versus-fastq, seed freq
    lastdb -Q1 $db sd-ccs-100.fq
    try lastal -Q1 -r1 -q2 -a1 -b1 -e44 -m100 -s0 $db sd-ccs-100.fq

    # incomplete sorting, lastal on one volume
    lastdb -i10 -s1 $db $dnaSeq
    try lastal -Q1 -e90 -a9 -f0 ${db}0 $fastq

    # multiple seeds, transition constraints
    lastdb -c -m 11101T011T11,111001010010111 $db $dnaSeq
    try lastal -s0 -f0 -e18 $db $dnaSeq

    # Iedera notation
    lastdb -c -m '#@#--##--#-#' $db $dnaSeq
    try lastal -s0 -f0 -e18 $db $dnaSeq

    # overlap alignment, tabular output ending in gaps
    lastdb -m1111110 $db $dnaSeq
    try lastal -T1 -Q1 -e60 -a9 -f0 $db $fastq

    # probabilistic overlap alignment
    try lastal -T1 -Q1 -e60 -a9 -j4 $db $fastq

    # fastq-versus-fasta gapless overlap alignment
    try lastal -T1 -Q1 -e60 -j1 -fTAB $db $fastq

    # expected counts
    try lastal -s0 -e18 -j7 $db $dnaSeq

    # named multi-seed, sparse query seeding
    lastdb -c -uMAM8 $db hg19-M.fa
    try lastal -e34 -k128 -f0 $db galGal3-M-32.fa

    # named score matrix, sparse query seeding
    try lastal -pHOXD70 -e4500 -k128 -f0 $db galGal3-M-32.fa

    # MAM4, gapless alignment culling
    lastdb -uMAM4 $db hg19-M.fa
    try lastal -e34 -C2 -f0 $db galGal3-M-32.fa

    # minimum seed length
    try lastal -e34 -f0 -l30 $db galGal3-M-32.fa

    # match-counting with min & max lengths
    lastdb -m1 $db $dnaSeq
    try lastal -j0 -l4 -L11 -s0 $db $dnaSeq

    lastdb -i10 $db tttttccccc.fa
    try lastal -e5 -f0 $db ttttt.fa | grep -v '^#'

    # tantan masking on DNA
    lastdb -cR01 $db galGal3-M-32.fa
    try lastal -e40 $db hg19-M.fa

    # tantan masking on protein
    lastdb -pcR01 $db Q2LCP8.fa
    try lastal -e100 $db Q5GS15.fa

    # tantan masking for translated alignment
    try lastal -F15 -pBLOSUM62 -e100 $db galGal3-M-32.fa

    # AT-rich DNA, tantan
    lastdb -cR02 $db at-rich.fa
    try lastal -pAT77 -e100 -s0 $db at-rich.fa

    # fastq + tantan
    lastdb -R01 $db $dnaSeq
    try lastal -Q1 -a15 -b3 -e80 $db nano.fq

    # fasta query versus fastq reference
    lastdb -Q1 $db sd-ccs-100.fq
    try lastal -a1 -D1000 $db galGal3-M-32.fa

    # prb query versus fastq reference
    try lastal -Q4 -a1 -D100 $db mouse_tss_prb.txt

    # fastq DNA versus protein
    lastdb -pcR00 $db Q2LCP8.fa
    try lastal -Q1 -pBL62 -F12 -D1000 $db sd-ccs-100.fq

    # BlastTab format
    try lastal -fBlastTab -pBL62 -b1 -F15 -D1e3 $db galGal3-M-32.fa

    # BlastTab+ format
    try lastal -fBlastTab+ -pBL62 -b1 -F15 -D1e3 $db galGal3-M-32.fa

    # DNA-versus-protein alignment without frameshifts
    try lastal -fTAB -pBL62 -b1 -F0 -D1e3 $db galGal3-M-32.fa

    # strand asymmetry
    lastdb $db hg19-M.fa
    try lastal -S1 -pBISF -Q1 -e120 -f0 -j4 $db bs100.fastq

    # culling
    try lastal -D1000 -fTAB -K2 $db galGal3-M-32.fa

    # minimizers
    lastdb -W3 $db galGal3-M-32.fa
    try lastal -W19 -fTAB $db hg19-M.fa

    # minimum-difference alignment
    try lastal -W1 -M -fTAB $db hg19-M.fa

    # asymmetric gap costs
    try lastal -fTAB -j4 -A2 -B2 $db hg19-M.fa
    try lastal -fTAB -j4 -Q1 -e90 -a7 -A12 -B4 $db $fastq

    # fastq-ignore
    try lastal -fTAB -j4 -Q0 -e90 -a7 -A12 -B4 -b9 -r6 -q18 $db $fastq

    # first alignments only
    try lastal -N2 $db hg19-M.fa

    # "U" nucleotide
    try "tr Tt Uu < hg19-M.fa | lastal -N2 $db"

    # -z %
    try lastal -W19 -z50% -fTAB $db hg19-M.fa

    # -x g
    try lastal -W19 -x1g -fTAB $db hg19-M.fa

    # non-negative score matrix
    try lastal -L9 -m0 -j1 -q0 -d5 -n6 -fTAB $db tttttccccc.fa

    # fastq-versus-fastq gapless overlap alignment
    lastdb -Q1 -uNEAR -cR01 $db $fastq
    try lastal -Q1 -T1 -j1 -s0 $db $fastq
} 2>&1 |
grep -v version |  # omit header lines with the LAST version number
diff -u last-test.out -

# Test: last-bisulfite, last-merge-batches, last-split, named seeds
lastdb -uBISF f hg19-M.fa
lastdb -uBISR r hg19-M.fa
../examples/last-bisulfite.sh f r bs100.fastq | grep -v '^#' | diff bs100.maf -
rm f.* r.*

./last-pair-test.sh
./last-postmask-test.sh
./last-split-test.sh
./last-train-test.sh
./maf-convert-test.sh
./maf-swap-test.sh

# Test: lastdb, lastal, last-split, maf-sort, maf-join
cd ../examples
./multiMito.sh | diff multiMito.maf -
