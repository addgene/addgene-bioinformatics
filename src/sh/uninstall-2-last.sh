DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd | sed s%/src/sh%% )"

pushd $DIR/local-tmp/bin
rm fastq-interleave
rm last-dotplot
rm last-map-probs
rm last-merge-batches
rm last-pair-probs
rm last-postmask
rm last-split
rm last-split8
rm last-train
rm lastal
rm lastal8
rm lastdb
rm lastdb8
rm maf-convert
rm maf-cut
rm maf-join
rm maf-sort
rm maf-swap
rm parallel-fasta
rm parallel-fastq
popd