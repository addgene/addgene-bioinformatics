DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd | sed s%/src/sh%% )"

pushd $DIR/local-tmp/bin
rm metaspades.py
rm plasmidspades.py
rm rnaspades.py
rm spades-bwa
rm spades-core
rm spades-corrector-core
rm spades-gbuilder
rm spades-gmapper
rm spades-hammer
rm spades-ionhammer
rm spades-kmercount
rm spades-truseq-scfcorrection
rm spades.py
rm spades_init.py
rm truspades.py
popd

pushd $DIR/local-tmp/share
rm -rf spades
popd
