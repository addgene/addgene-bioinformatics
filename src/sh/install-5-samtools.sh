sudo port install ncurses

git clone git@github.com:samtools/samtools.git

pushd samtools
autoreconf
./configure
make
make install prefix=~/Projects/Addgene/addgene-bioinformatics/local-tmp
popd

rm -rf samtools
rm -rf htslib
