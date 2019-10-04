sudo port install zlib
sudo port install bzip2
sudo port install lzma
sudo port install curl

git clone git@github.com:samtools/htslib.git

pushd htslib
autoreconf
./configure
make
make install prefix=~/Projects/Addgene/addgene-bioinformatics/local-tmp
popd
