DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd | sed s%/src/sh%% )"

sudo port install zlib
sudo port install bzip2
sudo port install lzma
sudo port install curl

git clone git@github.com:samtools/htslib.git

pushd htslib
autoreconf
./configure
make
make install prefix=$DIR/local-tmp
popd
