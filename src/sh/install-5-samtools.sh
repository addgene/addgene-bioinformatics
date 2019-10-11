DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd | sed s%/src/sh%% )"

sudo port install ncurses

git clone git@github.com:samtools/samtools.git

pushd samtools
autoreconf
./configure
make
make install prefix=$DIR/local-tmp
popd

rm -rf samtools
rm -rf htslib
