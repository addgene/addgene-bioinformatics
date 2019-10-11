DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd | sed s%/src/sh%% )"

git clone https://github.com/lh3/bwa.git

pushd bwa
make
mv bwa $DIR/local-tmp/bin
popd

rm -rf bwa
