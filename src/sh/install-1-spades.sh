DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd | sed s%/src/sh%% )"

wget http://cab.spbu.ru/files/release3.13.0/SPAdes-3.13.0.tar.gz

tar -xzf SPAdes-3.13.0.tar.gz

pushd SPAdes-3.13.0
PREFIX=$DIR/local-tmp ./spades_compile.sh
popd

rm SPAdes-3.13.0.tar.gz
rm -rf SPAdes-3.13.0
