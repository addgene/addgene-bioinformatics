DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd | sed s%/src/sh%% )"

hg clone http://last.cbrc.jp/last/

pushd last
make
make install prefix=$DIR/local-tmp
popd

rm -rf last
