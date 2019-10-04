hg clone http://last.cbrc.jp/last/

pushd last
make
make install prefix=~/Projects/Addgene/addgene-bioinformatics/local-tmp
popd

rm -rf last
