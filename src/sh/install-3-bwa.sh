git clone https://github.com/lh3/bwa.git

pushd bwa
make
mv bwa ~/Projects/Addgene/addgene-bioinformatics/local-tmp/bin
popd

rm -rf bwa
