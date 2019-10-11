DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

pushd $DIR
src/sh/install-1-spades.sh
src/sh/install-2-last.sh
src/sh/install-3-bwa.sh
src/sh/install-4-htslib.sh
src/sh/install-5-samtools.sh
src/sh/install-6-apc.sh
popd

echo
echo 'After installation:'
echo '  mv local-tmp local'
echo '  export PATH="<repository-home>/addgene-bioinformatics/local/bin:$PATH"'
echo
