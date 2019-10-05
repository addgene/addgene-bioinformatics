DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd | sed s%/src/sh%% )"

pushd $DIR/local-tmp/bin
rm bgzip
rm htsfile
rm tabix
popd
  
pushd $DIR/local-tmp/include
rm -rf htslib
popd

pushd $DIR/local-tmp/lib
rm libhts.1.9-393-g634aad4-dirty.dylib
rm libhts.2to3part12.dylib
rm libhts.a
rm libhts.dylib
rm -rf pkgconfig
popd

pushd $DIR/local-tmp/share/man/man1
rm bgzip.1
rm htsfile.1
rm tabix.1
popd

pushd $DIR/local-tmp/share/man/man5
rm faidx.5
rm sam.5
rm vcf.5
popd
