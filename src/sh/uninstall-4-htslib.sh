pushd /Users/raymondleclair/Projects/Addgene/addgene-bioinformatics/local-tmp/bin
rm bgzip
rm htsfile
rm tabix
popd
  
pushd /Users/raymondleclair/Projects/Addgene/addgene-bioinformatics/local-tmp/include
rm -rf htslib
popd

pushd /Users/raymondleclair/Projects/Addgene/addgene-bioinformatics/local-tmp/lib
rm libhts.1.9-393-g634aad4-dirty.dylib
rm libhts.2to3part12.dylib
rm libhts.a
rm libhts.dylib
rm -rf pkgconfig
popd

pushd /Users/raymondleclair/Projects/Addgene/addgene-bioinformatics/local-tmp/share/man/man1
rm bgzip.1
rm htsfile.1
rm tabix.1
popd

pushd /Users/raymondleclair/Projects/Addgene/addgene-bioinformatics/local-tmp/share/man/man5
rm faidx.5
rm sam.5
rm vcf.5
popd
