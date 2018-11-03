#! /bin/sh

try () {
    echo TEST "$@"
    eval "$@"
    echo
}

cd $(dirname $0)

PATH=../scripts:$PATH

{
    try maf-swap -h
    try maf-swap bs100.maf
    try maf-swap -n1 90089.maf
    try maf-swap -n3 ../examples/multiMito.maf
} 2>&1 |
diff -u $(basename $0 .sh).out -
