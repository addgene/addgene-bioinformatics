#! /bin/sh

# Regression tests

cd $(dirname $0)

PATH=../scripts:$PATH

r=maf-convert

maf1=SRR359290-1k.maf
maf2=bs100.maf

{
    $r -h
    head -n999 $maf1 | $r axt
    head -n999 $maf1 | $r blast
    head -n999 $maf1 | $r -l100 blast
    $r blast $maf2
    $r blasttab $maf2
    head -n999 $maf1 | $r chain
    $r html -l100 $maf2
    head -n999 $maf1 | $r -n html
    head -n999 $maf1 | $r psl
    head -n999 $maf1 | $r -p psl
    $r psl $maf2
    $r -j1e9 psl 102.maf
    $r psl 90089.maf
    $r -n sam $maf2
    head -n999 $maf1 | $r -r 'ID:1 PL:ILLUMINA SM:x' sam
    $r -d sam $maf1
    head -n999 $maf1 | $r -n tab
    head -n999 $maf1 | $r tab
} |
diff maf-convert-test.out -