#!/bin/sh

#need bcf against reference genome
REF=$1
BCF=$2
OUTPUT=$3

cat $REF | bcftools consensus $BCF > "$OUTPUT".fa
