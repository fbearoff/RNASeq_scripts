#!/bin/bash

sra_dir=$1

mkdir -p /scratch/frank/fastq/tmp

for sra in $sra_dir/*.sra; do
    base=${sra##*/}
    echo "Splitting sample ${base%.sra}"
    fasterq-dump $sra --split-files -O /scratch/frank/fastq -t /scratch/frank/fastq/tmp
done

for fq in /scratch/frank/fastq/*.fastq; do
    gzip $fq
done

