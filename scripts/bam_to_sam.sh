#!/bin/sh

for bam in ~/analyzed/ALS_RNASeq_call_variants/2-pass_output/*/*.bam; do
    sambamba view $bam -h -t 32 -o ~/analyzed/ALS_RNASeq_call_variants/2-pass_output/sams/${bam##*/}.sam
done
