#!/bin/sh

#call SNPs using standard bcftools pipeline

REF="/home/frank/GRCm38_construct/Mus_musculus.GRCm38.dna.primary_assembly.fa"
DATA_DIR="/home/frank/data/Sanger_sequences"

#Get system information
THREADS="$(grep -c ^processor /proc/cpuinfo)"

#Call SNPs
for sample in "$DATA_DIR"/*.bam; do
    echo "Calling SNPs on $sample"
    bcftools mpileup -Ou -f $REF $sample | bcftools call --threads "$THREADS" -mv -Ob -o "${sample%.bam}"_calls.bcf
    bcftools index "${sample%.bam}"_calls.bcf
done
