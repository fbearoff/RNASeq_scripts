#!/bin/bash

set -euo pipefail

sam_dir="/home/frank/analyzed/ALS_RNASeq_call_variants/2-pass_output/sams"
destination_dir="/home/frank/R_projects/ALS_RNASeq/DEXseq"
DEXseq_scripts_dir="/home/frank/R/x86_64-pc-linux-gnu-library/3.3/DEXSeq/python_scripts"

GTF="/home/frank/GRCm38_construct/Mus_musculus.GRCm38.88.gtf"
GFF="/home/frank/GRCm38_construct/Mus_musculus.GRCm38.88.gff"

#python "$DEXseq_scripts_dir"dexseq_prepare_annotation.py $GTF $GFF

mkdir -p $destination_dir || exit

for sample in $sam_dir/*.sam; do
    sample_short=${sample##*/}
    python "$DEXseq_scripts_dir"/dexseq_count.py -p yes -r pos $GFF $sample $destination_dir/${sample_short%%.*}.txt
done
