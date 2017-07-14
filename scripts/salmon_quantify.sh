#!/bin/bash

################################################################################
#Created by Frank Bearoff                                                      #
#Drexel University College of Medicine Department of Microbiology & Immunology #
################################################################################

read_files="/home/frank/data/first_rnaseq/6LibbysRNAs"
destination_directory="/home/frank/R_projects/rnaseq2"
read_type="SR"
transcriptome="/home/frank/GRCm38_construct/Mus_musculus.GRCm38.cdna.all.fa.gz"

# Create index of transcriptome
#k=31 is recommended value for reads of 75bp or longer
salmon index -t "$transcriptome" -i "$destination_directory"/index --type quasi -k 31

#quantify transcripts against index
#This script assumes single-end reads with -r flag
for sample in "$read_files"/*.gz;
do
echo "Processing sample ${sample#*/}"
salmon quant -i "$destination_directory"/index/ -l "$read_type" -r "${sample}" -o "$destination_directory"/quants/"${sample#*/}"_quant
done
