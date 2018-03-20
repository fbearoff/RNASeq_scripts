#!/bin/bash

################################################################################
#Created by Frank Bearoff                                                      #
#Drexel University College of Medicine Department of Microbiology & Immunology #
################################################################################

read_files="/home/frank/data/combined_ALS_RNASeq"
destination_directory="/home/frank/R_projects/combined_ALS_RNASeqs"
read_type="ISR"
transcriptome="/home/frank/GRCm38_construct/Mus_musculus.GRCm38.cdna.all.fa.gz"

# Create index of transcriptome
#k=31 is recommended value for reads of 75bp or longer
salmon index -t "$transcriptome" -i "$destination_directory"/index --type quasi -k 31

#quantify transcripts against index
for sample in "$read_files"/*.R1.fq.gz;
do
base_name="${sample##*/}"
echo "Processing sample ${base_name%%.R1*}"
salmon quant -i "$destination_directory"/index/ --gcBias -l $read_type -1"${sample}" -2 "${sample%%.R1*}.R2.fq.gz" -o "$destination_directory"/quants/"${base_name%%.R1*}"
done
