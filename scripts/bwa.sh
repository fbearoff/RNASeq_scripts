#!/bin/sh

#inputs
REF="/home/frank/GRCm38_construct/Mus_musculus.GRCm38.dna.primary_assembly.fa"
read1="/home/frank/data/SJL_sequence_Cox_lab/GT17_10963_WGS_SJL_USPD16082730_H7TKJCCXY_L2_1.fq.gz"
read2="/home/frank/data/SJL_sequence_Cox_lab/GT17_10963_WGS_SJL_USPD16082730_H7TKJCCXY_L2_2.fq.gz"
OUTPUT_FOLDER="/scratch/frank/SJL_sequence_Cox"

#Abandon script if any command fails
set -euo pipefail

#get system information
THREADS="$(grep -c ^processor /proc/cpuinfo)"
START="$(date)"

#Check for proper variables
if [ ! -f "$REF" ]; then
    echo "Reference assembly file not present, exiting!" && exit
fi

#make scratch directory
mkdir -p $OUTPUT_FOLDER ||exit

#index reference genome
if [ ! -f "$REF".sa ]; then
    bwa index -a bwtsw $REF
fi

#map reads
bwa mem -t $THREADS $REF $read1 $read2 |samblaster |samtools view -@ "$THREADS" -u |samtools sort -@ "$THREADS" -l 9 -o "$OUTPUT_FOLDER"/SJL.sorted2.bam

#move output to final location
mkdir -p "$HOME/analyzed" || exit
echo "Placing files in $HOME/analyzed/"
mv $OUTPUT_FOLDER $HOME/analyzed

echo "Started at $START"
echo -n "Finished at " && date
exit
