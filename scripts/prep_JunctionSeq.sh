#!/bin/bash

bam_folder="/scratch/frank/ALS_RNASeq2_call_variants/1-pass"
fastq_folder="/home/frank/data/ALS_RNASeq2/20180302blankenhorn"
samples_file="/home/frank/data/ALS_RNASeq2_samples.txt"
output_folder="/home/frank/analyzed/ALS_RNASeq2_call_variants/QC"
decoder_file="/home/frank/R_projects/ALS_RNASeq2/samples.txt"
fasta_ref="/home/frank/GRCm38_construct/Mus_musculus.GRCm38.dna.primary_assembly_space_stripped.fa"
gtf_file="/home/frank/GRCm38_construct/Mus_musculus.GRCm38.88.gtf"
qorts_location="/home/frank/bin/QoRTs.jar"

#Read in samples
declare -A aa_samples
while IFS=: read -r key value; do
    aa_samples[$key]+=$(echo "$value ")
done < "$samples_file"

echo "Samples are:"
for key in "${!aa_samples[@]}"; do
    echo "$key: ${aa_samples[$key]}"
done

read -p "Are these samples correct [Y/n]? " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]
    then  [[ "$0" = "$BASH_SOURCE" ]] && echo "Exiting!" && exit 1 || return 1
fi

mkdir -p "$output_folder" || exit

# for sample in ${aa_samples[@]}; do
#     echo "> Running QC on sample $sample"
#     java -jar $qorts_location QC \
#         --stranded \
#         --generatePlots \
#         --genomeFA $fasta_ref \
#         --rawfastq "$fastq_folder"/"$sample".R1.fq.gz,"$fastq_folder"/"$sample".R2.fq.gz \
#         "$bam_folder"/"$sample"/Aligned.sortedByCoord.out.bam \
#         $gtf_file \
#         "$output_folder"/"$sample"
# done
java -jar $qorts_location \
        mergeNovelSplices \
        --minCount 6 \
        --stranded \
        "$output_folder" \
        "$decoder_file" \
        "$gtf_file" \
        "$output_folder"
