#!/bin/bash

bam_folder="/home/frank/analyzed/ALS_RNASeq2_call_variants/1-pass"
fastq_folder="/home/frank/data/ALS_RNASeq2/20180302blankenhorn"
samples_file="/home/frank/data/ALS_RNASeq2_samples.txt"
output_folder="/scratch/frank/analyzed/ALS_RNASeq2_call_variants/QC"
decoder_file="/home/frank/R_projects/ALS_RNASeq2/samples.txt"
fasta_ref="/home/frank/GRCm38_construct/Mus_musculus.GRCm38.dna.primary_assembly_space_stripped.fa"
chrom_sizes="/home/frank/GRCm38_construct/Mus_musculus.GRCm38.dna.primary_assembly.fa.chrom_sizes"
gtf_file="/home/frank/GRCm38_construct/Mus_musculus.GRCm38.88.gtf"
qorts_location="/home/frank/bin/QoRTs.jar"

#Read in samples
declare -A aa_samples
while IFS=: read -r genotype sample; do
    aa_samples[$genotype]+=$(echo "$sample ")
done < "$samples_file"

echo "Samples are:"
for genotype in "${!aa_samples[@]}"; do
    echo "$genotype: ${aa_samples[$genotype]}"
done

read -p "Are these samples correct [Y/n]? " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]
    then  [[ "$0" = "$BASH_SOURCE" ]] && echo "Exiting!" && exit 1 || return 1
fi

mkdir -p "$output_folder" || exit

#Do QC on samples
for sample in ${aa_samples[@]}; do
    echo "> Running QC on sample $sample"
    java -jar $qorts_location QC \
        --stranded \
        --chromSizes "$chrom_sizes" \
        --title "$sample" \
        --generatePlots \
        --genomeFA $fasta_ref \
        --rawfastq "$fastq_folder"/"$sample".R1.fq.gz,"$fastq_folder"/"$sample".R2.fq.gz \
        "$bam_folder"/"$sample"/Aligned.sortedByCoord.out.bam \
        $gtf_file \
        "$output_folder"/"$sample"
done

#Generate novel junctions GFF for JunctionSeq
java -jar $qorts_location \
        mergeNovelSplices \
        --minCount 6 \
        --stranded \
        "$output_folder" \
        "$decoder_file" \
        "$gtf_file" \
        "$output_folder"

#generate summary wiggle tracks
mkdir -p "$output_folder"/mergedTracks
for genotype in "${!aa_samples[@]}"; do
    echo "Merging wiggle tracks for $genotype samples"
    java -jar "$qorts_location" \
        mergeWig \
        --calcMean \
        --infilePrefix "$output_folder"/ \
        --infileSuffix /QC.wiggle.fwd.wig.gz \
        --sampleList $(echo "${aa_samples[$genotype]}" |sed 's/\ /,/g' |sed 's/,*$//g') \
        --sizeFactorFile "$output_folder"/JS.GEO.size.factors.txt \
        --trackTitle "$genotype"_FWD \
        "$output_folder"/mergedTracks/"$genotype".fwd.wig.gz
    java -jar "$qorts_location" \
        mergeWig \
        --calcMean \
        --infilePrefix "$output_folder"/ \
        --infileSuffix /QC.wiggle.rev.wig.gz \
        --sampleList $(echo "${aa_samples[$genotype]}" |sed 's/\ /,/g' |sed 's/,*$//g') \
        --sizeFactorFile "$output_folder"/JS.GEO.size.factors.txt \
        --trackTitle "$genotype"_REV \
        "$output_folder"/mergedTracks/"$genotype".rev.wig.gz
done
