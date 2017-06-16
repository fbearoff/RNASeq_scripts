#!/bin/bash

#Created by Frank Bearoff fb99@drexel.edu
#Based on GATK pipeline https://software.broadinstitute.org/gatk/documentation/article.php?id=3891

#Dependencies (known working version)
    #GATK picard (2.9.2)
    #GATK (3.7-0)
    #STAR (STAR_2.5.3a)
    #seqtk (1.2-r101-dirty)
    #gnu-parallel (20170422)
    #samtools (1.4.1-1-g67999ae)
    #sambamba (0.6.6)
    #samblaster (0.1.24)
    #R (3.3.3)
    #igvtools (2.3.93)

#Notes
    #Designed for single end reverse strand reads
    #KNOWN_SITES needs to be matched to reference.dict
        #java -jar picard.jar UpdateSequenceDictionary \
        #I=original.vcf \
        #O=matched.vcf \
        #SD=reference.dict
    #KNOWN_SITES file needs to be sorted in same order as reference
        #java -jar picard.jar SortVcf \
        #I=original.vcf \
        #O=sorted.vcf \
        #SD=reference.dict
    #Directories cannot have a trailing slash for config
    #Supply SAMPLES_FILE as genotype:sample per line (full read file names e.g. *.R1.fq.gz)

##########CONFIG############
DATA_FOLDER="/home/frank/data/combined_rnaseq_1and2"  #location of read files
SAMPLES_FILE="/home/frank/data/combined_rnaseq_1and2/samples.txt"
OUTPUT_FOLDER="/scratch/frank/combined_rnaseq_1and2"
REF_LOCATION="/home/frank/GRCm38_construct/Mus_musculus.GRCm38.dna.primary_assembly.fa"    #Ensembl .fa
REF_ANNOTATION="/home/frank/GRCm38_construct/Mus_musculus.GRCm38.88.gtf"  #Ensembl .gtf
KNOWN_SITES="/home/frank/GRCm38_construct/Mus_musculus_all_chr_sorted_matched.vcf"      #ftp://ftp.ncbi.nih.gov/snp/organisms/mouse_10090/VCF/
PICARD_LOCATION="/home/frank/bin/picard.jar"
GENOMEANALYSISTK_LOCATION="/home/frank/bin/GenomeAnalysisTK.jar"
IGVTOOLS_LOCATION="/home/frank/tools//IGVTools/igvtools.jar"
##########CONFIG############

#Abandon script if any command fails
set -euo pipefail

#Gather system information
THREADS="$(grep -c ^processor /proc/cpuinfo)"
NJOBS="$(awk '( $1 == "MemTotal:" ) { printf "%i\n", $2/1048576/15 }' /proc/meminfo)" #Allocate 15G per job
START="$(date)"

#Check for proper variables
if [ ! -d "$DATA_FOLDER" ]; then
    echo "Data folder not present, exiting!" && exit
elif [ ! -f "$REF_LOCATION" ]; then
    echo "Reference assembly file not present, exiting!" && exit
elif [ ! -f "$REF_ANNOTATION" ]; then
    echo "Annotation file not present, exiting!" && exit
elif [ ! -f "$KNOWN_SITES" ]; then
    echo "Known sites file not present, exiting!" && exit
elif [ ! -f "$PICARD_LOCATION" ]; then
    echo "Picard JAR binary not present, exiting!" && exit
elif [ ! -f "$GENOMEANALYSISTK_LOCATION" ]; then
    echo "GATK JAR binary not present, exiting!" && exit
fi

#Read in samples
declare -A aa_samples
while IFS=: read -r key value; do
    aa_samples[$key]+=$(echo "$value ")
done < "$SAMPLES_FILE"

echo "Samples are:"
for key in "${!aa_samples[@]}"; do
    echo "$key: ${aa_samples[$key]}"
done

read -p "Are these samples correct [Y/n]? " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]
    then  [[ "$0" = "$BASH_SOURCE" ]] && echo "Exiting!" && exit 1 || return 1
fi

mkdir -p "$OUTPUT_FOLDER" || exit

#Prepare Reference
echo ">Preparing reference genome"
if [ ! -f "${REF_LOCATION%.fa}".dict ]; then
    java -jar "$PICARD_LOCATION" CreateSequenceDictionary R= "$REF_LOCATION" O= "${REF_LOCATION%.fa}".dict
elif [ ! -f "$REF_LOCATION".fai ]; then
    samtools faidx "$REF_LOCATION"
fi

#Reverse complement FASTA/Q:
mkdir "$OUTPUT_FOLDER"/reverse_complement || exit
cd "$DATA_FOLDER" || exit
echo ">Reverse complementing reads to $OUTPUT_FOLDER/reverse_complement"
parallel -j "$NJOBS" 'seqtk seq -r {} |gzip > '"$OUTPUT_FOLDER"'/reverse_complement/{}' ::: *.gz

#STAR generate index
echo ">Generating genome index"
mkdir "$OUTPUT_FOLDER"/genome || exit
STAR --runThreadN "$THREADS" --runMode genomeGenerate --genomeDir "$OUTPUT_FOLDER"/genome --genomeFastaFiles "$REF_LOCATION" --sjdbGTFfile "$REF_ANNOTATION" --outFileNamePrefix "$OUTPUT_FOLDER"/genome/

#1-Pass map reads against genome index
mkdir "$OUTPUT_FOLDER"/1-pass/ || exit
for sample in "$OUTPUT_FOLDER"/reverse_complement/*.gz; do
    echo ">Aligning sample ${sample##*/}"
    mkdir "$OUTPUT_FOLDER"/1-pass/"${sample##*/}"
    STAR --runThreadN "$THREADS" --genomeDir "$OUTPUT_FOLDER"/genome --readFilesIn "${sample}" --readFilesCommand gunzip -c --outFileNamePrefix "$OUTPUT_FOLDER"/1-pass/"${sample##*/}"/
done

#2-Pass map reads against genome index, mark duplicates, sort, and compress
mkdir "$OUTPUT_FOLDER"/2-pass_output || exit

for genotype in "${!aa_samples[@]}"; do
    cat $(echo "$OUTPUT_FOLDER/1-pass/$genotype*/SJ.out.tab") > "$OUTPUT_FOLDER/1-pass/$genotype.SJ.out.tab"
    for sample in $(echo "${aa_samples[$genotype]}"); do
        mkdir "$OUTPUT_FOLDER"/2-pass_output/"$sample"
        STAR --runThreadN "$THREADS" --genomeDir "$OUTPUT_FOLDER"/genome --readFilesIn "$OUTPUT_FOLDER"/reverse_complement/"$sample" --readFilesCommand gunzip -c --outFileNamePrefix "$OUTPUT_FOLDER"/2-pass_output/"$sample"/ --outSAMattrRGline ID:"${sample%%.*}" SM:"${sample%%.*}" LB:"${sample%%.*}" PL:illumina PU:nextseq --sjdbFileChrStartEnd "$OUTPUT_FOLDER/1-pass/$genotype.SJ.out.tab" --outStd SAM |samblaster|sambamba view -t "$THREADS" -S /dev/stdin -f bam -l 0 |sambamba sort -t "$THREADS" -l 9 /dev/stdin -o "$OUTPUT_FOLDER"/2-pass_output/"$sample"/"$sample".sorted.bam
    done
done

#Reassign MAPQ Scores and Split BAM
mkdir -p "$OUTPUT_FOLDER"/vc_BAMs/split || exit
cd "$OUTPUT_FOLDER"/reverse_complement || exit
echo ">Recalibrating and splitting BAMs for all samples"
parallel -j "$NJOBS" 'java -jar '"$GENOMEANALYSISTK_LOCATION"' -T SplitNCigarReads -R '"$REF_LOCATION"' -I '"$OUTPUT_FOLDER"'/2-pass_output/{}/{}.sorted.bam -o '"$OUTPUT_FOLDER"'/vc_BAMs/split/{}.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS' ::: *.gz

#Indel Realignment Target Creation
#http://gatkforums.broadinstitute.org/gatk/discussion/7156/howto-perform-local-realignment-around-indels
mkdir -p "$OUTPUT_FOLDER"/vc_BAMs/realign || exit
for sample in "$OUTPUT_FOLDER"/reverse_complement/*.gz; do
    echo ">Creating Realignment Targets for ${sample##*/}"
    java -jar "$GENOMEANALYSISTK_LOCATION" -nt "$THREADS" -T RealignerTargetCreator -R "$REF_LOCATION" -known "$KNOWN_SITES" -I "$OUTPUT_FOLDER"/vc_BAMs/split/"${sample##*/}".split.bam -o "$OUTPUT_FOLDER"/vc_BAMs/realign/"${sample##*/}".intervals
#Convert to BED for IGV
    awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' "$OUTPUT_FOLDER"/vc_BAMs/realign/"${sample##*/}".intervals > "$OUTPUT_FOLDER"/vc_BAMs/realign/"${sample##*/}".bed
    java -jar "$IGVTOOLS_LOCATION" index "$OUTPUT_FOLDER"/vc_BAMs/realign/"${sample##*/}".bed
done

#Indel Realignment
cd "$OUTPUT_FOLDER"/reverse_complement || exit
echo ">Performing INDEL realignment for all samples"
parallel -j "$NJOBS" 'java -Djava.io.tmpdir=/tmp -jar '"$GENOMEANALYSISTK_LOCATION"' -T IndelRealigner -R '"$REF_LOCATION"' -targetIntervals '"$OUTPUT_FOLDER"'/vc_BAMs/realign/{}.intervals -known '"$KNOWN_SITES"' -I '"$OUTPUT_FOLDER"'/vc_BAMs/split/{}.split.bam -o '"$OUTPUT_FOLDER"'/vc_BAMs/realign/{}.indelrealigner.bam' ::: *.gz

#Base Recalibration
#http://gatkforums.broadinstitute.org/gatk/discussion/2801/howto-recalibrate-base-quality-scores-run-bqsr
mkdir "$OUTPUT_FOLDER"/base_recalibration || exit
for sample in "$OUTPUT_FOLDER"/reverse_complement/*.gz; do
    echo ">Base Recalibration for sample ${sample##*/}"
    java -jar "$GENOMEANALYSISTK_LOCATION" -nct "$THREADS" -T BaseRecalibrator -R "$REF_LOCATION" -I "$OUTPUT_FOLDER"/vc_BAMs/realign/"${sample##*/}".indelrealigner.bam -knownSites "$KNOWN_SITES" -o "$OUTPUT_FOLDER"/base_recalibration/"${sample##*/}".recal_data.table
    java -jar "$GENOMEANALYSISTK_LOCATION" -nct "$THREADS" -T BaseRecalibrator -R "$REF_LOCATION" -I "$OUTPUT_FOLDER"/vc_BAMs/realign/"${sample##*/}".indelrealigner.bam -knownSites "$KNOWN_SITES" -BQSR "$OUTPUT_FOLDER"/base_recalibration/"${sample##*/}".recal_data.table -o "$OUTPUT_FOLDER"/base_recalibration/"${sample##*/}".post_recal_data.table
    java -jar "$GENOMEANALYSISTK_LOCATION" -T AnalyzeCovariates -R "$REF_LOCATION" -before "$OUTPUT_FOLDER"/base_recalibration/"${sample##*/}".recal_data.table -after "$OUTPUT_FOLDER"/base_recalibration/"${sample##*/}".post_recal_data.table -plots "$OUTPUT_FOLDER"/base_recalibration/"${sample##*/}".recalibration_plots.pdf
    java -jar "$GENOMEANALYSISTK_LOCATION" -nct "$THREADS" -T PrintReads -R "$REF_LOCATION" -I "$OUTPUT_FOLDER"/vc_BAMs/realign/"${sample##*/}".indelrealigner.bam -BQSR "$OUTPUT_FOLDER"/base_recalibration/"${sample##*/}".recal_data.table -o "$OUTPUT_FOLDER"/base_recalibration/"${sample##*/}".bam --bam_compression 9
done

#Call variants GVCF Pipeline
mkdir "$OUTPUT_FOLDER"/called_variants || exit
cd "$OUTPUT_FOLDER"/reverse_complement || exit
echo ">Calling variants"
parallel -j "$NJOBS" 'java -jar '"$GENOMEANALYSISTK_LOCATION"' -T HaplotypeCaller -R '"$REF_LOCATION"' -I '"$OUTPUT_FOLDER"'/base_recalibration/{}.bam --dontUseSoftClippedBases --dbsnp '"$KNOWN_SITES"' -ERC GVCF -stand_call_conf 20 -o '"$OUTPUT_FOLDER"'/called_variants/{}.g.vcf' ::: *.gz

#Genotype GVCF
for genotype in "${!aa_samples[@]}"; do
    base_command="java -jar $GENOMEANALYSISTK_LOCATION -T GenotypeGVCFs -R $REF_LOCATION --dbsnp $KNOWN_SITES -o $OUTPUT_FOLDER/called_variants/$genotype.vcf"
    for sample in $(echo "${aa_samples[$genotype]}"); do
        base_command=$base_command" -V $OUTPUT_FOLDER/called_variants/$sample.g.vcf"
    done
    eval "$base_command"
done

#Variant filtering
for genotype in "${!aa_samples[@]}"; do
    echo ">Filtering variants for ${genotype##*/}"
    java -jar "$GENOMEANALYSISTK_LOCATION" -T VariantFiltration -R "$REF_LOCATION" -V "$OUTPUT_FOLDER"/called_variants/"${genotype##*/}".vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o "$OUTPUT_FOLDER"/called_variants/"${genotype##*/}".filtered.vcf
done

#Move output to final destination
mkdir -p "$HOME"/analyzed/ || exit
echo ">Placing files in $HOME/analyzed/"
mv "$OUTPUT_FOLDER" "$HOME"/analyzed/

echo "Started at $START"
echo -n "Finished at " && date
exit
