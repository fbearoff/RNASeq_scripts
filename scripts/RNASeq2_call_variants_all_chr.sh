#!/bin/bash

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
    
#Notes
    #Designed for single end reverse strand reads
    #KNOWN_SITES needs to be matched to reference.dict
        #java -jar picard.jar UpdateSequenceDictionary
        #I=original.vcf
        #O=matched.vcf
        #SD=reference.dict
    #KNOWN_SITES file needs to be sorted in same order as reference
        #java -jar picard.jar SortVcf \
        #I=original.vcf \
        #O=sorted.vcf \
        #SD=reference.dict 
    #Directories cannot have a trailing slash for config
    #For second pass mapping, need to adjust sjdbFileChrStartEnd for sample names per genotype

##########CONFIG############
DATA_FOLDER="/home/frank/data/second_rnaseq/9FranksRNAs"  #location of read files
OUTPUT_FOLDER="/scratch/frank/RNASeq2_variants_all_chr"
REF_LOCATION="/home/frank/GRCm38_construct/Mus_musculus.GRCm38.dna.primary_assembly.fa"    #Ensembl .fa
REF_ANNOTATION="/home/frank/GRCm38_construct/Mus_musculus.GRCm38.88.gtf"  #Ensembl .gtf
KNOWN_SITES="/home/frank/GRCm38_construct/Mus_musculus_all_chr_sorted_matched.vcf"      #ftp://ftp.ncbi.nih.gov/snp/organisms/mouse_10090/VCF/
PICARD_LOCATION="/home/frank/bin/picard.jar"
GENOMEANALYSISTK_LOCATION="/home/frank/bin/GenomeAnalysisTK.jar"
THREADS="32"
##########CONFIG############

#Get start time
START="$(date)"

#Abandon script if any command fails
set -euo pipefail

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
parallel 'seqtk seq -r {} |gzip > '"$OUTPUT_FOLDER"'/reverse_complement/{}' ::: *.gz

#STAR generate index
echo ">Generating genome index"
mkdir "$OUTPUT_FOLDER"/genome || exit
STAR --runThreadN "$THREADS" --runMode genomeGenerate --genomeDir "$OUTPUT_FOLDER"/genome --genomeFastaFiles "$REF_LOCATION" --sjdbGTFfile "$REF_ANNOTATION" --outFileNamePrefix "$OUTPUT_FOLDER"/genome/

#1-Pass map reads against genome index
mkdir "$OUTPUT_FOLDER"/1-pass/ || exit
for sample in "$OUTPUT_FOLDER"/reverse_complement/*.gz;
do
echo ">Aligning sample ${sample##*/}"
mkdir "$OUTPUT_FOLDER"/1-pass/"${sample##*/}"
STAR --runThreadN "$THREADS" --genomeDir "$OUTPUT_FOLDER"/genome --readFilesIn "${sample}" --readFilesCommand gunzip -c --outFileNamePrefix "$OUTPUT_FOLDER"/1-pass/"${sample##*/}"/
done

#2-Pass map reads against genome index, mark duplicates, sort, and compress
mkdir "$OUTPUT_FOLDER"/2-pass_output || exit

#10.3 samples  #CHANGE FOR SAMPLE PREFIX IN COMMANDS
for sample in "$OUTPUT_FOLDER"/reverse_complement/103*.gz;
do
echo ">2-Pass aligning sample ${sample##*/}"
s="${sample##*/}"
mkdir "$OUTPUT_FOLDER"/2-pass_output/"${sample##*/}"
STAR --runThreadN "$THREADS" --genomeDir "$OUTPUT_FOLDER"/genome --readFilesIn "${sample}" --readFilesCommand gunzip -c --outFileNamePrefix "$OUTPUT_FOLDER"/2-pass_output/"${sample##*/}"/ --outSAMattrRGline ID:"${s%%.*}" SM:"${s%%.*}" LB:"${s%%.*}" PL:illumina PU:nextseq --sjdbFileChrStartEnd "$OUTPUT_FOLDER"/1-pass/1034.R1.fq.gz/SJ.out.tab "$OUTPUT_FOLDER"/1-pass/1035.R1.fq.gz/SJ.out.tab "$OUTPUT_FOLDER"/1-pass/1036.R1.fq.gz/SJ.out.tab --outStd SAM |samblaster|sambamba view -t "$THREADS" -S /dev/stdin -f bam -l 0 |sambamba sort -t "$THREADS" -l 9 /dev/stdin -o "$OUTPUT_FOLDER"/2-pass_output/"${sample##*/}"/"${sample##*/}".sorted.bam
done

#B6 samples    #CHANGE FOR SAMPLE PREFIX IN COMMANDS
for sample in "$OUTPUT_FOLDER"/reverse_complement/B6*.gz;
do
echo ">2-Pass aligning sample ${sample##*/}"
s="${sample##*/}"
mkdir "$OUTPUT_FOLDER"/2-pass_output/"${sample##*/}" 
STAR --runThreadN "$THREADS" --genomeDir "$OUTPUT_FOLDER"/genome --readFilesIn "${sample}" --readFilesCommand gunzip -c --outFileNamePrefix "$OUTPUT_FOLDER"/2-pass_output/"${sample##*/}"/ --outSAMattrRGline ID:"${s%%.*}" SM:"${s%%.*}" LB:"${s%%.*}" PL:illumina PU:nextseq --sjdbFileChrStartEnd "$OUTPUT_FOLDER"/1-pass/B64.R1.fq.gz/SJ.out.tab "$OUTPUT_FOLDER"/1-pass/B65.R1.fq.gz/SJ.out.tab "$OUTPUT_FOLDER"/1-pass/B66.R1.fq.gz/SJ.out.tab --outStd SAM |samblaster|sambamba view -t "$THREADS" -S /dev/stdin -f bam -l 0 |sambamba sort -t "$THREADS" -l 9 /dev/stdin -o "$OUTPUT_FOLDER"/2-pass_output/"${sample##*/}"/"${sample##*/}".sorted.bam
done

#FL10 samples    #CHANGE FOR SAMPLE PREFIX IN COMMANDS
for sample in "$OUTPUT_FOLDER"/reverse_complement/FL10*.gz;
do
echo ">2-Pass aligning sample ${sample##*/}"
s="${sample##*/}"
mkdir "$OUTPUT_FOLDER"/2-pass_output/"${sample##*/}" 
STAR --runThreadN "$THREADS" --genomeDir "$OUTPUT_FOLDER"/genome --readFilesIn "${sample}" --readFilesCommand gunzip -c --outFileNamePrefix "$OUTPUT_FOLDER"/2-pass_output/"${sample##*/}"/ --outSAMattrRGline ID:"${s%%.*}" SM:"${s%%.*}" LB:"${s%%.*}" PL:illumina PU:nextseq --sjdbFileChrStartEnd "$OUTPUT_FOLDER"/1-pass/FL101.R1.fq.gz/SJ.out.tab "$OUTPUT_FOLDER"/1-pass/FL102.R1.fq.gz/SJ.out.tab "$OUTPUT_FOLDER"/1-pass/FL103.R1.fq.gz/SJ.out.tab --outStd SAM |samblaster|sambamba view -t "$THREADS" -S /dev/stdin -f bam -l 0 |sambamba sort -t "$THREADS" -l 9 /dev/stdin -o "$OUTPUT_FOLDER"/2-pass_output/"${sample##*/}"/"${sample##*/}".sorted.bam
done
 
#Reassign MAPQ Scores and Split BAM
mkdir -p "$OUTPUT_FOLDER"/vc_BAMs/split || exit
cd "$OUTPUT_FOLDER"/reverse_complement || exit 
echo ">Recalibrating and splitting BAMs for all samples"
parallel 'java -jar '"$GENOMEANALYSISTK_LOCATION"' -T SplitNCigarReads -R '"$REF_LOCATION"' -I '"$OUTPUT_FOLDER"'/2-pass_output/{}/{}.sorted.bam -o '"$OUTPUT_FOLDER"'/vc_BAMs/split/{}.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS' ::: *.gz

#Indel Realignment Target Creation
#http://gatkforums.broadinstitute.org/gatk/discussion/7156/howto-perform-local-realignment-around-indels
mkdir "$OUTPUT_FOLDER"/vc_BAMs/realign || exit
for sample in "$OUTPUT_FOLDER"/reverse_complement/*.gz;
do
echo ">Creating Realignment Targets for ${sample##*/}"
java -jar "$GENOMEANALYSISTK_LOCATION" -nt "$THREADS" -T RealignerTargetCreator -R "$REF_LOCATION" -known "$KNOWN_SITES" -I "$OUTPUT_FOLDER"/vc_BAMs/split/"${sample##*/}".split.bam -o "$OUTPUT_FOLDER"/vc_BAMs/realign/"${sample##*/}".intervals
#Convert to BED for IGV
awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' "$OUTPUT_FOLDER"/vc_BAMs/realign/"${sample##*/}".intervals > "$OUTPUT_FOLDER"/vc_BAMs/realign/"${sample##*/}".bed
done

#Indel Realignment
cd "$OUTPUT_FOLDER"/reverse_complement || exit
echo ">Performing INDEL realignment for all samples"
parallel 'java -Djava.io.tmpdir=/tmp -jar '"$GENOMEANALYSISTK_LOCATION"' -T IndelRealigner -R '"$REF_LOCATION"' -targetIntervals '"$OUTPUT_FOLDER"'/vc_BAMs/realign/{}.intervals -known '"$KNOWN_SITES"' -I '"$OUTPUT_FOLDER"'/vc_BAMs/split/{}.split.bam -o '"$OUTPUT_FOLDER"'/vc_BAMs/realign/{}.indelrealigner.bam' ::: *.gz

#Base Recalibration
#http://gatkforums.broadinstitute.org/gatk/discussion/2801/howto-recalibrate-base-quality-scores-run-bqsr
mkdir "$OUTPUT_FOLDER"/base_recalibration || exit
for sample in "$OUTPUT_FOLDER"/reverse_complement/*.gz;
do
echo ">Base Recalibration for sample ${sample##*/}"
java -jar "$GENOMEANALYSISTK_LOCATION" -nct "$THREADS" -T BaseRecalibrator -R "$REF_LOCATION" -I "$OUTPUT_FOLDER"/vc_BAMs/realign/"${sample##*/}".indelrealigner.bam -knownSites "$KNOWN_SITES" -o "$OUTPUT_FOLDER"/base_recalibration/"${sample##*/}".recal_data.table
java -jar "$GENOMEANALYSISTK_LOCATION" -nct "$THREADS" -T BaseRecalibrator -R "$REF_LOCATION" -I "$OUTPUT_FOLDER"/vc_BAMs/realign/"${sample##*/}".indelrealigner.bam -knownSites "$KNOWN_SITES" -BQSR "$OUTPUT_FOLDER"/base_recalibration/"${sample##*/}".recal_data.table -o "$OUTPUT_FOLDER"/base_recalibration/"${sample##*/}".post_recal_data.table
java -jar "$GENOMEANALYSISTK_LOCATION" -T AnalyzeCovariates -R "$REF_LOCATION" -before "$OUTPUT_FOLDER"/base_recalibration/"${sample##*/}".recal_data.table -after "$OUTPUT_FOLDER"/base_recalibration/"${sample##*/}".post_recal_data.table -plots "$OUTPUT_FOLDER"/base_recalibration/"${sample##*/}".recalibration_plots.pdf
java -jar "$GENOMEANALYSISTK_LOCATION" -T PrintReads -R "$REF_LOCATION" -I "$OUTPUT_FOLDER"/vc_BAMs/realign/"${sample##*/}".indelrealigner.bam -BQSR "$OUTPUT_FOLDER"/base_recalibration/"${sample##*/}".recal_data.table -o "$OUTPUT_FOLDER"/base_recalibration/"${sample##*/}".bam
done

#Call variants
mkdir "$OUTPUT_FOLDER"/called_variants || exit
cd "$OUTPUT_FOLDER"/reverse_complement || exit
echo ">Calling variants"
parallel 'java -jar '"$GENOMEANALYSISTK_LOCATION"' -T HaplotypeCaller -R '"$REF_LOCATION"' -I '"$OUTPUT_FOLDER"'/base_recalibration/{}.bam --dbsnp '"$KNOWN_SITES"' -stand_call_conf 20 -o '"$OUTPUT_FOLDER"'/called_variants/{}.vcf' ::: *.gz

#Variant filtering
for sample in "$OUTPUT_FOLDER"/reverse_complement/*.gz;
do
echo ">Filtering variants for ${sample##*/}"
java -jar "$GENOMEANALYSISTK_LOCATION" -T VariantFiltration -R "$REF_LOCATION" -V "$OUTPUT_FOLDER"/called_variants/"${sample##*/}".vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o "$OUTPUT_FOLDER"/called_variants/"${sample##*/}".filtered.vcf 
done

#Move output to final destination
mkdir -p "$HOME"/analyzed/ || exit
echo ">Placing files in $HOME/analyzed/"
mv "$OUTPUT_FOLDER" "$HOME"/analyzed/

echo "Started at $START"  
echo -n "Finished at " && date
exit
