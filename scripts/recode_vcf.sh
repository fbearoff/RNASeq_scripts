#1/bin/bash
input="$1"
chr=$2
filter=$3

echo "Input file is '$input', Chromosome is $chr"

vcftools --vcf $input --chr $chr --out "${input%%.*}.Chr$chr" --recode --recode-INFO-all
