#1/bin/bash
input="$1"
chr=$2

echo $input $chr

vcftools --vcf $input --chr $chr --out "${input%%.*}.Chr$chr" --recode --recode-INFO-all
