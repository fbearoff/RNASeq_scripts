#!/bin/bash

for vcf in *.filtered.vcf; do
    bcftools query -i'%FILTER="PASS"' -f'%CHROM %POS %QUAL\n' $vcf > ${vcf%*.vcf}.lqual
done
