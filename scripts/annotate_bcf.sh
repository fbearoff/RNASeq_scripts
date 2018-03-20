#!/bin/sh

ANNOTATION="$1"
BCF_DIR="$2"

if [ "$#" -ne 2 ]; then
    echo "Annotate a directory of BCFs with rsIDs"
    echo "sh annotate_bcf.sh 'annotations.vcf.gz' 'BCF directory'"
    exit
fi

for bcf in $BCF_DIR/*.bcf; do
    bcf_short=${bcf##*/}
    bcftools annotate -a "$ANNOTATION" -c ID -Ob -o "$BCF_DIR"/"${bcf_short%%.*}"_rsid.bcf $bcf
    bcftools index "$BCF_DIR"/"${bcf_short%%.*}"_rsid.bcf
done
