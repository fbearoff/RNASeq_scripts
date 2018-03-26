#!/bin/sh

if [ "$#" -ne 1 ]; then
    echo "Convert a directory of BAMs to SAMs"
    echo "sh bam_to_sam.sh 'bam directory'"
    exit
fi

bam_dir="$1"

mkdir "$bam_dir"/sams

for bam in "$bam_dir"/*/*.bam; do
    bam_sn="${bam##*/}"
    sambamba view $bam -h -t 32 -o "$bam_dir"/sams/"${bam_sn%%.bam}".sam
done
