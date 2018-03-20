#!/bin/sh

INPUT_DIR="$1"
OUTPUT_DIR="$2"

if [ "$#" -ne 2 ]; then
    echo "Run FastQC on a directory of FastQ files"
    echo "sh fastQC.sh 'FastQ Dir' 'Output Dir'"
    exit
fi

mkdir -p "$OUTPUT_DIR" || exit

for file in $INPUT_DIR; do
    fastqc "$file" --outdir="$OUTPUT_DIR"
done
