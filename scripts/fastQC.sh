#!/bin/sh

INPUT_DIR="$1"
OUTPUT_DIR="$2"

mkidr -p "$OUTPUT_DIR" || exit

for file in "$INPUT_DIR"; do
    fastqc $file --outdir="$OUTPUT_DIR"
done
