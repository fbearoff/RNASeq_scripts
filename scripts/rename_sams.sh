#!/bin/bash

for f in *.sorted.bam; do
    echo "$f -> ${f%*.bam}.sam"
    mv -- "$f" "${f%*.bam}.sam"
done
