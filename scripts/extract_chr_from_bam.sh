#!//bin/bash

#run in folder with BAMs to extract

chr=$1

#check chromsome supplied as argument
if [ $# -eq 0 ]
  then
    echo "No chromosome supplied as argument, exiting!" && exit
fi

#extract desired chromosome to new BAM
mkdir ./$chr
for bam in *.bam; do
    samtools view -b $bam $chr > ./$chr/$bam
done

#index extracted BAMs
for bam in ./$chr/*.bam; do
    samtools index $bam
done
