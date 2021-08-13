#!/bin/bash
ReadsPath=/bionas1/bean/GBSplates/44/reads
while read parentInfo; do
#check if paired files exist
sample=`echo ${parentInfo}| cut -d' ' -f1`
echo ${sample}
#zcat ${ReadsPath}/${parent}*.fastq.gz | bgzip > ./../processed_data/${parent}.fastq.gz &
done <samples.txt