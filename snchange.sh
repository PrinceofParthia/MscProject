#!/bin/bash

for file in *.bam
do
    filename=$(echo $file | cut -d "." -f 1)
    samtools view -H $file |
    sed -e 's/SN:ENA|AL123456|AL123456.3/SN:AL123456.3/' |
    samtools reheader - $file > "fixed_${filename}.bam"
done

for file in SRR*.bam
do 
    rm $file
done
