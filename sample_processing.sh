#!/bin/bash

## This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International$
## Author: Francisco J. Romero-Campero
## Date: June 2021
## Email: fran@us.es

SAMPLEDIR=$1
i=$2
NUMSAM=$3
PIPER=$4

## Accessing sample folder
cd $SAMPLEDIR

## Sample quality control and read mapping to reference genome
fastqc sample_${i}.fq.gz
hisat2 --dta -x ../../genome/index -U sample_${i}.fq.gz -S sample_${i}.sam

## Generting sorted bam file
samtools sort -o sample_${i}.bam sample_${i}.sam
rm sample_${i}.sam
rm *.fq.gz
samtools index sample_${i}.bam

## Transcript assembly
stringtie -G ../../annotation/annotation.gtf -o sample_${i}.gtf -l sample_${i} sample_${i}.bam

## Preparing merge list file for transcriptome merging
echo ${SAMPLEDIR}/sample_${i}.gtf >> ../../results/merge_list.txt

## Gene Expression Quantification
stringtie -e -B -G ../../annotation/annotation.gtf -o sample_${i}.gtf sample_${i}.bam

## PARA IGV SE MANTIENE EL .BAM. Si se prefiere se cambia el .bam por .bw
## bamCoverage -bs 10 --normalizeUsing CPM --bam sample<N>.bam -o sample<N>.bw
## rm *.bam

## Write onto blackboard
echo "Finished processing sample_${i}" >> ../../results/blackboard

## Check if it is the last sample by reading the blackboard
cd ../../logs
BB=$(wc -l ../results/blackboard | awk '{ print $1 }')
if [ $BB -eq $NUMSAM ]
then
  sbatch $PIPER/transcriptome_merging.sh $FD
  echo "All samples have been processed"
fi





