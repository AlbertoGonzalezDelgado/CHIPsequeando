#!/bin/bash

#It is a script for ChIP-seq sample processing.
#Authors: Maria Ada Comparin, Ricardo Daniel de Arellano Casquete de Prado & Alberto Gonzalez Delgado.
#Date: December 2021

SAMPLEDIR=$1
i=$2
NUMREP=$3
FS=$4
TYPE=$5
PER_END=$6
FD=$7
NUMSAM=$8
UPSTREAM=$9
DOWNSTREAM=${10}
PVALUE=${11}
LENGTH=${12}
SIZE=${13}

## Accessing sample folder 
cd $SAMPLEDIR

## Sample quality control, mapping to reference genome and generating sorted bam file
if [ $PER_END -eq 1 ]
then
	if [ $TYPE -eq 1 ]
	then
		fastqc control_${i}.fq.gz
		bowtie2 -x ../../genome/index -U control_${i}.fq.gz -S control_${i}.sam
		samtools sort -o control_${i}.bam control_${i}.sam
		rm control_${i}.sam
		samtools index control_${i}.bam
	elif [ $TYPE -eq 2 ]
	then
		fastqc chip_${i}.fq.gz
		bowtie2 -x ../../genome/index -U chip_${i}.fq.gz -S chip_${i}.sam
		samtools sort -o chip_${i}.bam chip_${i}.sam
		rm chip_${i}.sam
		samtools index chip_${i}.bam
	fi
elif [ $PER_END -eq 2 ]
then
	if [ $TYPE -eq 1 ]
	then
		fastqc control_${i}_1.fq
		fastqc control_${i}_2.fq
		bowtie2 -x ../../genome/index -1 control_${i}_1.fq -2 control_${i}_2.fq -S control_${i}.sam
		samtools sort -o control_${i}.bam control_${i}.sam
		rm control_${i}.sam
		samtools index control_${i}.bam
	elif [ $TYPE -eq 2 ]
	then
		fastqc chip_${i}_1.fq
		fastqc chip_${i}_2.fq
		bowtie2 -x ../../genome/index -1 chip_${i}_1.fq -2 chip_${i}_2.fq -S chip_${i}.sam
		samtools sort -o chip_${i}.bam chip_${i}.sam
		rm chip_${i}.sam
		samtools index chip_${i}.bam
	fi
fi

echo "Quality control has been completed successfully"
echo "Mapping has been completed successfully"

## Write onto blackboard
if [ $TYPE -eq 1 ]
then
	echo "Finished processing control_${i}" >> ../../results/blackboard
elif [ $TYPE -eq 2 ]
then
	echo "Finished processing chip_${i}" >> ../../results/blackboard
else
	echo ""
	echo "Fatal error:"
	echo "Type of sample has not been specified properly"
	echo ""
  exit
fi

## Check if it is the last sample by reading the blackboard
cd ../../logs
BB=$(wc -l ../results/blackboard | awk '{ print $1 }')
if [ $BB -eq $NUMSAM ]
then
	sbatch $FS/peak_calling.sh $NUMREP $FD $FS $UPSTREAM $DOWNSTREAM $PVALUE $LENGTH $SIZE
	echo "All samples have been processed"
else
	echo "$BB sample(s) processed"
fi
