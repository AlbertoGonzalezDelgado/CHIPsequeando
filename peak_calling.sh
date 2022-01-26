#!/bin/bash

#It is a script for peak calling in ChIP-seq analysis. It will run a R script for getting target genes. Then it will get DNA motifs.
#Authors: Maria Ada Comparin, Ricardo Daniel de Arellano Casquete de Prado & Alberto Gonzalez Delgado.
#Date: December 2021
## Accessing results folder
NUMREP=$1
FD=$2
FS=$3
UPSTREAM=$4
DOWNSTREAM=$5
PVALUE=$6
LENGTH=$7
SIZE=$8

cd ../results

## Calling peaks
i=1
while [ $i -le $NUMREP ]
do
	macs2 callpeak -t ../samples/chip_${i}/chip_${i}.bam -c ../samples/control_${i}/control_${i}.bam -f BAM --outdir . -n peaks_${i}
	((i++))
done
echo "Peaks called"
echo ""
echo ""

## Merging data
#Narrowpeak
if [ $NUMREP -gt 1 ]
then
	intersectBed -a peaks_1_peaks.narrowPeak -b peaks_2_peaks.narrowPeak > peaks.narrowPeak
	echo "peaks_1_peaks.narrowPeak and peaks_2_peaks.narrowPeak have been merged"
	for k in `seq 3 $NUMREP`
	do
		intersectBed -a peaks.narrowPeak -b peaks_${k}_peaks.narrowPeak > peaks_2.narrowPeak
		rm peaks.narrowPeak
		mv peaks_2.narrowPeak peaks.narrowPeak
		echo "peaks.narrowPeak and peaks_${k}_peaks.narrowPeak have been merged"
	done
else
	mv peaks_1_peaks.narrowPeak peaks.narrowPeak
fi
echo "*_peaks.narrowPeak have been merged"
echo ""

## Searching for DNA motifs
echo ""
echo ""
echo "Searching for DNA motifs in $SIZE bp around TSS. Please, be patient. It will take long time"
findMotifsGenome.pl peaks.narrowPeak ../genome/genome.fa ./motifs_narrowPeak -size $SIZE -len $LENGTH
echo ""
echo "DNA motifs from peaks.narrowPeak have been found."

## Searching for target genes
cp $FS/target_genes_annotation.R .
echo ""
echo ""
echo "R scripts running"
Rscript target_genes_annotation.R peaks.narrowPeak $FD $UPSTREAM $DOWNSTREAM $PVALUE --verbose
rm target_genes_annotation.R

## Showing the analysis is finished
echo "          _______________________________________________ " >> ./finishing.txt
echo " ________|                                               |_______ " >> ./finishing.txt
echo " \       |               Analysis finished               |      / " >> ./finishing.txt
echo "  \      |               CONGRATULATIONS!!               |     / " >> ./finishing.txt
echo "  /      |_______________________________________________|     \ " >> ./finishing.txt
echo " /__________)                                        (__________\ " >> ./finishing.txt
echo ""





