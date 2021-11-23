#!/bin/bash

## This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0$
## Author: Francisco J. Romero-Campero
## Date: June 2021
## Email: fran@us.es

## Accessing results folder
FD=$1
cd ../results

## Merging sample transcriptomes
stringtie --merge -G ../annotation/annotation.gtf -o stringtie_merged.gtf merge_list.txt

## Comparing our assembly with the reference
gffcompare -r ../annotation/annotation.gtf -G -o comparison stringtie_merged.gtf

## Showing the analysis is finished
echo "Analysis finished. CONGRATULATIONS!" >> ./finishing.txt
