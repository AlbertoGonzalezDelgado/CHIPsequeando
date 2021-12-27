# CHIPsequeando

## What is CHIPsequeando?
CHIPsequeando is an automatic computational workflow specifically designed for the analysis of ChIP-seq data.

## How to install CHIPsequeando?
Download the code from Github to the folder you want. For example: 

```
cd
mkdir ChIPseq
cd ChIPseq
git clone ttps://github.com/AlbertoGonzalezDelgado/CHIPsequeando/ 
```

## Usage
On the folder: param_input there is a file where it is neccesary to write this parameters:
1. Working directory: This parameter specifies the directory where the output folder containing the results of the ChIP-seq data analysis will be generated. For example: *home/biohacker/* 
2. folder_name: This parameter specifies the name of the output folder. For example: *home/biohacker/work* 
3. genome: This parameter specifies the directory that contains the genome file that will be used in the ChIP-seq data analysis. For example: *home/biohacker/experiment/genome/genome.fq.gz* 
4. Annotation: This parameter specifies the directory that contains the genome annotation that will be used in the ChIP-seq data analysis.For example: *home/biohacker/experiment/annotation/annotation.fq.gz*  
5. Scripts:This parameter specifies the directory that contains the scripts of CHIPsequeando. For example: *home/ChIPseq/CHIPsequeando* 
6. Number of samples: This parameter specifies the number of samples 
7. Number of copies: This parameter specifies the number of copies per each sample
8. Chip_x: This parameter specifies the directory that contains the chips samples. It is neccesary to specified as chips_x as chips samples have to be analysed. For example: * 
* Chip_1: home/biohacker/experiment/samples//chip/chip_1.fq.gz  
* Chip_1: home/biohacker/experiment/samples//chip/chip_1.fq.gz *
10. Control_home/biohacker/experiment/annotation/annotation.fq.gz:
11.     
You must specify if your samples are paired end (number_chain_end = 2) or single end (number_chain_end = 1)
