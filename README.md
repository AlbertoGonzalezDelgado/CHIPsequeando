# CHIPsequeando

## What is CHIPsequeando?
CHIPsequeando is an automatic computational workflow specifically designed for the analysis of ChIP-seq data.

## How to install CHIPsequeando?
Download the code from Github to the folder you want. For example: 

{% filename %}

cd
mkdir ChIPseq
cd ChIPseq
git clone ttps://github.com/AlbertoGonzalezDelgado/CHIPsequeando/ 
{% endfilename %}
## Usage
On the param_input there is a file where it is neccesary to write this parameters:
1. Working directory: the directory of the folder you want 
You must specify if your samples are paired end (number_chain_end = 2) or single end (number_chain_end = 1)
