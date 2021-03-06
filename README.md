# CHIPsequeando

## What is CHIPsequeando?
CHIPsequeando is a computational pipeline designed for the automatic analysis of ChIP-seq data from *Arabidopsis thaliana*. From raw fastaq files, CHIPsequeando will provide a **list of target genes** that could be regulated by the transcription factor, **GO terms** results providing the functions regulated and **DNA motifs in TSS region** where the transcription factor is bound.

## How to install CHIPsequeando?
CHIPsequeando requires the following dependencies that should be installed previously:
* [FastQC](https://bioinformatics.babraham.ac.uk/projects/fastqc/)
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [Samtools](http://www.htslib.org/)
* [MACS2](https://pypi.org/project/MACS2/)
* [Slurm](https://slurm.schedmd.com/documentation.html)
* [HOMER](http://homer.ucsd.edu/homer/introduction/install.html)
* [R](https://www.r-project.org/)
  * [ChIPseeker](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html)
  * [TxDb.Athaliana.BioMart.plantsmart28](https://bioconductor.org/packages/release/data/annotation/html/TxDb.Athaliana.BioMart.plantsmart28.html)
  * [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
  * [org.At.tair.db](https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html)

Download the code from Github into the folder desired. For example: 
```
cd
mkdir ChIPseq
cd ChIPseq
git clone ttps://github.com/AlbertoGonzalezDelgado/CHIPsequeando/ 
```

## How to use CHIPsequeando
```
chip_data_process.sh <param_input>
```
To process ChIP-Seq data run the executable **chip_data_process.sh** with [param_input](param_input/param_input.txt) as unique input.
On the folder: param_input you can  find the file *param_input.txt* where it is neccesary to write this parameters:
1. **working_directory:** This parameter specifies the directory where the output folder containing the results will be generated. For example: */home/user/* 
2. **folder_name:** This parameter specifies the name of the output folder. For example: *work* 
3. **genome:** This parameter specifies the directory that contains the genome file that will be used in the ChIP-seq data analysis. This file has to be previously downloaded. For example: */home/user/experiment/genome/genome.fq.gz* 
4. **annotation:** This parameter specifies the directory that contains the genome annotation that will be used in the ChIP-seq data analysis. This file has to be previously downloaded. For example: */home/user/experiment/annotation/annotation.fq.gz*  
5. **scripts:** This parameter specifies the directory that contains the scripts of CHIPsequeando. For example: */home/user/ChIPseq/CHIPsequeando* 
6. **number_samples:** This parameter specifies the number of samples 
7. **number_copy:** This parameter specifies the number of copies per each sample.
8. **number_chain_end:** This parameter can take the values **1** when your data is **single-end** or **2** when your data is **paired-end** (two .fastq file per sample; for example *control1_1.fq.gz* and *control1_2.fq.gz*). 
9. **upstream_bp:** This parameter specifies the 5' region from TSS that will be used to annotate peaks and to get promoters. For example: *1000*.
10.**downstream_bp:** This parameter specifies the 3' region from TSS that will be used to annotate peaks and to get promoters. For example: *1000*.
11. **word_length:** This parameter specifies the length of the secuence that will be used to search DNA motifs inside the peaks annotated. It can be used several numbers separated by "," for searching for words using diferent lengths. For example: *4* or *6,8*.
12. **region_size:** This parameter specifies the width of DNA region around TSS that will be used to search for DNA motifs. For example: *200*
13. **p_value:** This parameter specifies the p-value used in Gene Ontology test. For example: *0.05*  
14. **chip_x:** This parameter specifies the directory that contains the chips samples. It is neccesary to specified as chips_x as chips samples have to be analysed. For example:
* *Chip_1: /home/user/experiment/samples/chip/chip_1.fq.gz* 
* *Chip_2: /home/user/experiment/samples/chip/chip_2.fq.gz*
15. **control_x:** his parameter specifies the directory that contains the control or input samples. It is neccesary to specified as controls_x as chips samples have to be analysed. For example:
* *Control_1: /home/user/experiment/samples//control/control_1.fq.gz* 
* *Control_1: /home/user/experiment/samples//control/control_2.fq.gz*

There is given a param input file example that could be helpful
[param_example_input](param_input/param_example_file_prr5.txt).

**Warning!** It is neccesary to write a space bar after each ":" in the param input file.

## Running CHIPsequeando
Once all parameters have been specified, the next step is to run the scripts. For example: 

```
/home/user/CHIPsequeando/chip_data_process.sh /home/user/CHIPsequeando/param_input/param_file_prr5.txt
```

## Results generated by CHIPsequeando
When CHIPsequeando pipeline is running, the next folders will be created in the output directory: 
1. **annotation:** containing the *annotation.gtf* file.
2. **genome:** containing the *genome.fa* file and the genome index. 
3. **logs:** containing the intermediary and final reports written in *.out* files. There will be one more files than copies in the experiment. 
4. **results:** containing:
   * *blackboard* file, where each sample processed is noted down.
   * *finishing.txt* file, showing the analysis has been completed succesfully.
   * *GO_terms_enrichment.csv*, containig the enriched GO terms.
   * *motifs* folders, with DNA motifs around TSS provided by HOMER.
   * *peaks* files, required by the R script for determining the target genes. Those include *peaks.narrowPeak* and *peaks.bed*.
   * *promoters.csv* containing the promoters found around TSS.
   * *Rplots.pdf*, consisting of two graphics representing distribution of peaks in genome and distribution of genomic loci relative to TSS.
   * *target_genes.txt*, a list of the TF target genes.
5. **samples:** containing for each copy a *chip* and a *control* folder. Everyone contains the bowtie2 results in *.bam* and *.bam.bai*, the quality analysis results (made by *fastqc*) in *.html* and in *.zip* formats, and the sample input in *fq.gz* format.

It is given an example of the results generated in an analysis that could be helpful [results_folder](results/).
