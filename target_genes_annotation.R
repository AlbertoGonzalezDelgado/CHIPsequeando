#It is a Rscript for getting target genes and Gene Ontology using ChIP-seq data.
#Authors: Maria Ada Comparin, Ricardo Daniel de Arellano Casquete de Prado & Alberto Gonzalez Delgado.
#Date: January 2022

## Loading parameters
args <- commandArgs(trailingOnly = TRUE)

peaks_file <- as.character(args[[1]])
FD <- as.character(args[[2]])
UPSTREAM <- as.numeric(args[[3]])
DOWNSTREAM <- as.numeric(args[[4]])
PVALUE <- as.numeric(args[[5]])

## Loading packages 
library(ChIPseeker)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
library("clusterProfiler")
library(org.At.tair.db)

## Reading peak files
peaks <- readPeakFile(peakfile = peaks_file, header=FALSE)

## Defining promoters around TSS
promoter <- getPromoters(TxDb=txdb, upstream=UPSTREAM, downstream=DOWNSTREAM)

## Peaks annotation
peakAnno <- annotatePeak(peak = peaks, tssRegion=c(-UPSTREAM, DOWNSTREAM), TxDb=txdb)

## Visualizing plots
plotAnnoPie(peakAnno)
plotDistToTSS(peakAnno, title="Distribution of genomic loci relative to TSS", ylab = "Genomic Loci (%) (5' -> 3')")

## Constructing data frame
annotation <- as.data.frame(peakAnno)
target.genes <- annotation$geneId[annotation$annotation == "Promoter"]
write(x = target.genes, file = "target_genes.txt", sep="")

## GO terms enrichment
go_terms <- clusterProfiler::enrichGO(gene = target.genes,OrgDb = org.At.tair.db,ont = "ALL",keyType = "TAIR", pvalueCutoff = PVALUE)
goenrichment <- as.data.frame(go_terms)
promoter <- as.data.frame(promoter)
write.csv(x = goenrichment,file = "GO_terms_enrichment.csv")
write.csv(x = promoter,file = "promoter.csv")