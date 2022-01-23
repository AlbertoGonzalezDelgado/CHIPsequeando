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
TYPE <- as.numeric(args[[6]])

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

## Peaks annotation for plots
if (TYPE==1)
{
  peakAnno1 <- annotatePeak(peak = peaks, tssRegion=c(-UPSTREAM, DOWNSTREAM), TxDb=txdb)
}
if (TYPE==2)
{
  peakAnno2 <- annotatePeak(peak = peaks, tssRegion=c(-UPSTREAM, DOWNSTREAM), TxDb=txdb)
}

## Visualizing plots
if (TYPE==1)
{
  plotAnnoPie(peakAnno1)
  plotDistToTSS(peakAnno1, title="Distribution of genomic loci relative to TSS (summitsbed)", ylab = "Genomic Loci (%) (5' -> 3')")
}
if (TYPE==2)
{
  plotAnnoPie(peakAnno2)
  plotDistToTSS(peakAnno2, title="Distribution of genomic loci relative to TSS (narrowPeak)", ylab = "Genomic Loci (%) (5' -> 3')")
}

## Constructing data frame
annotation <- as.data.frame(peakAnno)
target.genes <- annotation$geneId[annotation$annotation == "Promoter"]

if (TYPE==1)
{
  write(x = target.genes, file = "target_genes_summitsbed.txt", sep="")
}
if (TYPE==2)
{
  write(x = target.genes, file = "target_genes_narrowPeak.txt", sep="")
}

## GO terms enrichment
go_terms <- clusterProfiler::enrichGO(gene = target.genes,OrgDb = org.At.tair.db,ont = "ALL",keyType = "TAIR", pvalueCutoff = PVALUE)
goenrichment <- as.data.frame(go_terms)
promoter <- as.data.frame(promoter)

if (TYPE==1)
{
  write.csv(x = goenrichment,file = "GO_terms_enrichment_summitsbed.csv")
  write.csv(x = promoter,file = "promoter_summitsbed.csv")
}
if (TYPE==2)
{
  write.csv(x = goenrichment,file = "GO_terms_enrichment_narrowPeak.csv")
  write.csv(x = promoter,file = "promoter_narrowPeak.csv")
}
