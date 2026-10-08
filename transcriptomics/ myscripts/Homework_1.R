#### Homework 1: Changes in Gene Expression Between F0 and F4 ####
setwd("~/projects/eco_genomics_2026/transcriptomics/mydata")

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(vsn)  

#Import counts matrix 
countsTable <- read.table("salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
#Round counts matrix for DESeq2
countsTableRound <- round(countsTable)
#Import sample descriptions table
conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
#normalize counts table data
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ treatment)
  # 130580 genes and 38 samples total
#filter so that there are at least 15 counts per gene in 75% of samples
dds <- dds[rowSums(counts(dds) >= 15) >= 28,]

## Subset the DESeqDataSet to FO
dds_F0 <- subset(dds, select = generation == 'F0')
  # 25260 genes and 12 samples in F0
# Perform DESeq2 analysis on the subset
dds_F0 <- DESeq(dds_F0)
#Check DESeq2 results for F0
resultsNames(dds_F0)
#Name FO OWA vs AM result
F0_res_OWAvsAM <- results(dds_F0, name="treatment_OWA_vs_AM", alpha=0.05)
  #9.3% are up-regulated and 6.2% are down-regulated
#Name F0 OW vs AM
F0_res_OWvsAM <- results(dds_F0, name="treatment_OW_vs_AM", alpha=0.05)
  #14% are up-regulated and 8.1% are down-regulated
#Name F0 OA vs AM
F0_res_OAvsAM <- results(dds_F0, name="treatment_OA_vs_AM", alpha=0.05)
  #1.7% are up-regulated and 0.7% are down-regulated

## Subset the DESeqDataSet to F4
dds_F4 <- subset(dds, select = generation == 'F4')
dim(dds_F4)
# 25260 genes and 12 samples in F4
# Perform DESeq2 analysis on the subset
dds_F4 <- DESeq(dds_F4)
#Check DESeq2 results for F4
resultsNames(dds_F4)
#Name F4 OWA vs AM result
F4_res_OWAvsAM <- results(dds_F4, name="treatment_OWA_vs_AM", alpha=0.05)
summary(F4_res_OWAvsAM)
#0.41% are up-regulated and 0.55% are down-regulated
#Name F4 OW vs AM
F4_res_OWvsAM <- results(dds_F4, name="treatment_OW_vs_AM", alpha=0.05)
summary(F4_res_OWvsAM)
#0.44% are up-regulated and 0.17% are down-regulated
#Name F4 OA vs AM
F4_res_OAvsAM <- results(dds_F4, name="treatment_OA_vs_AM", alpha=0.05)
summary(F4_res_OAvsAM)
#0.42% are up-regulated and 0.2% are down-regulated