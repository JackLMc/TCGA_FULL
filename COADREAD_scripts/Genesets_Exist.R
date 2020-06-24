# A script to ensure that all genes in the Genesets I've chosen are in the Counts_cqn matrix
library(UsefulFunctions)
library(tidyverse)
library(ggpubr)

# Load data
load("./R_Data/Counts_clean.RData")

# The Genesets
## Cell type Genesets
CTGenesets <- read.csv("./Data/Genesets/Investigative_Genesets.csv", stringsAsFactors = T)

CTGenesets[!'%in%'(CTGenesets$HUGO.symbols, rownames(Counts_cqn)), ]
rownames(Counts_cqn)[grep("CD4", rownames(Counts_cqn))]


# Signature Genesets
SigGen <- read.csv("./Data/Genesets/Signature_Genesets.csv", stringsAsFactors = T)

SigGen[!'%in%'(SigGen$HUGO.symbols, rownames(Counts_cqn)), ]
rownames(Counts_cqn)[grep("CD4", rownames(Counts_cqn))]



# GO Terms



