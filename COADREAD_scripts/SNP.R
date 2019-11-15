# A R script to investigate associations of SNPs with Patient groups
library(tidyverse)
library(UsefulFunctions)


SNP <- read.delim("./Data/SNPs/matrix_genotypes.final.txt")

SNP1 <- column_to_rownames(SNP, var = "snpID") %>% 
  as.matrix() %>%
  t()
rownames(SNP) <- gsub("^X", "", rownames(SNP))




