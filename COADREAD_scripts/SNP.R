# A R script to investigate associations of SNPs with Patient groups
library(tidyverse)
library(UsefulFunctions)


SNP <- read.delim("./Data/SNPs/matrix_genotypes.final.txt")

SNP1 <- column_to_rownames(SNP, var = "snpID") %>% 
  as.matrix() %>%
  t()
rownames(SNP1) <- gsub("^X", "", rownames(SNP1))


manifest <- read.delim("./Data/SNPs/manifest.txt")



mapping <- read.delim("./Data/GDC_large_mapping_TCGA.txt")

manifest[manifest$filename %in% mapping$file_name,]



mapping$id <- gsub("-", ".", mapping$id)


install.packages("rjson")
library("rjson")

meta <- fromJSON(file = "~/Downloads/metadata.cart.2019-11-15.json")
head(meta[[2]])



mapping[mapping$id  %in% rownames(SNP1),]


head(rownames(SNP1))



