# A script to reduce the matrix genotypes file to only those SNPs that are in Toju's analysis
library(tidyverse)

## Read in required data
SNPs_to_read <- read.csv("./Data/SNPs/Toju/SNP_Toju.csv")$SNP # these are rs numbers... Need to do this on the annotation matrix...
get <- read.delim("./Data/SNPs/14_02/annotation/SNP_annotations.hg19.txt")[, c("Probe_Set_ID", "dbSNP_RS_ID")]

these_probes <- get[get$dbSNP_RS_ID %in% SNPs_to_read, ]

Toju_sig <- c("rs4577037", "rs4500045",
              "rs859", "rs4500045",
              "rs6692729", "rs6692729",
              "rs4500045", "rs13331952",
              "rs2291299", "rs4796105", "rs13331952")


these_probes[these_probes$dbSNP_RS_ID %in% Toju_sig, ]




## Subset input dataframe
SNP <- read.delim("./Data/SNPs/matrix_genotypes.final.txt")
head(SNP)[, 1:10]
SNP1 <- column_to_rownames(SNP, var = "snpID") %>%
  as.matrix() %>%
  t()
rownames(SNP1) <- gsub("^X", "", rownames(SNP1))

SNP2 <- SNP1[, colnames(SNP1) %in% these_probes$Probe_Set_ID]

rownames(SNP2) <- gsub("\\.", "-", rownames(SNP2))
SNP3 <- t(SNP2)

# Write out new input dataframe
write.table("./Output/SNP/matrix_genotypes.subset.txt", x = SNP3, row.names = T, quote = F, sep = "\t")


