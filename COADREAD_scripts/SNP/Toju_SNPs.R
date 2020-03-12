# A script to reduce the matrix genotypes file to only those SNPs that are in Toju's analysis
library(tidyverse)

## Read in required data
SNPs_to_read <- read.csv("./Data/SNPs/Toju/SNP_Toju.csv")$SNP # these are rs numbers... Need to do this on the final matrix...
get <- read.delim("./Data/SNPs/14_02/annotation/SNP_annotations.hg19.txt")[, c("Probe_Set_ID", "dbSNP_RS_ID")]

these_probes <- get[get$dbSNP_RS_ID %in% SNPs_to_read, ]


## Subset input dataframe
SNP <- read.delim("./Data/SNPs/matrix_genotypes.final.txt")

SNP1 <- column_to_rownames(SNP, var = "snpID") %>%
  as.matrix() %>%
  t()
rownames(SNP1) <- gsub("^X", "", rownames(SNP1))

SNP2 <- SNP1[, colnames(SNP1) %in% these_probes$Probe_Set_ID]


# Write out new input dataframe
write.table("./Output/SNP/matrix_genotypes.subset.txt", x = SNP2, row.names = F, quote = F, sep = "\t")

