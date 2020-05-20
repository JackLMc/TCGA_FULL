# A script to reduce the matrix genotypes file to only those SNPs that are in Toju's analysis
library(tidyverse)

## Read in required data
SNPs_to_read <- read.csv("./Data/SNPs/Toju/SNP_Toju.csv") # these are rs numbers... Need to do this on the annotation matrix...
write.table("./Output/SNP/SNPs_to_call.txt", x = SNPs_to_read$SNP, quote = F, row.names = F)


## Get the Bam list

setwd("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/6_TCGA_Full/TCGA_FULL/")
bam_files <- list.files(path = "/Volumes/2018/beggsa-tcgacolorectal/data/chr38/dna/", pattern = ".bam$")
write.table("./Output/SNP/bam_files.txt", x = bam_files, quote = F, row.names = F)


# Read in the bam_list
bam_files <- read.delim("./Output/SNP/bam_files.txt", header = F)

# Get the list of normal samples
converter <- read.delim("./Data/Important/GDC_large_mapping_TCGA.txt")

My_samp <- converter[converter$file_name %in% bam_files$V1, ] %>% droplevels()
normal_samps <- droplevels(subset(My_samp, cases.0.samples.0.sample_type == "Blood Derived Normal" | cases.0.samples.0.sample_type == "Solid Tissue Normal"))

bam_files_normals <- bam_files[bam_files$V1 %in% normal_samps$file_name, ]
write.table("./Output/SNP/bam_files_norm.txt", x = bam_files_normals, quote = F, row.names = F)

# 


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


