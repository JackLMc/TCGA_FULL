# A script to reduce the matrix genotypes file to only those SNPs that are in Toju's analysis
library(tidyverse)

## Create the SNPs_to_call.txt file
SNPs_to_read <- read.csv("./Data/SNPs/Toju/SNP_Toju.csv") # these are rs numbers... Need to do this on the annotation matrix...
write.table("./Output/SNP/SNPs_to_call.txt", x = SNPs_to_read$SNP, quote = F, row.names = F)


## Create the bam_file_norms.txt
library(tidyverse)
setwd("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/6_TCGA_Full/TCGA_FULL/")
bam_files <- list.files(path = "/Volumes/2018/beggsa-tcgacolorectal/data/chr38/dna/", pattern = ".bam$")

write.table("./Output/SNP/bam_files.txt", x = bam_files, quote = F, row.names = F)

### Read in the bam_list
bam_files <- read.delim("./Output/SNP/bam_files.txt", header = F)

### Get the list of normal samples
converter <- read.delim("./Data/Important/GDC_large_mapping_TCGA.txt")

My_samp <- converter[converter$file_name %in% bam_files$V1, ] %>% droplevels()
normal_samps <- droplevels(subset(My_samp, cases.0.samples.0.sample_type == "Blood Derived Normal" | cases.0.samples.0.sample_type == "Solid Tissue Normal"))

bam_files_normals <- bam_files[bam_files$V1 %in% normal_samps$file_name, ]
write.table("./Output/SNP/bam_files_norm.txt", x = bam_files_normals, quote = F, row.names = F)


## Make the Clusters.txt file
library(tidyverse)
converter <- read.delim("./Data/Important/GDC_large_mapping_TCGA.txt")
converter$Patient.ID <- gsub("-", ".", converter$cases.0.submitter_id)
normal_samps <- droplevels(subset(converter, cases.0.samples.0.sample_type == "Blood Derived Normal" | cases.0.samples.0.sample_type == "Solid Tissue Normal"))

bam_files_norm <- read.delim("./Output/SNP/bam_files_norm.txt", header = F)
My_samp <- normal_samps[normal_samps$file_name %in% bam_files_norm$V1, ] %>% droplevels()


dupes <- My_samp$Patient.ID[duplicated(My_samp$Patient.ID)]

dupes1 <- normal_samps[normal_samps$Patient.ID %in% dupes, ]


# droplevels(subset(My_samp, Patient.ID == "TCGA.A6.2682"))


pat_sub <- read.csv("./Output/Patient_Subtypes_09_03.csv")
clusters <- merge(My_samp[, c("file_name", "Patient.ID", "cases.0.samples.0.sample_type")], pat_sub[, c("Patient.ID", "Subtype")], all.x = T)

pat_sub[pat_sub$Patient.ID %in% clusters[is.na(clusters$Subtype), ]$Patient.ID, ] # They aren't in the file...

clusters <- merge(My_samp[, c("file_name", "Patient.ID")], pat_sub[, c("Patient.ID", "Subtype")])

library(reshape2)
dcast(clusters, Subtype ~., length)

write.table("./Output/SNP/clusters.txt", x = clusters, quote = F, row.names = F)

# 
# Get probe names
BiocManager::install("vcfR")

tmp_vcf <- readLines("./Data/SNPs/")
tmp_vcf_data <- read.table("test.vcf", stringsAsFactors = FALSE)



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


