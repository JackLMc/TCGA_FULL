# A script to clea and process clinical data from the COADREAD project of TCGA
## Packages
library(UsefulFunctions)
library(tidyverse)
library(ggpubr)

# Collate the patients I want to look for (taken from Counts_to_CIRC.R)
Patient_list <- read.delim("./Output/Patient_list.txt")
Patient_list <- levels(Patient_list$x)

#### Clinical data from the Broad GDAC FireHose ####
tcga_firehose <- read.delim("./Data/Clinical_Data/COADREAD.FireHose_MMR.txt", header = F)
tcga_FH <- column_to_rownames(tcga_firehose, var = "V1") %>% as.matrix() %>% t() %>% as.data.frame()
rownames(tcga_FH) <- NULL

# Subset my patients and MMR status
tcga_FH$Patient.ID <- toupper(tcga_FH$patient.bcr_patient_barcode)
tcga_FH$Patient.ID <- gsub("-", ".", tcga_FH$Patient.ID)
tcga_FH <- tcga_FH[tcga_FH$Patient.ID %in% Patient_list, ]
tcga_FH1 <- tcga_FH[, colnames(tcga_FH) %in% c("Patient.ID", 
                                               "patient.microsatellite_instability_test_results.microsatellite_instability_test_result.mononucleotide_and_dinucleotide_marker_panel_analysis_status"
                                              )]

# Rename columns in a way that columns can be added
colnames(tcga_FH1)[colnames(tcga_FH1) %in% c("patient.microsatellite_instability_test_results.microsatellite_instability_test_result.mononucleotide_and_dinucleotide_marker_panel_analysis_status")] <- c("MSI_STATUS")
tcga_FH1$MSI_STATUS <- toupper(tcga_FH1$MSI_STATUS)
tcga_FH1$MSI_STATUS <- ifelse((tcga_FH1$MSI_STATUS == "MSI-H"), "MSI-H", "MSS") # Rename MSI-L; MSS
tcga_FH1$MSI_STATUS <- as.factor(tcga_FH1$MSI_STATUS)
reshape2:: dcast(tcga_FH1, MSI_STATUS ~., length)

#### Clinical data from cBioPortal ####
tcga_cBio <- read.delim("./Data/Clinical_Data/coadread_cBio.tsv", header = T)

# Edit columns, ensure I only have the patients I want
tcga_cBio$Patient.ID <- gsub("-", ".", tcga_cBio$Patient.ID)
tcga_cBio <- tcga_cBio[tcga_cBio$Patient.ID %in% tcga_FH1$Patient.ID, ]
tcga_CB <- droplevels(subset(tcga_cBio, Sample.Type == "Primary"))
# tcga_CB$Patient.ID[duplicated(tcga_CB$Patient.ID)]

tcga_clinical <- merge(tcga_cBio, tcga_FH1, by = "Patient.ID")


# Finds two more patients...
# Write the data out.
write.csv("./Output/Clinical_Data_542.csv", x = tcga_clinical, row.names = F)

#### END ####
