# A script to clea and process clinical data from the COADREAD project of TCGA
## Packages
library(UsefulFunctions)
library(tidyverse)
library(ggpubr)

# Collate the patients I want to look for (taken from Counts_to_CIRC.R)
Patient_list <- read.delim("./Output/Patient_list.txt")
Patient_list <- levels(Patient_list$x)

# Clinical data from the Broad GDAC FireHose
tcga_firehose <- read.delim("./Data/Clinical_data/COADREAD.clin.merged.txt", header = F)

## Take these columns, edit later if I need more data, patient.bcr_patient_barcode NEEDS to be first
tcga_firehose2 <- tcga_firehose[tcga_firehose$V1 %in% c("patient.bcr_patient_barcode", 
                      "patient.microsatellite_instability_test_results.microsatellite_instability_test_result.mononucleotide_and_dinucleotide_marker_panel_analysis_status",
                      "patient.weight",
                      "patient.gender",
                      "patient.height"), ]

## Define a function which makes the first row the column name
header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

## Rename column
tcga_firehose3 <- header.true(tcga_firehose2)
colnames(tcga_firehose3)[colnames(tcga_firehose3) == "patient.bcr_patient_barcode"] <- "Parameter"

# Reorder data so that Patient Identifiers are a column
tcga_firehose4 <- tcga_firehose3 %>% gather(contains("tcga"), key = "Patient.ID", value = "value") %>%
  spread(., key = "Parameter", value = "value")

## Edit the data to be in the format of the other stuff
tcga_firehose4$Patient.ID <- toupper(tcga_firehose4$Patient.ID)
tcga_firehose4$Patient.ID <- gsub("-", ".", tcga_firehose4$Patient.ID)
colnames(tcga_firehose4) <- c("Patient.ID", "Gender", "Height", "MSI_STATUS", "Weight")

tcga_firehose4$MSI_STATUS <- toupper(tcga_firehose4$MSI_STATUS)
tcga_firehose4$MSI_STATUS <- ifelse((tcga_firehose4$MSI_STATUS == "MSI-H"), "MSI-H", "MSS") # Rename MSI-L; MSS

# Only take the patients I need - Count data available.
tcga_FH <- tcga_firehose4[tcga_firehose4$Patient.ID %in% Patient_list, ]
tcga_FH$MSI_STATUS <- as.factor(tcga_FH$MSI_STATUS)


reshape2:: dcast(tcga_FH, MSI_STATUS ~., length)

# Finds two more patients...
# Write the data out.
write.csv("./Output/Clinical_Data_542.csv", x = tcga_FH, row.names = F)

#### END ####
