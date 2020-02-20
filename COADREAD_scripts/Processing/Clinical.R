# A script to clea and process clinical data from the COADREAD project of TCGA
## Packages
library(UsefulFunctions)
library(tidyverse)
library(ggpubr)

# Collate the patients I want to look for (taken from Counts_to_CIRC.R)
Patient_list <- read.delim("./Output/Patient_list.txt")
Patient_list <- levels(Patient_list$x)

### Method of determining MMR status
# Take old data from TCGA_pub_clinical.csv 
tcga_pub_clinical <- read.csv("./Data/Clinical_data/tcga_pub_clinical.csv")

# From a paper: 
# Pan-cancer immunogenomic perspective on the tumor microenvironment based on PD-L1 and CD8 T cell infiltration
Paper_clinical <- read.csv("./Data/Clinical_data/Paper_Clinical.csv")
Paper_clinical$Patient.ID <- samptopat(Paper_clinical$TCGA.ID)
Paper_clinical$Patient.ID <- gsub("-", ".", Paper_clinical$Patient.ID)

# # uneeded, but look at shared pats.
# capp <- tcga_pub_clinical[tcga_pub_clinical$Patient.ID %in% Paper_clinical$Patient.ID,]

## find the patients in the paper clinical data
clin_app <- Paper_clinical[Paper_clinical$Patient.ID %in% Patient_list, ] # Finds 332 Patients
COADREAD <- droplevels(subset(Paper_clinical, Cancer.type == "COAD" | Cancer.type == "READ"))

## Find the remaining patients
missing_pats <- Patient_list[!'%in%'(Patient_list, COADREAD$Patient.ID)]
pats_tcga_pub <- tcga_pub_clinical[tcga_pub_clinical$Patient.ID %in% missing_pats,]

# Collate the data of interest
head(tcga_pub_clinical)
paper_data <- clin_app[, c("Patient.ID", "Microsatellite.instability")]
cbio_data <- pats_tcga_pub[, c("Patient.ID", "MSI_STATUS")]
colnames(paper_data) <- c("Patient.ID", "MSI_STATUS")
full_Data <- rbind(cbio_data, paper_data)

FD <- full_Data[!duplicated(full_Data), ]
FD1 <- FD[!is.na(FD$MSI_STATUS),] %>% droplevels()
FD1$MSI_STATUS <- ifelse((FD1$MSI_STATUS == "MSI-H"), "MSI-H", "MSS") ## Calls MSI-L patients MSS

count(FD1, vars = c("MSI_STATUS"))

write.csv("./Output/Clinical_Data_540.csv", x = FD1, row.names = F)

#### END ####

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
