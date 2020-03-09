# A script to clea and process clinical data from the COADREAD project of TCGA
## Packages
library(UsefulFunctions)
library(tidyverse)
library(ggpubr)

# Collate the patients I want to look for (taken from Counts_to_CIRC.R)
Patient_list <- read.delim("./Output/Patient_list.txt")
Patient_list <- levels(Patient_list$x)

# Clinical data from the Broad GDAC FireHose
tcga_firehose <- read.delim("./Data/Clinical_Data_Firehose/COADREAD.clin.merged.txt", header = F)
tcga_FH <- column_to_rownames(tcga_firehose, var = "V1") %>% as.matrix() %>% t()
rownames(tcga_FH) <- NULL

# Remove all rows with NA across the board
na_omitted_tcga_clin <- tcga_FH[rowSums(is.na(tcga_FH)) != 629, ] %>% as.data.frame()
tcga_clin <- na_omitted_tcga_clin

## Edit the data to be in the format of the other stuff
tcga_clin$Patient.ID <- toupper(tcga_clin$patient.bcr_patient_barcode)
tcga_clin$Patient.ID <- gsub("-", ".", tcga_clin$Patient.ID)

# Only take the patients I need - Count data available.
tcga_FH <- tcga_clin[tcga_clin$Patient.ID %in% Patient_list, ]
tcga_FH_choose <- tcga_FH[, colSums(is.na(tcga_FH)) < 406] # 0.75 * 542, # 75% of patients have data

## 629 columns to choose from...
write.csv("./Output/tcga_FH_choose.csv", x = tcga_FH_choose, row.names = F)

## GO OFF AND CHOOSE THEM

## Take these columns, edit later if I need more data, patient.bcr_patient_barcode NEEDS to be first
CoI <- c("Patient.ID", 
         "patient.microsatellite_instability_test_results.microsatellite_instability_test_result.mononucleotide_and_dinucleotide_marker_panel_analysis_status",
         "patient.weight", "patient.gender", "patient.height",
         "patient.age_at_initial_pathologic_diagnosis",
         "patient.anatomic_neoplasm_subdivision", "patient.colon_polyps_present",
         "patient.days_to_birth", "patient.days_to_death", 
         "patient.days_to_last_followup", "patient.ethnicity",
         "patient.follow_ups.follow_up.vital_status", "patient.histological_type",
         "patient.lymph_node_examined_count", "patient.race_list.race", "patient.stage_event.pathologic_stage",
         "patient.tumor_samples.tumor_sample.tumor_locations.tumor_location.site_of_disease_description",
         "patient.venous_invasion", "patient.vital_status", "patient.year_of_initial_pathologic_diagnosis")

tcga_FH1 <- tcga_FH[, colnames(tcga_FH) %in% CoI]


# CoI[!'%in%'(CoI, colnames(tcga_FH1))]

head(tcga_FH1)
colnames(tcga_FH1) <- c("Patient_age_at_diagnosis", "Neoplasm_subdivision", "Polyps",
                        "Days_to_birth", "Days_to_death", "Days_to_last_followup",
                        "Ethnicity", "Followup_vital_status", "Gender", "Height", "Histological_Type",
                        "Lymph_count", "Race", "Pathologic_stage", "Venous_invasion",
                        "Vital_status", "Weight", "Year_of_diagnosis", "MSI_STATUS", "Site_of_disease", "Patient.ID")

tcga_FH1$MSI_STATUS <- toupper(tcga_FH1$MSI_STATUS)
tcga_FH1$MSI_STATUS <- ifelse((tcga_FH1$MSI_STATUS == "MSI-H"), "MSI-H", "MSS") # Rename MSI-L; MSS
tcga_FH1$MSI_STATUS <- as.factor(tcga_FH1$MSI_STATUS)


reshape2:: dcast(tcga_FH1, MSI_STATUS ~., length)


# Finds two more patients...
# Write the data out.
write.csv("./Output/Clinical_Data_542.csv", x = tcga_FH1, row.names = F)

#### END ####
