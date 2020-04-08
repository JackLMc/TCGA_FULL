# A script to clean and process clinical data from the COADREAD project of TCGA
## Packages
library(UsefulFunctions)
library(tidyverse)
library(ggpubr)

# Collate the patients I want to look for (taken from Counts_to_CIRC.R)
Patient_list <- read.delim("./Output/Patient_list.txt")
Patient_list <- levels(Patient_list$x)

#### Clinical data from the GDC BCR Biotab ####
## Follow-up folder - more recent data?
# annot_COAD <- read.delim("./Data/Clinical_Data/BCR_BioTab_GDC/Followup/annotations_COAD.txt")
# annot_READ <- read.delim("./Data/Clinical_Data/BCR_BioTab_GDC/Followup/annotations_READ.txt")
followup_coad <- read.delim("./Data/Clinical_Data/BCR_BioTab_GDC/Followup/nationwidechildrens.org_clinical_follow_up_v1.0_coad.txt")[-1:-2,]
followup_read <- read.delim("./Data/Clinical_Data/BCR_BioTab_GDC/Followup/nationwidechildrens.org_clinical_follow_up_v1.0_read.txt")[-1:-2,]

colnames_to_remove <- colnames(followup_coad)[!'%in%'(colnames(followup_coad), colnames(followup_read))] # Ensure the data has the same columns
followup_coad1 <- followup_coad[, !'%in%'(colnames(followup_coad), colnames_to_remove)]
followup <- rbind(followup_coad1, followup_read)

followup$Patient.ID <- gsub("-", ".", followup$bcr_patient_barcode) # Get the barcode in the right format
followup1 <- followup[followup$Patient.ID %in% Patient_list, ] # Subset the patients we have


Patient_subtypes <- read.csv("./Output/Patient_Subtypes_09_03.csv")[, c("Patient.ID", "Subtype")] # Merge the subtype data in
followup3 <- merge(Patient_subtypes, followup1, by = "Patient.ID") %>% droplevels()

followup3 <- lapply(followup3, function(column) as.character(column)) %>% bind_rows() %>% as.data.frame() # Make all columns characters
followup3$form_completion_date <- as.Date(followup3$form_completion_date)     

# Change some numerical columns to numerics
## last_contact_days_to, death_days_to, 
followup3$last_contact_days_to <- as.numeric(as.character(followup3$last_contact_days_to))
followup3$death_days_to <- as.numeric(as.character(followup3$death_days_to))

remove_the_pats <- followup3[is.na(followup3$last_contact_days_to) & is.na(followup3$death_days_to),]$Patient.ID
followup3 <- followup3[!'%in%'(followup3$Patient.ID, remove_the_pats), ] # TCGA.CL.4957 removes this patient - has NA data for both days to contact and death

# NA being introduced by coercing [Not Available] to a numeric
pats_with_followup <- followup3[!is.na(followup3$last_contact_days_to), ]
pats_with_followup$Patient.ID <- as.factor(pats_with_followup$Patient.ID)

pats_followup <- list() # Find the most recent data (some have two inputs)
c <- 1
for(i in levels(pats_with_followup$Patient.ID)){
  print(i)
  work <- droplevels(subset(pats_with_followup, Patient.ID == i))
  this <- work[which.max(work$last_contact_days_to), ]
  this$Patient.ID <- as.character(this$Patient.ID)
  pats_followup[[i]] <- this
  
  c <- c + 1}
pats_follow <- bind_rows(pats_followup) 

duplicated(pats_follow$Patient.ID) # They're all unique
pats_follow$Patient.ID <- as.factor(pats_follow$Patient.ID)

# The Patients who have passed - no followup data
pats_nofollow <- followup3[is.na(followup3$last_contact_days_to), ]
pats_with_two_deaths <- pats_nofollow[duplicated(pats_nofollow$Patient.ID),]$Patient.ID
choose_these <- pats_nofollow[pats_nofollow$Patient.ID %in% pats_with_two_deaths, ]
# TCGA.CM.5862 has two form completions a year apart
# TCGA.DM.A0XD has three inputs into the data

# Remove both the patients
remove_pats <- pats_nofollow[!'%in%'(pats_nofollow$Patient.ID, pats_with_two_deaths), ]

# Choose the data from these
merge_these <- droplevels(subset(choose_these, form_completion_date == "2012-1-16" | bcr_followup_barcode == "TCGA-DM-A0XD-F17210"))
## Choosing the latest completion date, or a random one from the other patient

complete_nofollow <- rbind(merge_these, remove_pats)

pats_cleaned <- rbind(complete_nofollow, pats_follow)
pats_DaA <- pats_cleaned[!is.na(pats_cleaned$last_contact_days_to) & !is.na(pats_cleaned$death_days_to),]$Patient.ID
pats_cleaned$last_contact_days_to <- ifelse((pats_cleaned$Patient.ID %in% pats_DaA & pats_cleaned$vital_status == "Dead"), NA, pats_cleaned$last_contact_days_to)

# Some patients are alive and dead in the dataframe
duplicated_pats <- pats_cleaned[pats_cleaned$Patient.ID %in% pats_cleaned[duplicated(pats_cleaned$Patient.ID), ]$Patient.ID, ]
duplicated_pats$Patient.ID <- as.factor(duplicated_pats$Patient.ID)


alive_Dead_data <- list()
c <- 1
for(i in levels(duplicated_pats$Patient.ID)){
  print(i)
  work <- droplevels(subset(duplicated_pats, Patient.ID == i))
  work$today_minus_form <-  Sys.Date() - work$form_completion_date
  these_data <- work[which.min(work$today_minus_form), ]
  these_data1 <- these_data[, !'%in%'(colnames(these_data), c("today_minus_form"))]
  alive_Dead_data[[i]] <- these_data1
  c <- c + 1
}

duplicated_cleaned <- bind_rows(alive_Dead_data)
Rem_dup_pats <- pats_cleaned[!'%in%'(pats_cleaned$Patient.ID, pats_cleaned[duplicated(pats_cleaned$Patient.ID), ]$Patient.ID), ]

## This is the final dataframe that is cleaned
Patients_Cleaned_Followup <- rbind(Rem_dup_pats, duplicated_cleaned)
rownames(Patients_Cleaned_Followup) <- Patients_Cleaned_Followup$Patient.ID

## Create Kaplan-Meier plot (survival plot)
# BiocManager::install("survcomp")
library(survcomp)
set.seed(123)
survival.time <- as.integer(ifelse(is.na(Patients_Cleaned_Followup$last_contact_days_to),
                                   Patients_Cleaned_Followup$death_days_to, Patients_Cleaned_Followup$last_contact_days_to))
censor <- ifelse(!is.na(as.integer(Patients_Cleaned_Followup$death_days_to)), 1, 0)


subtype <- as.factor(Patients_Cleaned_Followup$Subtype)

df <- data.frame(survival.time, censor, subtype)
# df1 <-droplevels(subset(df, subtype != "MSS"))

km.coxph.plot(formula.s = Surv(survival.time, censor) ~ subtype, data.s = df, x.label = "Time (days)",
              y.label = "Probability of survival", main.title="",
              leg.text = c("MSI-H", "MSS", "MSS-hiCIRC"), leg.pos = "topright", leg.inset = 0,
              .col = c("#56B4E9", "#009E73", "#999999"),
              .lty = c(1, 1, 1, 1, 1), show.n.risk = TRUE, n.risk.step = 1000, n.risk.cex = 0.85,
              verbose = FALSE)

# No difference in overall survival, apparently you don't see the difference in MSI-H because they generally do a lot worse if they relapse
## A better proxy would be Relapse free survival - not sure how I'm going to get this...

Patients_Cleaned_Followup


# What other data is around?
## Clinical Folder
# annot_COAD <- read.delim("./Data/Clinical_Data/BCR_BioTab_GDC/Clinical/annotations_COAD.txt")
# annot_READ <- read.delim("./Data/Clinical_Data/BCR_BioTab_GDC/Clinical/annotations_READ.txt")
clinical_coad <- read.delim("./Data/Clinical_Data/BCR_BioTab_GDC/Clinical/nationwidechildrens.org_clinical_patient_coad.txt")[-1:-2, ]
clinical_read <- read.delim("./Data/Clinical_Data/BCR_BioTab_GDC/Clinical/nationwidechildrens.org_clinical_patient_read.txt")[-1:-2, ]

colnames_to_remove <- colnames(clinical_coad)[!'%in%'(colnames(clinical_coad), colnames(clinical_read))]
clinical_coad1 <- clinical_coad[, !'%in%'(colnames(clinical_coad), colnames_to_remove)]
clinical <- rbind(clinical_coad1, clinical_read)
clinical$Patient.ID <- gsub("-", ".", clinical$bcr_patient_barcode)
clinical1 <- clinical[clinical$Patient.ID %in% Patient_list, ]

rownames(clinical1) <- clinical1$Patient.ID
Patient_subtypes <- read.csv("./Output/Patient_Subtypes_09_03.csv")[, c("Patient.ID", "Subtype")]

clinical3 <- merge(Patient_subtypes, clinical1, by = "Patient.ID")


tcga_FH$Patient.ID[!'%in%'(tcga_FH$Patient.ID, clinical3$Patient.ID)]
droplevels(subset(tcga_FH, Patient.ID == "TCGA.F5.6810"))[, 1:10]

head(clinical3)








#### END ####
