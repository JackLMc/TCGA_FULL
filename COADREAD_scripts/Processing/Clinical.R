library(UsefulFunctions)
library(tidyverse)
library(ggpubr)

my_comparisons <- list(c("MSS-hiCIRC", "MSI-H"),
                       c("MSS-hiCIRC", "MSS"),
                       c("MSI-H", "MSS"))
cbcols <- c("MSS-hiCIRC" = "#999999",
            "MSI-H" = "#56B4E9",
            "MSS" = "#009E73",
            "MSI-L" = "#E69F00")

# Collate the patients I want to look for
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
clin_app <- Paper_clinical[Paper_clinical$Patient.ID %in% Patient_list, ] # Finds 380 Patients
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
FD1$MSI_STATUS <- ifelse((FD1$MSI_STATUS == "MSI-H"), "MSI-H", "MSS")


write.csv("./Output/Clinical_Data_614.csv", x = FD1, row.names = F)

# Collated survival stats ----
clin_out_paper<- clin_app[, c("Patient.ID", "Overall.survival..days.", "Survival.event..1.death.")]
clin_out_cBio <- tcga_pub_clinical[, c("Patient.ID", "OS_MONTHS", "OS_STATUS")]

clin_out_paper$OS_MONTHS <- clin_out_paper$Overall.survival..days./12
clin_out_paper <- clin_out_paper[!is.na(clin_out_paper$Survival.event..1.death.), ]
clin_out_paper$died <- ifelse((clin_out_paper$Survival.event..1.death. == 1), T, F)

clin_out_paper <- clin_out_paper[, c("Patient.ID", "OS_MONTHS", "died")]


clin_out_cBio <- droplevels(subset(clin_out_cBio, OS_STATUS != "NC"))
clin_out_cBio$died <- ifelse((clin_out_cBio$OS_STATUS == "LIVING"), F, T)
clin_out_cBio <- clin_out_cBio[, c("Patient.ID", "OS_MONTHS", "died")]

clin_out <- rbind(clin_out_paper, clin_out_cBio)
pat_sub <- read.csv("./Output/Patient_Subtypes.csv")


clinical_outcome <- merge(pat_sub[, c("Patient.ID", "Subtype")], clin_out, by = "Patient.ID")

# Survival
library(survival)
os.surv <- Surv(clinical_outcome$OS_MONTHS, clinical_outcome$died)
fit1 <- survfit(os.surv ~ Subtype, data = clinical_outcome)
library(survminer)
survp <- ggsurvplot(fit1, pval = T, pvalmethod = T, palette = c("#56B4E9",  "#009E73", "#999999"),
                    risk.table = F)

#### end ####

