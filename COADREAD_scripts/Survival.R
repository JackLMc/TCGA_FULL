# A script to investigate survival between patient subgroups
## Packages
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

source("./COADREAD_scripts/Processing/Clinical.R")

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
pat_sub <- read.csv("./Output/Patient_Subtypes_13_02.csv")


clinical_outcome <- merge(pat_sub[, c("Patient.ID", "Subtype")], clin_out, by = "Patient.ID")

# Survival
library(survival)
os.surv <- Surv(clinical_outcome$OS_MONTHS, clinical_outcome$died)
fit1 <- survfit(os.surv ~ Subtype, data = clinical_outcome)
library(survminer)
pdf("./Figures/Clinical/Survival.pdf")
ggsurvplot(fit1, pval = T, pvalmethod = T, palette = c("#56B4E9",  "#009E73", "#999999"),
                    risk.table = F)
dev.off()


#### END ####
