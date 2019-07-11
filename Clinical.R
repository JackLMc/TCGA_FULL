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


# ## Investigate from clinical carts
# pat_sub <- read.csv("./Output/Patient_Subtypes.csv")
# exposure <- read.delim("./Data/clinical.cart.2019-03-15/exposure.tsv")
# expo <- droplevels(subset(exposure, weight != "--"))
# expo$Patient.ID <- gsub("-", ".", expo$submitter_id)
# 
# 
# CC <- merge(expo[, c("Patient.ID", "bmi", "weight", "height")], pat_sub)
# CC$weight <- as.numeric(as.character(CC$weight))
# 
# CC$bmi
# CC$bmi <- as.numeric(as.character(CC$bmi))
# CC$height <- as.numeric(as.character(CC$height))
# 
# CC1 <- droplevels(subset(CC, height != "--"))
# CC1$height <- as.numeric(as.character(CC1$height))
# 
# 
# ggplot(CC1, aes(x = Subtype, y = height)) +
#   geom_boxplot(alpha = 0.5, width = 0.2) + 
#   geom_violin(aes(Subtype, fill = Subtype),
#               scale = "width", alpha = 0.8) +
#   scale_fill_manual(values = cbcols) +
#   labs(x = "MSI Status", y = "weight") +
#   theme_bw() +
#   theme(axis.text = element_text(size = 16)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme(legend.direction = "horizontal", legend.position = "top") + 
#   stat_compare_means(comparisons = my_comparisons,
#                      label = "p.signif", method = "wilcox.test")
# 
# 
# ## Survival stats - collated.
# Paper_clinical <- read.csv("./Data/Clinical_data/Paper_Clinical.csv")
# Paper_clinical$Patient.ID <- samptopat(Paper_clinical$TCGA.ID)
# Paper_clinical$Patient.ID <- gsub("-", ".", Paper_clinical$Patient.ID)
# 
# pat_sub <- read.csv("./Output/Patient_Subtypes.csv")
# 
# 
# survival <- merge(Paper_clinical, pat_sub, by = "Patient.ID")
# 
# surviv1 <- survival[!is.na(survival$Overall.survival..days.), ] %>% droplevels()
# surviv1 <- surviv1[, c("Patient.ID", "Overall.survival..days.", "Survival.event..1.death.", "Microsatellite.instability")]
# surviv1$OS_MONTHS <- surviv1$Overall.survival..days./30.42
# surviv1$OS_STATUS <- ifelse((surviv1$Survival.event..1.death. == 1), "DECEASED", "LIVING")
# surviv1 <- surviv1[, c("Patient.ID", "OS_MONTHS", "OS_STATUS", "Microsatellite.instability")]
# colnames(surviv1)[colnames(surviv1) == "Microsatellite.instability"] <- "MSI_STATUS"
# 
# tcga_pub <- tcga_pub_clinical[!is.na(tcga_pub_clinical$OS_MONTHS), ] %>% droplevels()
# tcga_pub <- tcga_pub[, c("Patient.ID", "OS_STATUS", "OS_MONTHS", "MSI_STATUS")]
# tcga_pub <- tcga_pub[!'%in%'(tcga_pub$Patient.ID, surviv1$Patient.ID),]
# 
# att <- rbind(tcga_pub, surviv1)
# att$MSI_STATUS <- ifelse((att$MSI_STATUS == "MSI-H"), "MSI-H", "MSS")
# clin6 <- att
# clin6 <- clin6[!duplicated(clin6),]
# 
# 
# clin6$died <- ifelse((clin6$OS_STATUS == "LIVING"), F, T)
# clin6$OS_YEARS <- clin6$OS_MONTHS / 12
# clin7 <- merge(clin6, pat_sub, by = "Patient.ID")
# 
# library(survival)
# clin7 <- clin7[!duplicated(clin7), ]
# os.surv <- Surv(clin7$OS_MONTHS, clin7$died)
# 
# 
# fit1 <- survfit(os.surv ~ MSI_STATUS, data = clin7)
# library(survminer)
# 
# survp <- ggsurvplot(fit1, pval = T, pvalmethod = T, palette = c("#56B4E9",  "#009E73", "#999999"), 
#                     risk.table = F)
# ggsave("./Figures/Clinical/survival_comb.pdf", plot = survp$plot, 
#        height = 6, width = 6)
# 
# library(reshape2)
# dcast(clin7, Subtype ~ OS_STATUS, length)
# 
# 
# # pat_sub$Patient.ID[!'%in%'(pat_sub$Patient.ID, att$Patient.ID)] %>% droplevels()
# 
# 
# head(survival)
# vir <- survival
# vir$Viral <- ifelse((is.na(vir$Virus.detection.)), "NC", as.character(vir$Virus.detection.))
# vir$Viral <- as.factor(vir$Viral)
# 
# dcast(Subtype ~ Viral, vir)
# 
# dcast(vir, Subtype ~ Viral, length)
# 
# head(vir)
# 
# # PANCAN stuff
# this <- read.delim("./Data/Clinical_data/clinical_PANCAN_patient_with_followup.tsv")
# this$Patient.ID <- gsub("-", ".", this$bcr_patient_barcode)
# 
# pat_sub <- read.csv("./Output/Patient_Subtypes.csv")
# that <- merge(this, pat_sub, by = "Patient.ID")
# 
# emptycols <- sapply(that, function (k) all(is.na(k)))
# that1 <- that[!emptycols]
# 
# head(that1)
# 
# View(colnames(that1))
# 
# str(that1$weight)
# 
# that2 <- that1[!is.na(that1$weight),]
# that2 <- droplevels(subset(that2, weight != "[Not Available]"))
# that2 <- droplevels(subset(that2, weight != "."))
# levels(that2$weight)
# 
# that2$weight <- as.numeric(as.character(that2$weight))
# 
# ggplot(that2, aes(x = Subtype, y = weight)) +
#   geom_boxplot(alpha = 0.5, width = 0.2) + 
#   geom_violin(aes(Subtype, fill = Subtype),
#               scale = "width", alpha = 0.8) +
#   scale_fill_manual(values = cbcols) +
#   labs(x = "MSI Status", y = "weight") +
#   theme_bw() +
#   theme(axis.text = element_text(size = 16)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme(legend.direction = "horizontal", legend.position = "top") + 
#   stat_compare_means(comparisons = my_comparisons,
#                      label = "p.signif", method = "wilcox.test")
# 
# 
# 
# emptycols <- sapply(that1, function (k) all(k == "[Not Available]"))
# that1 <- that1[!emptycols]
# emptycols <- sapply(that1, function (k) all(k == ""))
# that1 <- that1[!emptycols]
# 
# 
# View(colnames(that1))
# 
# dcast(that2, Subtype ~ history_of_neoadjuvant_treatment)
# 
# 
# that3 <- droplevels(subset(that2, preoperative_pretreatment_cea_level != "[Not Available]"))
# that3$preoperative_pretreatment_cea_level <- as.numeric(as.character(that3$preoperative_pretreatment_cea_level))
# 
# ggplot(that3, aes(x = Subtype, y = preoperative_pretreatment_cea_level)) +
#   geom_boxplot(alpha = 0.5, width = 0.2) + 
#   geom_violin(aes(Subtype, fill = Subtype),
#               scale = "width", alpha = 0.8) +
#   scale_fill_manual(values = cbcols) +
#   labs(x = "MSI Status", y = "weight") +
#   theme_bw() +
#   theme(axis.text = element_text(size = 16)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme(legend.direction = "horizontal", legend.position = "top") + 
#   stat_compare_means(comparisons = my_comparisons,
#                      label = "p.signif", method = "wilcox.test")
