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
paper_data <- clin_app[, c("Patient.ID", "Microsatellite.instability")]
cbio_data <- pats_tcga_pub[, c("Patient.ID", "MSI_STATUS")]
colnames(paper_data) <- c("Patient.ID", "MSI_STATUS")
full_Data <- rbind(cbio_data, paper_data)

FD <- full_Data[!duplicated(full_Data), ]
FD1 <- FD[!is.na(FD$MSI_STATUS),] %>% droplevels()
FD1$MSI_STATUS <- ifelse((FD1$MSI_STATUS == "MSI-H"), "MSI-H", "MSS")


write.csv("./Output/Clinical_Data_614.csv", x = FD1, row.names = F)






pat_sub <- read.csv("./Output/Patient_Subtypes.csv")
ehh <- read.delim("./Data/clinical.cart.2019-03-15/exposure.tsv")
eh <- droplevels(subset(ehh, weight != "--"))
eh$Patient.ID <- gsub("-", ".", eh$submitter_id)


try <- merge(eh[, c("Patient.ID", "bmi", "weight", "height")], pat_sub)
try$weight <- as.numeric(as.character(try$weight))

try$bmi
try$bmi <- as.numeric(as.character(try$bmi))
try$height <- as.numeric(as.character(try$height))

try1 <- droplevels(subset(try, height != "--"))
try1$height <- as.numeric(as.character(try1$height))


ggplot(try1, aes(x = Subtype, y = height)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(aes(Subtype, fill = Subtype),
              scale = "width", alpha = 0.8) +
  scale_fill_manual(values = cbcols) +
  labs(x = "MSI Status", y = "weight") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top") + 
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test")





clinical_foll <- read.delim("./Data/clinical.cart.2019-03-15/clinical.tsv")
clinical_foll$Patient.ID <- gsub("-", ".", clinical_foll$submitter_id)


this <- merge(clinical_foll, pat_sub, by = "Patient.ID")

this1 <- droplevels(subset(this, age_at_diagnosis != "--"))
this1$age_at_diagnosis <- as.numeric(as.character(this1$age_at_diagnosis))

ggplot(this1, aes(x = Subtype, y = age_at_diagnosis)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(aes(Subtype, fill = Subtype),
              scale = "width", alpha = 0.8) +
  scale_fill_manual(values = cbcols) +
  labs(x = "MSI Status", y = "weight") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top") + 
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test")
