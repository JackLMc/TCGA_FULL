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
            "MSS" = "#009E73")

Clin542 <- read.csv("./Output/Clinical_Data_542.csv")
pat_sub <- read.csv("./Output/Patient_Subtypes_13_02.csv")

WD <- merge(pat_sub, Clin542, by = "Patient.ID")




head(WD)
library(reshape2)
dcast(WD, MSI_STATUS ~ Pathologic_stage)

WD$Age_days <- abs(WD$Days_to_birth)
WD$Age_years <- abs(WD$Age_days/365.25)
WD$Age_years <- floor(WD$Age_years)

WD1 <- WD

# Collate the sites
WD1$Collated_site <- ifelse((is.na(WD1$Neoplasm_subdivision)), as.character(WD1$Site_of_disease), as.character(WD1$Neoplasm_subdivision))

## Rename to cecum, ascending, hepatic and transverse to right sided, all others to left - Mik et al., 2017 paper
WD1$Collated_site <- ifelse((WD1$Collated_site == "cecum" | WD1$Collated_site == "ascending colon" |
                               WD1$Collated_site == "hepatic flexure" | WD1$Collated_site == "transverse colon"), "right",
                                          ifelse((is.na(WD1$Collated_site)), NA, "left"))

WD1$Histological_Type <- ifelse((grepl("mucinous", WD1$Histological_Type)), "mucinous adenocarcinoma", "adenocarcinoma")

WD1$Pathologic_stage <- ifelse((WD1$Pathologic_stage == "stage i" | WD1$Pathologic_stage == "stage ia"), "1",
                               ifelse((WD1$Pathologic_stage == "stage ii" | WD1$Pathologic_stage == "stage iia" | WD1$Pathologic_stage == "stage iib"), "2",
                                      ifelse((WD1$Pathologic_stage == "stage iii" | WD1$Pathologic_stage == "stage iiia" |
                                                WD1$Pathologic_stage == "stage iiib" | WD1$Pathologic_stage == "stage iiic"), "3", 
                                             ifelse((WD1$Pathologic_stage == "<NA>"), NA, "4")))) %>% as.integer()

WD1$Pathologic_stage

range(WD1$Age_years, na.rm = T)
iR1 <- factorthese(WD1, colnames(WD1)[!'%in%'(colnames(WD1), c("Age_years", "Lymph_count", "Pathologic_stage"))])
str(iR1)
# Set MSS as the baseline
iR1$Subtype <- relevel(iR1$Subtype, "MSS")

iR2 <- column_to_rownames(iR1, var = "Patient.ID")

library(nnet)
test <- multinom(Subtype ~ 
                   Gender + Histological_Type + Pathologic_stage +
                   Collated_site +
                   Age_years + Venous_invasion + Lymph_count + Polyps, data = iR2)



z <- summary(test)$coefficients/summary(test)$standard.errors

p <- (1 - pnorm(abs(z), 0, 1))*2
p


dcast(WD1, Subtype ~ Pathologic_stage)
dcast(WD1, Subtype ~ Collated_site)
dcast(WD1, Subtype ~ Venous_invasion)
dcast(WD1, Subtype ~ Gender)




View(head(p)[, 1:2])




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



here <- merge(pat_sub, Clin542, by = "Patient.ID")
here$height_sq <- (here$Height/100) * (here$Height/100)
here$bmi <- here$Weight / here$height_sq

Bmi_only <- here[, c("Patient.ID", "Subtype", "bmi")]
 BMI <- na.omit(Bmi_only)
BMI <-droplevels(subset(BMI, bmi < 100))

View(BMI)

droplevels(subset(here, bmi > 40))


head(here)
library(tidyverse)
library(ggpubr)
ggplot(BMI, aes(x = Subtype, y = bmi)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(aes(Subtype, fill = Subtype),
              scale = "width", alpha = 0.8) +
  scale_fill_manual(values = cbcols) +
  labs(x = "MSI Status", y = "BMI") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top") + 
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test")

