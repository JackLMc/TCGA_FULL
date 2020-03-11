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
pat_sub <- read.csv("./Output/Patient_Subtypes_09_03.csv")
WD <- merge(pat_sub, Clin542, by = "Patient.ID")

library(reshape2)

# Remove columns where NA values make up half of patients
WD1 <- WD[, colSums(is.na(WD)) < 300]

## Rename to cecum, ascending, hepatic and transverse to right sided, all others to left - Mik et al., 2017 paper
# levels(WD1$Patient.Primary.Tumor.Site)
WD1$Cancer_Sided <- ifelse((WD1$Patient.Primary.Tumor.Site == "Cecum" | WD1$Patient.Primary.Tumor.Site == "Ascending Colon" |
                               WD1$Patient.Primary.Tumor.Site == "Hepatic Flexure" | WD1$Patient.Primary.Tumor.Site == "Transverse Colon"), "Right",
                                          ifelse((is.na(WD1$Patient.Primary.Tumor.Site)), NA, "Left"))

# levels(WD1$Cancer.Type.Detailed)
WD1$Histological_Type <- ifelse((grepl("Mucinous", WD1$Cancer.Type.Detailed)), "Mucinous Adenocarcinoma", "Adenocarcinoma")



# levels(WD1$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code)
WD1$Stage <- ifelse((WD1$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code == "Stage I" | WD1$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code == "Stage IA"), "Stage I",
                               ifelse((WD1$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code == "Stage II" | WD1$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code == "Stage IIA" | WD1$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code == "Stage IIB"), "Stage II",
                                      ifelse((WD1$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code == "Stage III" | WD1$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code == "Stage IIIA" |
                                                WD1$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code == "Stage IIIB" | WD1$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code == "Stage IIIC"), "Stage III", 
                                             ifelse((WD1$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code == "<NA>"), NA, "Stage IV"))))


range(WD1$Diagnosis.Age, na.rm = T)

## Collate clinical data on MSI_STATUS
### Age
droplevels(subset(WD1, MSI_STATUS == "MSS"))$Diagnosis.Age %>% mean(., na.rm = T)
droplevels(subset(WD1, MSI_STATUS == "MSS"))$Diagnosis.Age %>% sd(., na.rm = T)

droplevels(subset(WD1, MSI_STATUS == "MSI-H"))$Diagnosis.Age %>% mean(., na.rm = T)
droplevels(subset(WD1, MSI_STATUS == "MSI-H"))$Diagnosis.Age %>% sd(., na.rm = T)

### BMI

WD1$Height_sq <- (WD1$Patient.Height/100) * (WD1$Patient.Height/100)
WD1$BMI <- WD1$Patient.Weight / WD1$Height_sq

droplevels(subset(WD1, MSI_STATUS == "MSS"))$BMI %>% mean(., na.rm = T)
droplevels(subset(WD1, MSI_STATUS == "MSS"))$BMI %>% sd(., na.rm = T)

droplevels(subset(WD1, MSI_STATUS == "MSI-H"))$BMI %>% mean(., na.rm = T)
droplevels(subset(WD1, MSI_STATUS == "MSI-H"))$BMI %>% sd(., na.rm = T)


### Gender
droplevels(subset(WD1, MSI_STATUS == "MSS")) %>% count(., Sex)
droplevels(subset(WD1, MSI_STATUS == "MSI-H")) %>% count(., Sex)

## Colon or Rectum
droplevels(subset(WD1, MSI_STATUS == "MSS")) %>% count(., Tumor.Tissue.Site)
droplevels(subset(WD1, MSI_STATUS == "MSI-H")) %>% count(., Tumor.Tissue.Site)

### Histological Type
droplevels(subset(WD1, MSI_STATUS == "MSS")) %>% count(., Histological_Type)
droplevels(subset(WD1, MSI_STATUS == "MSI-H")) %>% count(., Histological_Type)

### Sided
droplevels(subset(WD1, MSI_STATUS == "MSS")) %>% count(., Cancer_Sided)
droplevels(subset(WD1, MSI_STATUS == "MSI-H")) %>% count(., Cancer_Sided)


### Stage
droplevels(subset(WD1, MSI_STATUS == "MSS")) %>% count(., Stage)
droplevels(subset(WD1, MSI_STATUS == "MSI-H")) %>% count(., Stage)


### Venous Invasion
droplevels(subset(WD1, MSI_STATUS == "MSS")) %>% count(., Vascular.invasion.indicator)
droplevels(subset(WD1, MSI_STATUS == "MSI-H")) %>% count(., Vascular.invasion.indicator)




colnames(WD1)

# ## Tumour Tissue Source
# droplevels(subset(WD1, MSI_STATUS == "MSS")) %>% count(., Tissue.Source.Site)
# droplevels(subset(WD1, MSI_STATUS == "MSI-H")) %>% count(., Tissue.Source.Site)

### Lymph Count
droplevels(subset(WD1, MSI_STATUS == "MSS"))$Lymph.Node.s..Examined.Number %>% mean(., na.rm = T)
droplevels(subset(WD1, MSI_STATUS == "MSS"))$Lymph.Node.s..Examined.Number %>% sd(., na.rm = T)

droplevels(subset(WD1, MSI_STATUS == "MSI-H"))$Lymph.Node.s..Examined.Number %>% mean(., na.rm = T)
droplevels(subset(WD1, MSI_STATUS == "MSI-H"))$Lymph.Node.s..Examined.Number %>% sd(., na.rm = T)



## Check other bits
colnames(WD1)
droplevels(subset(WD1, MSI_STATUS == "MSS")) %>% count(., Lymphovascular.invasion.indicator)
droplevels(subset(WD1, MSI_STATUS == "MSI-H")) %>% count(., Lymphovascular.invasion.indicator)

Neoadjuvant.Therapy.Type.Administered.Prior.To.Resection.Text

###  Multinomial linear regression for Subtype
# Calculate BMI
multi <- WD1

Bmi_only <- multi[, c("Patient.ID", "Subtype", "BMI")] %>% na.omit
BMI <-droplevels(subset(Bmi_only, BMI < 100))

library(tidyverse)
library(ggpubr)
ggplot(BMI, aes(x = Subtype, y = BMI)) +
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


### Check multinomial linear regression
head(multi)

multi$BMI_group <- ifelse((multi$BMI < 18.5), "Underweight", 
                          ifelse((multi$BMI > 18.5 & multi$BMI < 24.9), "Healthy", 
                                 ifelse((multi$BMI > 24.9 & multi$BMI < 29.9), "Overweight", 
                                        ifelse((is.na(multi$BMI)), NA, "Obese"))))

head(multi)
iR1 <- factorthese(multi, colnames(multi)[!'%in%'(colnames(multi), c("Diagnosis.Age", "Lymph.Node.s..Examined.Number", "BMI"))])

# Set MSS as the baseline
iR1$Subtype <- relevel(iR1$Subtype, "MSS")
iR2 <- column_to_rownames(iR1, var = "Patient.ID")

library(nnet)
test <- multinom(Subtype ~ 
                   Sex + Histological_Type + Stage +
                   Cancer_Sided +
                   Diagnosis.Age + Vascular.invasion.indicator + Lymph.Node.s..Examined.Number + BMI_group + Patient.s.Vital.Status, data = iR2)

z <- summary(test)$coefficients/summary(test)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1))*2
p


dcast(multi, Subtype ~ Stage)
dcast(multi, Subtype ~ Cancer_Sided)
dcast(multi, Subtype ~ Vascular.invasion.indicator)
dcast(multi, Subtype ~ Sex)
dcast(multi, Subtype ~ BMI_group)


# Survival
clinical_outcome <- multi[, c("Patient.ID", "Subtype", "Overall.Survival..Months.", "Overall.Survival.Status")]
clinical_outcome$died <- ifelse((clinical_outcome$Overall.Survival.Status == "LIVING"), F, 
                                ifelse((is.na(clinical_outcome$Overall.Survival.Status)), NA, T))

clinical_outcome <- clinical_outcome[!is.na(clinical_outcome$Overall.Survival..Months.), ]
clinical_outcome$Overall.Survival..Months. <- as.numeric(clinical_outcome$Overall.Survival..Months.)

library(survival)
os.surv <- Surv(clinical_outcome$Overall.Survival..Months., clinical_outcome$died)
fit1 <- survfit(os.surv ~ Subtype, data = clinical_outcome)
library(survminer)
pdf("./Figures/Clinical/Survival.pdf")
ggsurvplot(fit1, pval = T, pvalmethod = T, palette = c("#56B4E9",  "#009E73", "#999999"),
           risk.table = T)
dev.off()


#### END ####


