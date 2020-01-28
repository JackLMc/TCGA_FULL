# A script to clean and process mutation data from the TCGA-COADREAD project
## Packages
library(UsefulFunctions)
library(tidyverse)
library(ggpubr)
library(maftools)

# Read in the public mutation data and merge with the clinical data ----
mutation_maf <- read.maf("./Data/Mutations/mc3.v0.2.8.PUBLIC.maf", verbose = F, isTCGA = T)

mutation <- muta@data # THESE ARE SOMATIC SNPs
mutation$IMPACT <- as.factor(mutation$IMPACT)
mutation$Patient.ID <- samptopat(mutation$Tumor_Sample_Barcode)
mutation$Patient.ID <- gsub("-", ".", mutation$Patient.ID)

# Collate the patients I want to look for (taken from 1_FPKM.R)
Patient_list <- read.delim("./Output/Patient_list.txt")
Patient_list <- levels(Patient_list$x)
tcga_mut <- mutation[mutation$Patient.ID %in% Patient_list, ]
tcga_mut <- factorthese(tcga_mut, c("Feature_type", "Feature", "BIOTYPE",
                                    "VARIANT_CLASS", "Variant_Classification",
                                    "Consequence", "Variant_Type", "Variant"))


POLE_muts <- droplevels(subset(tcga_mut, Hugo_Symbol == "POLE" | Hugo_Symbol == "POLD1"))
POLE_pats <- POLE_muts[duplicated(POLE_muts$Patient.ID), ]$Patient.ID %>% as.data.frame()
write.csv(POLE_pats, file = "./Output/POLE_mutants.csv", row.names = F)

KRAS_muts <- droplevels(subset(tcga_mut, Hugo_Symbol == "KRAS"))
KRAS_pats <- KRAS_muts[duplicated(KRAS_muts$Patient.ID), ]$Patient.ID %>% as.data.frame()
write.csv(KRAS_pats, file = "./Output/KRAS_mutants.csv", row.names = F)

rm(list = setdiff(ls(), c("mutation_maf", "tcga_mut")))
save.image("./R_Data/Mutation_clean.RData")

#### END ####



