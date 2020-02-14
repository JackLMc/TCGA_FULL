# A script to investigate recurrent mutations across the patient subgroups from the TCGA-COADREAD project
## Packages
library(UsefulFunctions)
library(tidyverse)
library(ggpubr)
library(maftools)

my_comparisons <- list(c("MSS-hiCIRC", "MSI-H"),
                       c("MSS-hiCIRC", "MSS"),
                       c("MSI-H", "MSS"))
cbcols <- c("MSS-hiCIRC" = "#999999",
            "MSI-H" = "#56B4E9",
            "MSS" = "#009E73")

load("./R_Data/Mutation_clean.RData")

# Recurrent Mutations ----
devtools::install_github(repo = "PoisonAlien/TCGAmutations")

pat_sub <- read.csv("Output/Patient_Subtypes_13_02.csv")
pat_sub$Patient.ID <- gsub("\\.", "-", pat_sub$Patient.ID) %>% as.factor()

pat_sub1 <- droplevels(subset(pat_sub, Subtype == "MSS-hiCIRC"))
pat_sub2 <- droplevels(subset(pat_sub, Subtype == "MSS"))
pat_sub3 <- droplevels(subset(pat_sub, Subtype == "MSI-H"))

no_pats1 <- round(0.5 * nlevels(pat_sub1$Patient.ID)) # Apply this to the gene plot later - only take mutations present in half of patients
no_pats2 <- round(0.5 * nlevels(pat_sub2$Patient.ID))
no_pats3 <- round(0.5 * nlevels(pat_sub3$Patient.ID))

patients1 <- levels(pat_sub1$Patient.ID)
patients2 <- levels(pat_sub2$Patient.ID)
patients3 <- levels(pat_sub3$Patient.ID)

mutation_maf1 <- subsetMaf(mutation_maf, tsb = patients1, isTCGA = T, mafObj = T)
mutation_maf2 <- subsetMaf(mutation_maf, tsb = patients2, isTCGA = T, mafObj = T)
mutation_maf3 <- subsetMaf(mutation_maf, tsb = patients3, isTCGA = T, mafObj = T)

dev.off()
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
set.seed(123)
pdf("./Figures/Mutation/Recurrent/hiCIRC_muts.pdf")
geneCloud(mutation_maf1, minMut = no_pats1, random.order = F)
dev.off()

pdf("./Figures/Mutation/Recurrent/MSS_muts.pdf")
geneCloud(mutation_maf2, minMut = no_pats2, random.order = F)
dev.off()

pdf("./Figures/Mutation/Recurrent/MSI_muts.pdf")
geneCloud(mutation_maf3, minMut = no_pats3, random.order = F)
dev.off()

# save.image("./R_Data/Mutations.RData")
# load("./R_Data/Mutations.RData")

## Pick from geneCloud plot - clustered mutations?
lollipopPlot(mutation_maf1, gene = "APC")
lollipopPlot(mutation_maf2, gene = "APC")
lollipopPlot(mutation_maf3, gene = "APC")
dev.off()

lollipopPlot(mutation_maf1, gene = "TP53")
lollipopPlot(mutation_maf2, gene = "TP53")
lollipopPlot(mutation_maf3, gene = "TP53")

## Oncoplots
oncoplot(maf = mutation_maf1, top = 10)
dev.off()

oncoplot(maf = mutation_maf2, top = 10)
dev.off()

oncoplot(maf = mutation_maf1, top = 10)


# # All pretty much equal across subtypes
# # transitions and transversions
# laml.titv = titv(maf = mutation_maf1, plot = FALSE, useSyn = TRUE)
# #plot titv summary
# plotTiTv(res = laml.titv)
# 
# laml.titv = titv(maf = mutation_maf2, plot = FALSE, useSyn = TRUE)
# #plot titv summary
# plotTiTv(res = laml.titv)
# 
# laml.titv = titv(maf = mutation_maf3, plot = FALSE, useSyn = TRUE)
# #plot titv summary
# plotTiTv(res = laml.titv)
# dev.off()

somaticInteractions(maf = mutation_maf1, top = 10, pvalue = c(0.05, 0.1))
dev.off()
somaticInteractions(maf = mutation_maf2, top = 10, pvalue = c(0.05, 0.1))
dev.off()
somaticInteractions(maf = mutation_maf3, top = 10, pvalue = c(0.05, 0.1))
dev.off()



## Where do the recurrent mutations hit in the patient groups?
clinical_mutation <- merge(tcga_mut, pat_sub[, c("Patient.ID", "Subtype")], by = "Patient.ID")
TP53 <- droplevels(subset(clinical_mutation, Hugo_Symbol == "TP53"))
droplevels(subset(pat_sub, Subtype == "MSS-hiCIRC"))$Patient.ID

# 58/96 (have mutations in TP53)


TP53 %>% 
  group_by(Subtype) %>%
  summarise(no_pats = length(unique(Patient.ID)))
clinical_mutation  %>% 
  group_by(Subtype) %>%
  summarise(no_pats = length(unique(Patient.ID)))

hiCIRCs <- droplevels(subset(clinical_mutation, Hugo_Symbol == "TP53" & Subtype == "MSS-hiCIRC"))
hiCIRCs %>% 
  group_by(Exon_Number) %>%
  summarise(no_pats = length(unique(Patient.ID)))

MSS <- droplevels(subset(clinical_mutation, Hugo_Symbol == "TP53" & Subtype == "MSS"))
MSS %>% 
  group_by(Exon_Number) %>%
  summarise(no_pats = length(unique(Patient.ID)))


APC <- droplevels(subset(clinical_mutation, Hugo_Symbol == "APC"))

APC %>% 
  group_by(Subtype) %>%
  summarise(no_pats = length(unique(Patient.ID)))
clinical_mutation  %>% 
  group_by(Subtype) %>%
  summarise(no_pats = length(unique(Patient.ID)))

hiCIRCs <- droplevels(subset(clinical_mutation, Hugo_Symbol == "APC" & Subtype == "MSS-hiCIRC"))
hiCIRCs %>% 
  group_by(IMPACT) %>%
  summarise(no_pats = length(unique(Patient.ID)))

MSS <- droplevels(subset(clinical_mutation, Hugo_Symbol == "APC" & Subtype == "MSS"))
MSS %>% 
  group_by(IMPACT) %>%
  summarise(no_pats = length(unique(Patient.ID)))
