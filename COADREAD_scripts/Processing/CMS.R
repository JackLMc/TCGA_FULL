# A script to call the CMS types of the colorectal cancer TCGA cohort
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

# source("Clinical.R") # Run to gain the clinical dataframe that's in Output (Clin_540)

load("./R_Data/Counts_clean1.RData")
# BiocManager::install("Biobase")
# devtools::install_github("Lothelab/CMScaller")
library(Biobase)
library(CMScaller)
library(biomaRt)
library(UsefulFunctions)
library(tidyverse)
library(ggpubr)
library(reshape2)

ensembl_DB <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# find1 <- listAttributes(ensembl_DB)$name %>% grepl(pattern = "entrez")
# listAttributes(ensembl_DB)[find1, ]

Gene_Map <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                  filters = "ensembl_gene_id", values = Counts_cleaned$Gene, ensembl_DB) # Only finds 

colnames(Counts_cleaned)[colnames(Counts_cleaned) == "Gene"] <- "ensembl_gene_id"

Counts_merged <- merge(Counts_cleaned, Gene_Map, by = "ensembl_gene_id")  
Counts_ENTREZ <- Counts_merged[, !'%in%'(colnames(Counts_merged), "ensembl_gene_id")]

Counts_long <- Counts_ENTREZ %>% 
  gather(-entrezgene_id, value = "Count", key = "Patient.ID")

Counts_totalled <- dcast(Counts_long, entrezgene_id ~ Patient.ID, sum, value.var = "Count")

Counts <- Counts_totalled[!is.na(Counts_totalled$entrezgene_id), ]
Counts$entrezgene_id <- as.factor(Counts$entrezgene_id)

Counts <- droplevels(subset(Counts, entrezgene_id != "")) # Remove genes which don't have a entrezgene_id
row.names(Counts) <- NULL
Counts <- column_to_rownames(Counts, var = "entrezgene_id")


### CMS prediction of TCGA primary colorectal cancers
## Need just counts, not cqn normalised counts...

res <- CMScaller(Counts, RNAseq=TRUE, doPlot=TRUE)
head(res)

CMS_groups <- rownames_to_column(res, var = "Patient.ID")
CMS_groups <- CMS_groups[, c("Patient.ID", "prediction")]
colnames(CMS_groups)[colnames(CMS_groups) == "prediction"] <- "CMS"
head(CMS_groups)

Patient_list <- read.csv("./Output/Patient_Subtypes_09_03.csv")$Patient.ID %>% as.character()
CMS_groups <- CMS_groups[CMS_groups$Patient.ID %in% Patient_list, ]

write.table("./Output/CMS_groups.txt", x = CMS_groups, quote = F, row.names = F)


pat_sub <- read.csv("./Output/Patient_Subtypes_09_03.csv")

### Camera Gene Set Analysis with CMS informative gene sets

colnames(Counts) %in% pat_sub$Patient.ID

Counts1 <- Counts[, colnames(Counts) %in% levels(pat_sub$Patient.ID)]
cam <- CMSgsa(emat=Counts1, class=pat_sub$Subtype, RNAseq=T)



