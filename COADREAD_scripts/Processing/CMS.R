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
res <- CMScaller(Counts, RNAseq = T, doPlot = T)

CMS_groups <- rownames_to_column(res, var = "Patient.ID")
CMS_groups <- CMS_groups[, c("Patient.ID", "prediction")]
colnames(CMS_groups)[colnames(CMS_groups) == "prediction"] <- "CMS"

Patient_list <- read.csv("./Output/Patient_Subtypes_09_03.csv")$Patient.ID %>% as.character()
CMS_groups <- CMS_groups[CMS_groups$Patient.ID %in% Patient_list, ]

write.csv("./Output/CMS_groups.csv", x = CMS_groups, row.names = F)


pat_sub <- read.csv("./Output/Patient_Subtypes_09_03.csv")

colnames(Counts) %in% pat_sub$Patient.ID
Counts1 <- Counts[, colnames(Counts) %in% levels(pat_sub$Patient.ID)]
cam <- CMSgsa(emat=Counts1, class=pat_sub$Subtype, RNAseq=T)


CMS_groups <- read.csv("./Output/CMS_groups.csv")
Groups <- merge(CMS_groups, pat_sub, by = "Patient.ID")

library(reshape2)
dcast(Groups, CMS ~ Subtype)

CMS4 <- droplevels(subset(Groups, CMS == "CMS2" & Subtype != "MSI-H" | CMS == "CMS4" & Subtype !="MSI-H"))
dcast(CMS4, CMS ~ Subtype)

CMS4$Combined <- ifelse((CMS4$CMS == "CMS2" & CMS4$Subtype == "MSS-hiCIRC"), "CMS2hi", 
                        ifelse((CMS4$CMS == "CMS2" & CMS4$Subtype == "MSS"), "CMS2lo", 
                               ifelse((CMS4$CMS == "CMS4" & CMS4$Subtype == "MSS"), "CMS4lo", "CMS4hi")))



#### Load up the data ####
load("./R_Data/Counts_clean.RData")
library(tidyverse)
# install.packages(c("devtools", "curl"))
library(devtools)
# install_github("ebecht/MCPcounter", ref = "master", subdir = "Source")

library(MCPcounter)
library(UsefulFunctions)
library(ggpubr)

my_comparisons <- list(c("MSS-hiCIRC", "MSI-H"),
                       c("MSS-hiCIRC", "MSS"),
                       c("MSI-H", "MSS"))
cbcols <- c("MSS-hiCIRC" = "#999999",
            "MSI-H" = "#56B4E9",
            "MSS" = "#009E73",
            "MSI-L" = "#E69F00")

MCP_estimate <- MCPcounter.estimate(Counts_cqn, featuresType = "HUGO_symbols") %>% 
  as.data.frame() %>%
  rownames_to_column(., var = "Cell_Population") %>% 
  gather(contains("TCGA"), key = "Patient.ID", value = "Estimate")

MCP_estimate$Cell_Population <- as.factor(MCP_estimate$Cell_Population)



MCP <- merge(MCP_estimate, CMS4[, c("Patient.ID", "Combined")], by = "Patient.ID")
MCP <- factorthese(MCP, c("Patient.ID", "Combined"))


for(i in levels(MCP$Cell_Population)){
  print(i)
  work <- droplevels(subset(MCP, Cell_Population == i))
  temp_plot <- ggplot(work, aes(x = Combined, y = Estimate)) +
    geom_boxplot(alpha = 0.5, width = 0.2) + 
    geom_violin(aes(Combined, fill = Combined),
                scale = "width", alpha = 0.8) +
    scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00","#009E73" )) +
    labs(x = "Subtype", y = paste(i, " - MCP Estimate")) +
    theme_bw() +
    theme(axis.text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.direction = "horizontal", legend.position = "top") + 
    stat_compare_means(comparisons = list(c("CMS2hi", "CMS2lo"),
                                          c("CMS2hi", "CMS4hi"),
                                          c("CMS2lo", "CMS4lo"),
                                          c("CMS4hi", "CMS4lo")),
                       label = "p.signif", method = "wilcox.test")
  filen <- paste0(i, ".pdf")
  ggplot2:: ggsave(filen, plot = temp_plot, device = "pdf",
                   path = "./Figures/Gene_Sets/Enrichment/MCPcounter/CMS/All",
                   height = 6, width = 6)}














