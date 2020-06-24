# A script to determine mcp counter abundances of immune cells
# Load data and required packages
load("./R_Data/Counts_clean.RData")
library(tidyverse)
# install.packages(c("devtools", "curl"))
library(devtools)
install_github("ebecht/MCPcounter",ref="master", subdir="Source")

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



head(Counts_cqn)
## Estimate MCPcounter abundances
MCP_estimate <- MCPcounter.estimate(Counts_cqn, featuresType = "HUGO_symbols") %>% 
  as.data.frame() %>%
  rownames_to_column(., var = "Cell_Population") %>% 
  gather(contains("TCGA"), key = "Patient.ID", value = "Estimate")

MCP_estimate$Cell_Population <- as.factor(MCP_estimate$Cell_Population)
pat_sub <- read.csv("./Output/Patient_Subtypes_09_03.csv")

MCP <- merge(MCP_estimate, pat_sub[, c("Patient.ID", "Subtype")], by = "Patient.ID")
MCP <- factorthese(MCP, c("Patient.ID", "Subtype"))


for(i in levels(MCP$Cell_Population)){
  print(i)
  work <- droplevels(subset(MCP, Cell_Population == i))
  temp_plot <- ggplot(work, aes(x = Subtype, y = Estimate)) +
    geom_boxplot(alpha = 0.5, width = 0.2) + 
    geom_violin(aes(Subtype, fill = Subtype),
                scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) +
    labs(x = "Subtype", y = paste(i, " - MCP Estimate")) +
    theme_bw() +
    theme(axis.text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.direction = "horizontal", legend.position = "top") + 
    stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif", method = "wilcox.test")
  filen <- paste0(i, ".pdf")
  ggplot2:: ggsave(filen, plot = temp_plot, device = "pdf",
                   path = "./Figures/Gene_Sets/Enrichment/MCPcounter",
                   height = 6, width = 6)}



# Estimate my Th17 cells?
CTGenesets <- read.csv("./Exploratory_Data/Genesets/Castro_collated.csv")
colnames(CTGenesets) <- c("Cell population", "HUGO symbols")

head(CTGenesets)

MCP_estimate <- MCPcounter.estimate(Counts_cqn, featuresType = "HUGO_symbols", genes = CTGenesets) %>% 
  as.data.frame() %>%
  rownames_to_column(., var = "Cell_Population") %>% 
  gather(contains("TCGA"), key = "Patient.ID", value = "Estimate")

MCP <- merge(MCP_estimate, pat_sub[, c("Patient.ID", "Subtype")], by = "Patient.ID")
str(MCP)
MCP <- factorthese(MCP, c("Patient.ID", "Subtype", "Cell_Population"))


for(i in levels(MCP$Cell_Population)){
  print(i)
  work <- droplevels(subset(MCP, Cell_Population == i))
  temp_plot <- ggplot(work, aes(x = Subtype, y = Estimate)) +
    geom_boxplot(alpha = 0.5, width = 0.2) + 
    geom_violin(aes(Subtype, fill = Subtype),
                scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) +
    labs(x = "Subtype", y = paste(i, " - MCP Estimate")) +
    theme_bw() +
    theme(axis.text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.direction = "horizontal", legend.position = "top") + 
    stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif", method = "wilcox.test")
  filen <- paste0(i, ".pdf")
  ggplot2:: ggsave(filen, plot = temp_plot, device = "pdf",
                   path = "./Figures/Gene_Sets/Enrichment/MCPcounter/Bespoke",
                   height = 6, width = 6)}



