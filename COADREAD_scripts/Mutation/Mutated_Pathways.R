# A script to investigate mutated pathways in the TCGA-COADREAD project
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



#### Mutated pathways ####
if (!requireNamespace("BiocManager", quietly = T))
  install.packages("BiocManager")

BiocManager::install("seq2pathway")
devtools:: install_github("https://github.com/cran/WGCNA")
library(seq2pathway)


browseVignettes("seq2pathway")

