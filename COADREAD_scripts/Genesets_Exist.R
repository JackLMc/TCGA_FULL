# A script to ensure that all genes in the Genesets I've chosen are in the Counts_cqn matrix
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



CTGenesets <- read.csv("./Exploratory_Data/Genesets/Investigative_Genesets.csv", stringsAsFactors = T)
SigGen <- read.csv("./Exploratory_Data/Genesets/Signature_Genesets.csv")


load("./R_Data/Counts_clean.RData")

pat_sub <- read.csv("./Output/Patient_Subtypes_09_03.csv")