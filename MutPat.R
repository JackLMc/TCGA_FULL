library(UsefulFunctions)
library(tidyverse)
library(ggpubr)

my_comparisons <- list(c("MSS-hiCIRC", "MSI-H"),
                       c("MSS-hiCIRC", "MSS"),
                       c("MSI-H", "MSS"))
cbcols <- c("MSS-hiCIRC" = "#999999",
            "MSI-H" = "#56B4E9",
            "MSS" = "#009E73")


library(maftools)
muta <- read.maf("./Data/Mutations/mc3.v0.2.8.PUBLIC.maf")

mutation <- muta@data # THIS IS SOMATIC SNPs


head(mutation)


