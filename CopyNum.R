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


# Looking for the files
thousand.folders <- list.dirs(path = "./Data/CopyNum", full.names = T)
filelist1 <- sapply(thousand.folders[-1], function(x){
  list.files(x, pattern = ".focal_score_by_genes.txt$", full.names = T)})
filelist = unlist(filelist1)

# Read in files and combine
lists <- lapply(filelist, read.delim, header = T, check.names = F)

# for (i in names(lists)){
#   p <- gsub("./Data/CopyNum/", "", i)
#   colnames(lists[[i]]) <- c("Gene", p)
# }

names(lists) <- gsub("./Data/CopyNum/", "", names(lists))



# Patients I have CIRC scores and microsatellite status for.
# pat_sub <- read.csv("./Data/Important/patient_subtypes.csv")
converter <- read.delim("./Data/CopyNum/sample.tsv")[, c("sample_id", "sample_submitter_id", "case_id", "case_submitter_id")]

converter$Patient.ID <- gsub("-", ".", converter$case_submitter_id)
# converter1 <- converter[converter$Patient.ID %in% pat_sub$Patient.ID, ]

# clin <- read.delim("./Data/Clinical/clinical.tsv")
# clin$Patient.ID <- gsub("-", ".", clin$submitter_id)
# clin1 <- clin[clin$Patient.ID %in% pat_sub$Patient.ID, ]

# lists1 <- lists[names(lists) %in% converter$File.ID]

multi_join <- function(list_of_loaded_data, join_func, ...){
  require("dplyr")
  output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
  return(output)
}

combined_df <- multi_join(lists, full_join)

colnames(combined_df)[colnames(combined_df) == "Gene Symbol"] <- "Gene.Symbol"
colnames(combined_df)[colnames(combined_df) == "Gene ID"] <- "Gene.ID"

GDC <- read.delim("./Data/GDC_large_mapping_TCGA.txt")

levels(GDC$data_type)


temp_df <- combined_df %>% gather(key = "File.ID", value = "CopyNum", -Gene.Symbol, -Gene.ID, -Cytoband)

temp_df[grepl("31a3877d", temp_df$File.ID), ]

head(temp_df)

converter$File.ID <- as.factor(converter$case_id)

clin <- read.delim("~/Downloads/clinical.cart.2019-03-28/clinical.tsv")[, c("case_id", "submitter_id")]
head(clin)

clin$File.ID <- clin$case_id
clin$Patient.ID <- gsub("-", ".", clin$submitter_id)

temp_df1 <- merge(converter2[, c("Patient.ID", "File.ID")], temp_df, by = "File.ID")
library(reshape2)


tail(temp_df)
