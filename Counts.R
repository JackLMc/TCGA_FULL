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
thousand.folders <- list.dirs(path = "./Data/Counts", full.names = T)
filelist1 <- sapply(thousand.folders[-1], function(x){
  list.files(x, pattern = "htseq.counts.gz$", full.names = T)})
filelist = unlist(filelist1)

# Read in files and combine
lists <- lapply(filelist, read.delim, header = F)
listsDF <- lists
lists <- listsDF

for (i in names(lists)){
  p <- gsub("./Data/Counts/", "", i)
  colnames(lists[[i]]) <- c("Gene", p)
}

names(lists) <- gsub("./Data/Counts/", "", names(lists))


# Patients I have CIRC scores for.
pat_sub <- read.csv("./Data/patient_subtypes.csv")
converter <- read.delim("./Data/Sample_Map.tsv")
converter$Patient.ID <- gsub("-", ".", converter$Case.ID)
converter1 <- converter[converter$Patient.ID %in% pat_sub$Patient.ID, ]
converter1$File.Name <- gsub(".htseq.counts.gz", "", converter1$File.Name)
# clin <- read.delim("./Data/Clinical/clinical.tsv")
# clin$Patient.ID <- gsub("-", ".", clin$submitter_id)
# clin1 <- clin[clin$Patient.ID %in% pat_sub$Patient.ID, ]

lists1 <- lists[names(lists) %in% converter1$File.ID]

multi_join <- function(list_of_loaded_data, join_func, ...){
  require("dplyr")
  output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
  return(output)
}

combined_df <- multi_join(lists1, full_join)

temp_df <- combined_df %>% gather(key = "File.ID", value = "Count", -Gene)
converter2 <- droplevels(subset(converter1, Sample.Type != "Solid Tissue Normal" &
                                  Sample.Type != "Blood Derived Normal"))


levels(converter2$Sample.Type)
levels(converter2$Sample.ID)

head(converter2)
this <- data.frame(stringsAsFactors = F)

converter2$Patient.ID <- as.factor(converter2$Patient.ID)
c <- 1
for(i in levels(converter2$Patient.ID)){
  work <- droplevels(subset(converter2, Patient.ID == i))
  this[c, "Patient.ID"] <- i
  this[c, "Sample_lev"] <- nlevels(work$Sample.ID)
  this[c, "File_ID_lev"] <- nlevels(work$File.ID)
  this[c, "File_Name_lev"] <- nlevels(work$File.Name)
  c <- c + 1
  }

temp_df1 <- merge(converter2[, c("Patient.ID", "File.ID")], temp_df, by = "File.ID")
Counts <- temp_df1[, c("Patient.ID", "Gene", "Count")] %>% spread(key = "Gene", value = "Count")
library(reshape2)
try <- dcast(temp_df1, Patient.ID ~ Gene, sum, value.var = "Count")

View(head(try))

?dcast
head(temp_df1)


for(i in levels(temp_df1$Patient.ID)
    
    ){
  work <- droplevels(subset(temp_df1))
}

