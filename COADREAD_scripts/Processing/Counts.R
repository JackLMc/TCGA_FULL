library(UsefulFunctions)
library(tidyverse)
library(ggpubr)

my_comparisons <- list(c("MSS-hiCIRC", "MSI-H"),
                       c("MSS-hiCIRC", "MSS"),
                       c("MSI-H", "MSS"))
cbcols <- c("MSS-hiCIRC" = "#999999",
            "MSI-H" = "#56B4E9",
            "MSS" = "#009E73")

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


## Patients
pat_sub <- read.csv("./Output/Patient_Subtypes.csv")
converter <- read.delim("./Data/Important/gdc_sample_sheet_counts.tsv")
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



# this <- data.frame(stringsAsFactors = F)
# converter2$Patient.ID <- as.factor(converter2$Patient.ID)
# c <- 1
# for(i in levels(converter2$Patient.ID)){
#   work <- droplevels(subset(converter2, Patient.ID == i))
#   this[c, "Patient.ID"] <- i
#   this[c, "Sample_lev"] <- nlevels(work$Sample.ID)
#   this[c, "File_ID_lev"] <- nlevels(work$File.ID)
#   this[c, "File_Name_lev"] <- nlevels(work$File.Name)
#   c <- c + 1
#   }


temp_df1 <- merge(converter1[, c("Patient.ID", "File.ID")], temp_df, by = "File.ID")
library(reshape2)
Counts <- dcast(temp_df1, Gene ~ Patient.ID, sum, value.var = "Count")
# droplevels(subset(temp_df1, Patient.ID == "TCGA.A6.2672" & Gene == "ENSG00000000003.13")) # Matches the figure in "Try"

# library_sizes <- colSums(Counts[!'%in%'(colnames(Counts), "Gene")]) %>% as.data.frame() %>%
#   rownames_to_column(., var = "Patient.ID")
# colnames(library_sizes) <- c("Patient.ID", "Library_size")
# write.csv("./Output/Library_Sizes.csv", x = library_sizes, row.names = F)

Counts_cleaned <- droplevels(Counts[!'%in%'(Counts$Gene, 
                                            c("__alignment_not_unique", "__no_feature",
                                              "__not_aligned", "__too_low_aQual", "__ambiguous")), ])

save.image(file = "./R_Data/Counts.RData")


