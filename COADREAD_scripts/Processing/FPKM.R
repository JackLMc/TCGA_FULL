# A script to clean and process the FPKM files from the COADREAD project on TCGA
## Packages
library(UsefulFunctions)
library(tidyverse)
library(ggpubr)

# Looking for the files
thousand.folders <- list.dirs(path = "./Data/FPKM", full.names = T)
filelist1 <- sapply(thousand.folders[-1], function(x){
  list.files(x, pattern = "FPKM.txt.gz$", full.names = T)})
filelist = unlist(filelist1)

# Read in files and combine
lists <- lapply(filelist, read.delim, header = F)
listsDF <- lists
lists <- listsDF

for (i in names(lists)){
  p <- gsub("./Data/FPKM/", "", i)
  colnames(lists[[i]]) <- c("Gene", p)}

names(lists) <- gsub("./Data/FPKM/", "", names(lists))

# GDC sample sheet
converter <- read.delim("./Data/Important/gdc_sample_sheet_FPKM.tsv")
converter$Patient.ID <- gsub("-", ".", converter$Case.ID)
# converter1 <- converter[converter$Patient.ID %in% pat_sub$Patient.ID, ]
converter$File.Name <- gsub(".FPKM.txt.gz", "", converter$File.Name)
# clin <- read.delim("./Data/Clinical/clinical.tsv")
# clin$Patient.ID <- gsub("-", ".", clin$submitter_id)
# clin1 <- clin[clin$Patient.ID %in% pat_sub$Patient.ID, ]

lists1 <- lists[names(lists) %in% converter$File.ID]
converter <- converter[converter$File.ID %in% names(lists), ]

multi_join <- function(list_of_loaded_data, join_func, ...){
  require("dplyr")
  output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
  return(output)
}

combined_df <- multi_join(lists1, full_join, by = "Gene")

temp_df <- combined_df %>% gather(key = "File.ID", value = "FPKM", -Gene)
converter2 <- droplevels(subset(converter,  Sample.Type == "Primary Tumor"))

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

temp_df1 <- merge(converter2[, c("Patient.ID", "File.ID")], temp_df, by = "File.ID")
library(reshape2)
dupes <- temp_df1[duplicated(temp_df1[, !'%in%'(colnames(temp_df1), c("FPKM"))]), ]

FPKMs <- dcast(temp_df1, Gene ~ Patient.ID, sum, value.var = "FPKM") # this does nothing to the data:
# nlevels(temp_df1$Gene) * nlevels(temp_df1$File.ID) 
# dim(temp_df1)


## Should these be used to talk about library sizes?
# Patient_list <- colnames(FPKMs)[!'%in%'(colnames(FPKMs), "Gene")]
# write.table(Patient_list, "./Output/Patient_list.txt", row.names = F)

# Gaining second dataframe (Symbols)
# biocLite("Homo.sapiens", dependencies = T)
library(Homo.sapiens)
FPKMs$Gene <- gsub("\\..*", "", FPKMs$Gene) #Removes version from Ensembl gene ID
geneid <- FPKMs$Gene
genes <- select(Homo.sapiens, keys = geneid, columns = c("SYMBOL", "TXCHROM", "ENTREZID"), 
                keytype = "ENSEMBL")
genes <- genes[!duplicated(genes$SYMBOL),]

colnames(FPKMs)[colnames(FPKMs) %in% "Gene"] <- "ENSEMBL"

FPKM <- merge(genes[, c("ENSEMBL", "SYMBOL")], FPKMs, by = "ENSEMBL") %>% 
  gather(contains("TCGA"), key = "Patient.ID", value = "FPKM")
library(reshape2)
FPKM1 <- dcast(FPKM, SYMBOL ~ Patient.ID, mean, value.var = "FPKM")
FPKM2a <- FPKM1[!is.na(FPKM1$SYMBOL), ]
rownames(FPKM2a) <- NULL
FPKM2 <- FPKM2a

BiocManager::install("cqn")

FPKM2[!'%in%'(colnames(FPKM2), c("SYMBOL"))] <- log2(FPKM2[!'%in%'(colnames(FPKM2), c("SYMBOL"))] + 1)
library(cqn)

?cqn

FPKM3 <- FPKM2 %>%
  column_to_rownames(., var = "SYMBOL") %>%
  as.matrix()




rm(list = setdiff(ls(), c("FPKM", "FPKMs", "FPKM1", "FPKM2", "FPKM2a", "FPKM3")))
save.image("./R_Data/FPKM_clean.RData")

write.table("./Output/Patient_list.txt", x = as.factor(colnames(FPKM3)), row.names = F)

#### END ####