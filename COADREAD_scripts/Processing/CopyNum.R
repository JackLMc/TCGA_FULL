# A script to investigate copy number changes across Subtypes
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
# converter <- read.delim("./Data/CopyNum/sample.tsv")[, c("sample_id", "sample_submitter_id", "case_id", "case_submitter_id")]
# 
# converter$Patient.ID <- gsub("-", ".", converter$case_submitter_id)
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

combined_df <- multi_join(lists, full_join, by = c("Gene Symbol", "Gene ID", "Cytoband"))

colnames(combined_df)[colnames(combined_df) == "Gene Symbol"] <- "Gene.Symbol"
colnames(combined_df)[colnames(combined_df) == "Gene ID"] <- "Gene.ID"

##  
COADREAD <- read.delim("./Data/CopyNum/biospecimen.cart.2020-03-12/aliquot.tsv")#[, c("sample_submitter_id", "aliquot_id")]

# converter <- read.delim("./Data/Important/gdc_sample_sheet_FPKM.tsv")[, c("Case.ID", "Sample.ID", "Sample.Type")]
# converter1 <- droplevels(subset(converter, Sample.Type != "Solid Tissue Normal" &
#                                   Sample.Type != "Blood Derived Normal"))
# 
# COADREAD1 <- COADREAD[COADREAD$sample_submitter_id %in% converter1$Sample.ID, ] %>% droplevels()

# issue is there are multiple aliquot_ids for each sample ID

temp_df <- combined_df %>% gather(key = "aliquot_id", value = "CopyNum", -Gene.Symbol, -Gene.ID, -Cytoband)
temp_df1 <- merge(COADREAD[, c("aliquot_id", "sample_submitter_id")], temp_df, by = "aliquot_id")

# sumbitter id ends in 01 = Primary tumour
temp_df2 <- separate(temp_df1, col = "sample_submitter_id", into = c("this", "that", "there", "Sample.Type.Code"), sep = "-")
temp_df2$Patient.ID <- paste(temp_df2$this, temp_df2$that, temp_df2$there, sep = ".")
temp_df3 <- droplevels(subset(temp_df2, Sample.Type.Code == "01A" | Sample.Type.Code == "01B"))

# Take only the patients I have CIRC scores for
pat_sub <- read.csv("./Output/Patient_Subtypes_09_03.csv")
temp_df4 <- temp_df3[temp_df3$Patient.ID %in% pat_sub$Patient.ID, ] %>% droplevels()


# try <- temp_df4[!duplicated(temp_df4), ]
temp_df4$Patient.ID <- as.factor(temp_df4$Patient.ID)
temp_df4$Sample.Type.Code <- as.factor(temp_df4$Sample.Type.Code)




output <- data.frame(stringsAsFactors = F)
c <- 1
for(i in levels(temp_df4$Patient.ID)){
  print(i)
  work <- droplevels(subset(temp_df4, Patient.ID == i))
  output[c, "sample"] <- i
  output[c, "num_aliquo"] <- nlevels(work$aliquot_id)
  c <- c + 1
}

patients_with_two <- output[output$num_aliquo != 1, ]$sample
head(patients_with_two)



temp_df3 <- temp_df2[, c("Patient.ID", "Sample.Type.Code", "Gene.Symbol", "Gene.ID", "Cytoband", "CopyNum")]


## Choose which aliquot to use
choosing <- temp_df5[temp_df5$Patient.ID %in% patients_with_two, ] %>% droplevels()
choosing$Patient.ID <- as.factor(choosing$Patient.ID)
choosing$File.ID <- paste(choosing$Patient.ID, choosing$Sample.Type.Code, sep = "-")



data2 <- choosing[duplicated(choosing), ]


data1 <- choosing[, c("File.ID", "Gene.Symbol", "CopyNum")] %>% droplevels()

head(data2)


data2 <- spread(data1, key = "Gene.Symbol", value = "CopyNum")




output <- data.frame()
c <- 1
for(i in levels(choosing$Patient.ID)){
  print(i)
  one_pat <- droplevels(subset(choosing, Patient.ID == i))
  output[c, "Patient.ID"] <- i
  output[c, "Percentage_Match"] <- sum(one_pat[1, ] == one_pat[2, ])
  c <- c + 1
}
one_pat[1, ]

head(output)




temp_df2 <- temp_df1[!duplicated(temp_df1), ]



# Check all patients are 01A
length(temp_df1$sample_submitter_id[grepl("10A", temp_df1$sample_submitter_id)])
# length(temp_df1$sample_submitter_id[grepl("10B", temp_df1$sample_submitter_id)])
# length(temp_df1$sample_submitter_id[grepl("11A", temp_df1$sample_submitter_id)])
# length(temp_df1$sample_submitter_id[grepl("11B", temp_df1$sample_submitter_id)])
length(temp_df1$sample_submitter_id[grepl("01A", temp_df1$sample_submitter_id)])



temp_df2 <- temp_df1[, !'%in%'(colnames(temp_df1), c("aliquot_id", "sample_submitter_id", "Gene.ID"))]

library(reshape2)
head(temp_df2)
COPY <- spread(temp_df2, key = "Patient.ID", value = "CopyNum")

## Should these be used to talk about library sizes?
# Patient_list <- colnames(COPY_cleaned)[!'%in%'(colnames(COPY_cleaned), "Gene")]
# write.table(Patient_list, "./Output/Patient_list.txt", row.names = F)

# Gaining second dataframe (Symbols)
# biocLite("Homo.sapiens", dependencies = T)
library(Homo.sapiens)
COPY$Gene <- gsub("\\..*", "", COPY$Gene) #Removes version from Ensembl gene ID
geneid <- COPY$Gene

genes <- select(Homo.sapiens, keys = geneid, columns = c("SYMBOL", "TXCHROM", "ENTREZID"), 
                keytype = "ENSEMBL")
genes <- genes[!duplicated(genes$SYMBOL),]

colnames(COPY)[colnames(COPY) %in% "Gene"] <- "ENSEMBL"
FPKM <- merge(genes[, c("ENSEMBL", "SYMBOL")], COPY, by = "ENSEMBL") %>% 
  gather(contains("TCGA"), key = "Patient.ID", value = "FPKM")
library(reshape2)
FPKM1 <- dcast(FPKM, SYMBOL ~ Patient.ID, sum, value.var = "FPKM")



tail(temp_df)
