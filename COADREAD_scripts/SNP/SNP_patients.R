# A R script to find the file IDs for the SNP enrichment analysis of the TCGA COADREAD project
library(tidyverse)
library(UsefulFunctions)

set.seed(123)

# install.packages("rjson")
library("rjson")

meta <- fromJSON(file = "./Data/SNPs/metadata.cart.2019-11-15.json")

# finding out what I need from the meta_data
# test <- as.data.frame(meta[[2]]$annotations)
# cbind(test, meta[[2]]$file_id)
# available_data <- list()
# i <- 1
# for(i in seq_along(meta)){
#   work <- meta[[i]]
#   print(work$file_id)
#   available_data[[work$file_id]] <- names(work)
#   }
# common <- Reduce(intersect, available_data)
# 
# common
# 
# meta[[1]]$archive

# Need to take out the associated_entities and file_id for each element of the meta list
## This gains patient ID which I can link to my groups
list_of_annots <- list()
for(i in seq_along(meta)){
  work <- meta[[i]]
  print(work$file_id)
  annot <- as.data.frame(work$associated_entities)
  annot1 <- cbind(annot, work$file_id)  %>% mutate_all(as.character)
  list_of_annots[[work$file_id]] <- annot1
}
meta_df <- bind_rows(list_of_annots, .id = "column_label")

required <- c("work$file_id", "entity_submitter_id")
meta_df1 <- meta_df[, colnames(meta_df) %in% required]
colnames(meta_df1)[colnames(meta_df1) == "work$file_id"] <- "file_id"

# meta_df1

# Need to find the sample types for each of the files, to avoid using primary tumours - I want germline SNPs, nothing from the primary tumour
list_of_cases <- list()
for(i in seq_along(meta)){
  work <- meta[[i]]
  print(work$file_id)
  cases <- work$cases
  list_of_cases[[work$file_id]] <- cases
}

# Extracting data from the cases
list_of_samples <- list()
for(i in names(list_of_cases)){
  work <- list_of_cases[[i]]
  print(i)
  l <- work[[1]]$samples[[1]]
  samp <- data.frame(as.list(unlist(l))) %>%
    mutate_all(as.character)
  rownames(samp) <- NULL
  samp$file_id <- i
  # samp1 <- select(samp, -contains("portion"))
  list_of_samples[[i]] <- samp
}

meta_sample <- bind_rows(list_of_samples)
meta_data_samples <- meta_sample[, c("file_id", "submitter_id", "sample_type")]

META <- merge(meta_df1, meta_data_samples, by = "file_id")
META$Patient.ID <- samptopat(META$submitter_id)
META$Patient.ID <- gsub("-", ".", META$Patient.ID)
META$file_id <- gsub("-", ".", META$file_id)
META <- factorthese(META, colnames(META))

normals <- droplevels(subset(META, sample_type == "Blood Derived Normal"))
hiCIRC <- read.csv("./Output/Patient_Subtypes_13_02.csv")

files <- merge(normals, hiCIRC, by = "Patient.ID") %>% droplevels()

dim(files)

CIRC <- droplevels(subset(files, Subtype == "MSS-hiCIRC"))
MSS <- droplevels(subset(files, Subtype == "MSS"))
MSI <- droplevels(subset(files, Subtype == "MSI-H"))

nlevels(CIRC$Patient.ID)
nlevels(MSS$Patient.ID)
nlevels(MSI$Patient.ID)

## Find any duplicates
files$Patient.ID <- as.factor(files$Patient.ID)
files$file_id <- as.factor(files$file_id)
Dupes <- data.frame()
c <- 1
for(i in levels(files$Patient.ID)){
  print(i)
  working <- droplevels(subset(files, Patient.ID == i))
  number <- nlevels(working$file_id)
  Dupes[c, "Patient.ID"] <- i
  Dupes[c, "Number"] <- number
  c <- c + 1
}

Pats_with_dupes <- droplevels(subset(Dupes, Number != 1))$Patient.ID

## Check whether the files for the duplicated patients are the same
SNP <- read.delim("./Data/SNPs/matrix_genotypes.final.txt")

SNP1 <- column_to_rownames(SNP, var = "snpID") %>%
  as.matrix() %>%
  t()
rownames(SNP1) <- gsub("^X", "", rownames(SNP1))

SNP_normals <- SNP1[rownames(SNP1) %in% normals$file_id, ]

file_id_for_dupes <- normals[normals$Patient.ID %in% Pats_with_dupes, ] %>% droplevels()
dupes_file <- file_id_for_dupes$file_id %>% as.character()

SNP_dupes <- SNP_normals[rownames(SNP_normals) %in% dupes_file, ]
dim(SNP_dupes)

file_id_for_dupes$Patient.ID <- as.factor(file_id_for_dupes$Patient.ID)

output <- data.frame()
c <- 1
for(i in levels(file_id_for_dupes$Patient.ID)){
  print(i)
  one_pat <- droplevels(subset(file_id_for_dupes, Patient.ID == i))$file_id %>% as.character()
  reduced_dupe <- SNP_dupes[rownames(SNP_dupes) %in% one_pat, ]
  output[c, "Patient.ID"] <- i
  output[c, "Percentage_Match"] <- sum(reduced_dupe[1, ] == reduced_dupe[2, ])/906600
  c <- c + 1
}

mean(output$Percentage_Match)

# The patients with more than one file in the SNP array do not have concordant results, later date. Boris said it shouldn't matter too much.


# Some patients have duplications, take first instance of these patients
undup_hiCIRC <- files[!duplicated(files$Patient.ID),]$file_id %>% as.character
working_SNPs <- SNP_normals[rownames(SNP_normals) %in% undup_hiCIRC, ]

df1 <- rownames_to_column(as.data.frame(working_SNPs), var = "file_id")
df2 <- merge(df1, files[, c("file_id", "Patient.ID")], by = "file_id") %>% column_to_rownames(., var = "Patient.ID")

head(df2)[, 1:6]
this_bit <- rownames_to_column(df2, var = "Patient.ID") %>% .[, 1:2]
TB <- merge(this_bit, hiCIRC, by = "Patient.ID")[, c("Patient.ID", "file_id", "Subtype")]

write.csv(TB, file = "./Output/file_id_pats_SNP.csv", row.names = F)

clusters <- TB[, c("file_id", "Subtype")]
clusters$file_id <- gsub("\\.", "-", clusters$file_id)
write.table("./Output/clusters_13_02.txt", x = clusters, row.names = F, quote = F, sep = "\t")

clusters_MSS <- droplevels(subset(clusters, Subtype != "MSI-H"))
write.table("./Output/clusters_13_02.MSS.txt", x = clusters_MSS, row.names = F, quote = F, sep = "\t")


save.image("./R_Data/SNP_patient_find.RData")

#### END ####