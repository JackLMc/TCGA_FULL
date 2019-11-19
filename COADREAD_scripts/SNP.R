# A R script to investigate associations of SNPs with Patient groups
library(tidyverse)
library(UsefulFunctions)


SNP <- read.delim("./Data/SNPs/matrix_genotypes.final.txt")

SNP1 <- column_to_rownames(SNP, var = "snpID") %>% 
  as.matrix() %>%
  t()
rownames(SNP1) <- gsub("^X", "", rownames(SNP1))

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

# Need to find the sample types for each of the files, to avoid using primary tumours
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
  samp1 <- select(samp, -contains("portion"))
  list_of_samples[[i]] <- samp1
  }

meta_sample <- bind_rows(list_of_samples)
meta_data_samples <- meta_sample[, c("file_id", "submitter_id", "sample_type")]


META <- merge(meta_df1, meta_data_samples, by = "file_id")
META$Patient.ID <- samptopat(META$submitter_id)
META$Patient.ID <- gsub("-", ".", META$Patient.ID)
META$file_id <- gsub("-", ".", META$file_id)
META <- factorthese(META, colnames(META))

normals <- droplevels(subset(META, sample_type == "Blood Derived Normal"))
SNP_normals <- SNP1[rownames(SNP1) %in% normals$file_id, ]

hiCIRC <- read.csv("./Output/Patient_Subtypes.csv")

length(normals$Patient.ID)
nlevels(normals$Patient.ID)
## 623 (blood and matched normal)
## 512 are blood
## 111 are matched normal

hiCIRC_normals <- merge(normals, hiCIRC, by = "Patient.ID")
CIRC <- droplevels(subset(hiCIRC_normals, Subtype == "MSS-hiCIRC"))
MSS <- droplevels(subset(hiCIRC_normals, Subtype == "MSS"))
MSI <- droplevels(subset(hiCIRC_normals, Subtype == "MSI-H"))

nlevels(CIRC$Patient.ID)
nlevels(MSS$Patient.ID)
nlevels(MSI$Patient.ID)


# file_id_for_dupes <- droplevels(subset(META, Patient.ID == "TCGA.A6.6781"))$file_id %>% as.character()
# file_id_for_dupes <- gsub("-", ".", file_id_for_dupes)
# length(file_id_for_dupes)
# 
# SNP_dupes <- SNP_normals[rownames(SNP_normals) %in% file_id_for_dupes, ]
# dim(SNP_dupes)
# # View(head(SNP_dupes)[, 1:6])
# SNP_dupes[SNP_dupes[1, ] != SNP_dupes[2, ]]

# Some patients have duplications, take first instance of these patients
undup_hiCIRC <- hiCIRC_normals[!duplicated(hiCIRC_normals$Patient.ID),]$file_id %>% as.character
head(SNP_normals)[, 1:6]


working_SNPs <- SNP_normals[rownames(SNP_normals) %in% undup_hiCIRC, ]
df1 <- rownames_to_column(as.data.frame(working_SNPs), var = "file_id")
df2 <- merge(df1, hiCIRC_normals[, c("file_id", "Patient.ID")], by = "file_id") %>% column_to_rownames(., var = "Patient.ID")
mat1 <- df2[, !('%in%'(colnames(df2), "file_id")) ] %>% as.matrix %>% t() %>% as.data.frame()

str(mat1)
count(mat1, "TCGA.CM.6679")

df3 <- factorthese(mat1, colnames(mat1))


str(df3)
# df3 <- factorthese(df3, colnames(df3))
str(df3)

#


## Other information...
list_of_diagnoses <- list()
for(i in names(list_of_cases)){
  work <- list_of_cases[[i]]
  print(i)
  l <- work[[1]]$diagnoses[[1]]
  diag <- data.frame(as.list(unlist(l))) %>%
    mutate_all(as.character)
  rownames(diag) <- NULL
  diag$file_id <- i
  diag1 <- select(diag, -contains("portion"))
  list_of_diagnoses[[i]] <- diag1
}









