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
working_SNPs <- SNP_normals[rownames(SNP_normals) %in% undup_hiCIRC, ]


df1 <- rownames_to_column(as.data.frame(working_SNPs), var = "file_id")
df2 <- merge(df1, hiCIRC_normals[, c("file_id", "Patient.ID")], by = "file_id") %>% column_to_rownames(., var = "Patient.ID")
mat1 <- df2[, !('%in%'(colnames(df2), "file_id")) ] %>% as.matrix %>% as.data.frame()
df3 <- factorthese(mat1, colnames(mat1))

these_patients <- hiCIRC[hiCIRC$Patient.ID %in% rownames(mat1), c("Patient.ID", "Subtype")]
write.csv("./Output/Patients_for_SNPs.csv", row.names = F, x = these_patients)


# Table for 
df4 <- rownames_to_column(df3, var = "Patient.ID")
df5 <- merge(df4, these_patients, by = "Patient.ID")


rm(list=setdiff(ls(), c("df3", "df4", "df5", "mat1", "these_patients")))

save.image("./Output/SNP.RData")
load("./Output/SNP.RData")
library(tidyverse)

list_of_freq <- list()
c <- 1
for(i in levels(df5$Subtype)){
  print(i)
  work <- df5[levels(df5$Subtype) %in% i, ]
  work1 <- work[, !('%in%'(colnames(work), "Subtype"))]
  rownames(work1) <- NULL
  work2 <- column_to_rownames(work1, var = "Patient.ID")
  
  for(j in colnames(work2)){
    print(j)
    work3 <- work2[[j]]
    partA <- table(work3) %>% as.data.frame %>% mutate_all(as.character)
    partA$SNP <- j
    partA$Subtype <- i
    list_of_freq[[i]] <- partA
    colnames(list_of_freq[[i]]) <- c("Factor_levels", "Frequency", "SNP", "Subtype")
  }
  c <- c + 1
}

head(list_of_freq)

# Creating a frequency table - change to df3 if you want to leave the SNPs which are -1 in...
list_of_freq <- list()
for(i in colnames(df3)){
  print(i)
  work <- df3[[i]]
  partA <- table(work) %>% as.data.frame %>% mutate_all(as.character)
  partA$SNP <- i
  list_of_freq[[i]] <- partA
  colnames(list_of_freq[[i]]) <- c("Factor_levels", "Frequency", "SNP")
}




head(list_of_freq)
DF <- bind_rows(list_of_freq)
DF$Frequency <- as.numeric(DF$Frequency)

# save.image("./Data/SNPs/")

# -1 Factor_level is likely those which haven't been called successfully. 
patients_with_minus <- droplevels(subset(DF, Factor_levels == -1))$Patient.ID
this <- df3[, rownames(df3) %in% patients_with_minus]
View(this[, 1:10])

minus_1 <- this[-1, ] %>% head()[, 1:6]
minus <- numericthese(minus_1, colnames(minus_1))
minus[colSums(minus) != -257]
rownames(minus_1)


df3a <- df3[!('%in%'(rownames(df3), rownames(minus_1))), ]
dim(df3a)

# Creating a frequency table
list_of_freq <- list()
for(i in colnames(df3a)){
  print(i)
  work <- df3a[[i]]
  partA <- table(work) %>% as.data.frame %>% mutate_all(as.character)
  partA$SNP <- i
  list_of_freq[[i]] <- partA
  colnames(list_of_freq[[i]]) <- c("Factor_levels", "Frequency", "SNP")
}
DF <- bind_rows(list_of_freq)

DF1 <- droplevels(subset(DF, Factor_levels != -1))
DF2 <- merge(DF1, hiCIRC[, c("Patient.ID", "Subtype")], by = "Patient.ID")
DF2$Unique <- paste(DF2$Subtype, DF2$Factor_levels, sep = "_") %>% as.factor()

DF3 <- data.frame()
c <- 1
for(i in levels(DF2$Unique)){
  print(i)
  work <- droplevels(subset(DF2, Unique == i))
  summed <- sum(work$Frequency)
  DF3[c, "Unique"] <- i
  DF3[c, "Total_Frequency"] <- summed
  c <- c + 1
  }

DF4 <- separate(DF3, col = "Unique", into = c("Subtype", "Allele"), sep = "_") %>% 
  spread(., key = "Allele", value = "Total_Frequency", sep = "_")

colnames(DF4) <- gsub("Allele_", "", colnames(DF4))

DF4$Total <- rowSums(DF4[, !('%in%'(colnames(DF4), "Subtype"))]) 
col_total <- colSums(DF4[, !('%in%'(colnames(DF4), c("Subtype")))])

totals <- as.data.frame(col_total) %>% t
Subtype <- "Total"
add_on <- cbind(Subtype, totals) %>% as.data.frame()
add_on$`0` <- as.numeric(as.character(add_on$`0`))
add_on$`1` <- as.numeric(as.character(add_on$`1`))
add_on$`2` <- as.numeric(as.character(add_on$`2`))
add_on$`Total` <- as.numeric(as.character(add_on$`Total`))
rownames(add_on) <- NULL


DF5 <- rbind(DF4, add_on)
DF6 <- droplevels(subset(DF5, Subtype != "Total"))
rownames(DF6) <- NULL

DF7 <- column_to_rownames(DF6, var = "Subtype") 
DF8 <- DF7[rownames(DF7) != "MSI-H", ] # remove the MSI-H patients to see if that affects it.

chisq.test(DF7[, !('%in%'(colnames(DF8), "Total"))])


DF7[, !('%in%'(colnames(DF7), "Total"))]


# P-value is highly significant, if I remove the genotypes which are -1 from the few patients (2.2e-16)
# Remove the SNPs from the analysis all together, p-value is the same if you remove all the SNPs (2.2e-16)


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









