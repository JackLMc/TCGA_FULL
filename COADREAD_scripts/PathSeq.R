# PathSeq analysis ----
# Preprocessing ----
# Looking for the files
## Mount the drive
library(tidyverse)
library(UsefulFunctions)
# thousand.folders <- list.dirs(path = "/Volumes/2018/beggsa-tcgacolorectal/download_rest/bacterial_project/results", full.names = T)
# filelist1 <- sapply(thousand.folders, function(x){
#   list.files(x, pattern = ".pathseq.txt$", full.names = T)})
# filelist = unlist(filelist1)
pathseq <- read.delim("./Data/PathSeq/Raw_Data_20190711/pathseq.txt", stringsAsFactors = F)

# # Read in files and combine
# lists <- lapply(filelist, read.delim, header = T, stringsAsFactors = F)
# listsDF <- lists
# lists <- listsDF
# names(lists) <- filelist
# names(lists) <- gsub("/Volumes/2018/beggsa-tcgacolorectal/download_rest/bacterial_project/results/", "", names(lists))
# names(lists) <- gsub("pathseq.txt", "bam", names(lists))

# Colnames = c(tax_id, taxonomy, type, name,
# kingdom, score, score_normalized, reads, unambiguous, reference_length)
GDC_convert <- read.delim("./Data/GDC_large_mapping_TCGA.txt")[, c("file_name", "cases.0.submitter_id", "cases.0.samples.0.sample_type",
                                                                   "data_category",
                                                                   "cases.0.samples.0.submitter_id",
                                                                   "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id")]

colnames(GDC_convert) <- c("File.Name", "Patient.ID", "Sample.Type", "Data.Type", "TCGA_Submitter_ID", "TCGA.ID")
GDC_convert$Patient.ID <- gsub("-", ".", GDC_convert$Patient.ID)
GDC_convert <- GDC_convert[GDC_convert$File.Name %in% names(lists), ]

# GDC_convert$Sample.ID <- gsub("_hg19.*$|_Illumina.*$|_gdc.*$", "", GDC_convert$File.Name)


# colnames(GDC_convert)[colnames(GDC_convert) == "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"] <- "file.ID"
GDC_convert_t <- droplevels(subset(GDC_convert, Sample.Type == "Metastatic" |
                                Sample.Type == "Primary Tumor" | Sample.Type == "Recurrent Tumor"))

GDC_convert_n <- droplevels(subset(GDC_convert, Sample.Type == "Blood Derived Normal" |
                                     Sample.Type == "Solid Tissue Normal"))

GDC_convert_t$Patient.ID <- as.factor(GDC_convert_t$Patient.ID)
GDC_convert_t$File.Name<- as.factor(GDC_convert_t$File.Name)
GDC_convert_t <- droplevels(GDC_convert_t)

# lists1 <- lists[names(lists) %in% GDC_convert$File.Name]
# lists1a <- lists1[sapply(lists1, function(x) dim(x)[1]) > 0]
# lists2 <- lapply(names(lists1a), 
#                   function(n, x){
#                     x[[n]]$File.Name <- n
#                     return (x[[n]])},
#                  lists1)
# 
# multi_join <- function(list_of_loaded_data, join_func, ...){
#   require("dplyr")
#   output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
#   return(output)}
# 
# pathseq <- multi_join(lists2, full_join, by = c("tax_id", "taxonomy", "type", "name",
#                                                     "kingdom", "score", "score_normalized",
#                                                     "reads", "unambiguous", "reference_length", "File.Name"))
pathseq1 <- merge(pathseq, GDC_convert, by = "File.Name") %>% droplevels()
pathseq1$Patient.ID <- as.factor(pathseq1$Patient.ID)
pathseq1$File.Name <- as.factor(pathseq1$File.Name)
nlevels(pathseq1$File.Name)
# pathseq1a <- droplevels(subset(pathseq1, Sample.Type == "Primary Tumor"))

# Num_Files <- data.frame()
# c <- 1
# for(i in levels(pathseq1a$Patient.ID)){
#   print(i)
#   work <- droplevels(subset(pathseq1a, Patient.ID == i))
#   num_files <- nlevels(work$File.ID)
#   Num_Files[c, "Patient.ID"] <- i
#   Num_Files[c, "Number"] <- num_files
#   c <- c + 1
# }
# 
# more_than_one <- droplevels(subset(Num_Files, Number > 1))$Patient.ID
# info_more <- GDC_convert_t[GDC_convert_t$Patient.ID %in% more_than_one, ]
# info_more$Sample.ID <- gsub("_hg19.*$|_Illumina.*$|_gdc.*$", "", info_more$File.Name)
# info_more1 <- subset(info_more, !grepl("_gapfillers", info_more$Sample.ID), drop = T) %>% droplevels() 
# info_more1$File.Name <- as.factor(info_more1$File.Name)



# Use Submitter IDs 
pathseq2 <- factorthese(pathseq1, c("Patient.ID", "kingdom", "name", "taxonomy", "type", "File.Name")) %>% droplevels()
pathaa <- pathseq2[grep("gapfillers", pathseq2$File.Name), ] %>% droplevels()

df <- data.frame()
c <- 1
for(i in levels(pathaa$TCGA_Submitter_ID)){
  print(i)
  work <- droplevels(subset(pathseq2, TCGA_Submitter_ID == i))
  levs <- nlevels(work$File.Name)
  df[c, "TCGA_Submitter_ID"] <- i
  df[c, "Number"] <- levs
  c <- c + 1
  }

pats_with_gapfillers <- droplevels(subset(df, Number > 1))$TCGA_Submitter_ID
GDC_convert[GDC_convert$TCGA_Submitter_ID %in% pats_with_gapfillers, ] %>% droplevels()

missing_pats <- droplevels(subset(df, Number == 1))$TCGA_Submitter_ID

head(missing_patient_files)

## NEED TO FIGURE OUT HOW TO DEAL WITH THE GAPFILLERS, SHOULD MERGE IN, BUT ALSO... WE DON'T HAVE ALL THE DATA

# missing_patient_files <- GDC_convert[GDC_convert$TCGA_Submitter_ID %in% missing_pats, ] %>% droplevels()
# write.csv(missing_patient_files, file = "./Output/Patients_missing_pathseq.csv", row.names = F)



pathseq2$Gap_filler <- gsub("_hg19.*$|_Illumina.*$|_gdc.*$", "", pathseq2$File.Name) %>% as.factor()
# pathseq2a <- subset(pathseq2, !grepl("_gapfillers", pathseq2$Sample.ID), drop = T) %>% droplevels() 
pathseq2$Sample.ID <- gsub("_hg19.*$|_Illumina.*$|_gdc.*$|_gapfillers.*$", "", pathseq2$File.Name) %>% as.factor()

pathseq2$Sample.ID <- as.factor(pathseq2$Sample.ID)
pathseq2$Gap_filler <- as.factor(pathseq2$Gap_filler)



df <- data.frame()
c <- 1
for(i in levels(pathseq2$Sample.ID)){
  print(i)
  work <- droplevels(subset(pathseq2, Sample.ID == i))
  levs <- nlevels(work$Gap_filler)
  df[c, "Sample.ID"] <- i
  df[c, "Number"] <- levs
  c <- c + 1
}

head(df)

df1 <- droplevels(subset(df, Number == 1))

pathseq2[grep("gapfillers", pathseq2$Gap_filler), ]

droplevels(subset(pathseq2, Sample.ID == "TCGA-A6-2674-01A-02D-1953-10"))

levels(pathseq2$Sample.ID)

pathseq2a$File.Name <- as.factor(pathseq2a$File.Name)

# nlevels(pathseq2$File.Name) 
# nlevels(pathseq2a$File.Name) # This is a lot of patients to lose. Need to rerun the analysis after merging the patient files
# Droplevels for the ones I have CIRC scored
pat_sub <- read.csv("./Output/Patient_Subtypes.csv")
pathseq2b <- droplevels(pathseq2a[pathseq2a$Patient.ID %in% pat_sub$Patient.ID, ])
pathseq3 <- pathseq2b[, c("Patient.ID", "name", "score",
                         "score_normalized", "reads", "kingdom", "type", "unambiguous", "Sample.ID")]
pathseq3$Sample.ID <- as.factor(pathseq3$Sample.ID)

head(pathseq3)


pathseq2[grep("gapfill", pathseq2$File.Name),]

# Fill in the missing data for microbes not found in those samples
microbe_taxa <- pathseq1[, c("kingdom", "type", "name")]
microbe_taxa <- microbe_taxa[!duplicated(microbe_taxa), ]

missing_microbes <- list()
for(i in levels(pathseq3$Sample.ID)){
  print(i)
  work <- droplevels(subset(pathseq3, Sample.ID == i))
  pat_microbes <- as.character(levels(work$name))
  all_microbes <- as.character(levels(pathseq3$name))
  see <- all_microbes[!(all_microbes %in% pat_microbes)]
  missing <- ifelse((see != 0), c(see), "None")
  missing_microbes[[i]] <- see
  }

missing_microbes1 <- lapply(missing_microbes, as.data.frame, stringsAsFactors = F)
missing_microbes2 <- lapply(missing_microbes1, setNames, c("name"))

missing_microbes3 <- lapply(names(missing_microbes2), 
                 function(n, x){
                   x[[n]]$Sample.ID <- n
                   return (x[[n]])},
                 missing_microbes2)

missing_microbes3 <- mapply(cbind, missing_microbes3, "score" = 0, SIMPLIFY = F)
missing_microbes3 <- mapply(cbind, missing_microbes3, "score_normalized" = 0, SIMPLIFY = F)
missing_microbes3 <- mapply(cbind, missing_microbes3, "reads" = 0, SIMPLIFY = F)
missing_microbes3 <- mapply(cbind, missing_microbes3, "unambiguous" = 0, SIMPLIFY = F)


missing_data_a <- multi_join(missing_microbes3[1:500], full_join, by = c("Sample.ID", "name", "score", "score_normalized",
                                                "reads", "unambiguous"))

missing_data_b <- multi_join(missing_microbes3[501:1082], full_join, by = c("Sample.ID", "name", "score", "score_normalized",
                                                                         "reads", "unambiguous"))

missing_data <- rbind(missing_data_a, missing_data_b)
missing_data1 <- merge(missing_data, microbe_taxa, by = c("name"))
samp_patient <- pathseq3[, c("Patient.ID", "Sample.ID")]
samp_pat <- samp_patient[!duplicated(samp_patient), ]
missing_data2 <- merge(missing_data1, samp_pat, by = "Sample.ID")
pathseq4 <- rbind(pathseq3, missing_data2)

# save.image("./R_Data/PathSeq.RData")

# LOAD IT UP ----
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

load("./R_Data/PathSeq.RData")

# Overall abundances ----
pathseq_tum <- pathseq4[!('%in%'(pathseq4$Sample.ID, GDC_convert_n$Sample.ID)),]
superk <- droplevels(subset(pathseq_tum, type == "superkingdom"))

SK <- merge(superk, pat_sub, by = "Patient.ID")
for(i in levels(SK$kingdom)){
  print(i)
  work <- droplevels(subset(SK, kingdom == i))
  temp_plot <- ggplot(work, aes(x = Subtype, y = score_normalized)) +
    geom_boxplot(alpha = 0.5, width = 0.2) + 
    geom_violin(aes(Subtype, fill = Subtype),
                scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) +
    labs(x = "Subtype", y = paste(i, "Normalised score")) +
    theme_bw() +
    theme(axis.text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.direction = "horizontal", legend.position = "top") + 
    stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif", method = "wilcox.test")
  filen <- paste0(i, ".pdf")
  ggplot2:: ggsave(filen, plot = temp_plot, device = "pdf",
                   path = "./Figures/PathSeq/Superkingdom",
                   height = 6, width = 6)
  }

# Diversity indexes -----
# Diversity for tumour versus normal
pathseq5 <- merge(pathseq4, GDC_convert[, c("Sample.ID", "Sample.Type")], by = "Sample.ID")

pathseq5a <- droplevels(subset(pathseq5, kingdom != "root" & Sample.Type != "Recurrent Tumor"))

my_comparisons1 <- list(c("Blood Derived Normal","Primary Tumor"),
                        c("Blood Derived Normal", "Solid Tissue Normal"),
                        c("Primary Tumor", "Solid Tissue Normal"))

library(vegan)
for(i in levels(pathseq5a$kingdom)){
  print(i)
  work <- droplevels(subset(pathseq5a, kingdom == i & type == "species"))
  thesecols <- c("Sample.ID", "name", "reads")
  div1 <- work[, colnames(work) %in% thesecols]
  div2 <- spread(div1, key = "name", value = "reads") %>% column_to_rownames(., var = "Sample.ID")
  div3 <- diversity(div2, index = "shannon") %>% as.data.frame() %>% rownames_to_column(., var = "Sample.ID")
  colnames(div3)[colnames(div3) == "."] <- "shannon"
  div4 <- merge(div3, GDC_convert, by = "Sample.ID")
  
  temp_plot <- ggplot(div4, aes(x = Sample.Type, y = shannon)) +
    geom_boxplot(alpha = 0.5, width = 0.2) + 
    geom_violin(aes(Sample.Type, fill = Sample.Type),
                scale = "width", alpha = 0.8) +
    # scale_fill_manual(values = cbcols) +
    labs(x = "Sample.Type", y = paste(i, "Shannon diversity")) +
    theme_bw() +
    theme(axis.text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.direction = "horizontal", legend.position = "top") + 
    stat_compare_means(comparisons = my_comparisons1,
                       label = "p.signif", method = "wilcox.test")
  filen <- paste0(i, ".pdf")
  ggplot2:: ggsave(filen, plot = temp_plot, device = "pdf",
                   path = "./Figures/PathSeq/Diversity/Sample_Type",
                   height = 6, width = 6)
}

# Diverstiy across MMR
samp_pat_sub <- merge(samp_pat, pat_sub, by = "Patient.ID")
pathseq4a <- merge(pathseq4, GDC_convert, by = c("Sample.ID", "Patient.ID"))
pathseq4a1 <- droplevels(subset(pathseq4a, kingdom != "root" & Sample.Type != "Recurrent Tumor"))
pathseq4b <- merge(pathseq4a1, samp_pat_sub, by =  c("Sample.ID", "Patient.ID"))
pathseq6 <- pathseq4b
pathseq6$Unique <- paste(pathseq6$Subtype, pathseq6$Sample.Type, sep = "_") %>% as.factor()

library(vegan)
for(i in levels(pathseq6$kingdom)){
  print(i)
  work <- droplevels(subset(pathseq6, kingdom == i & type == "species"))
  thesecols <- c("Sample.ID", "name", "reads")
  div1 <- work[, colnames(work) %in% thesecols]
  div2 <- spread(div1, key = "name", value = "reads") %>% column_to_rownames(., var = "Sample.ID")
  div3 <- diversity(div2, index = "shannon") %>% as.data.frame() %>% rownames_to_column(., var = "Sample.ID")
  colnames(div3)[colnames(div3) == "."] <- "shannon"
  
  MMR <- work[, colnames(work) %in% c("Sample.ID", "Unique")]
  div4 <- merge(div3, MMR, by = "Sample.ID")
  
  temp_plot <- ggplot(div4, aes(x = Unique, y = shannon)) +
    geom_boxplot(alpha = 0.5, width = 0.2) + 
    geom_violin(aes(Unique, fill = Unique),
                scale = "width", alpha = 0.8) +
   # scale_fill_manual(values = cbcols) +
    labs(x = "Subtype", y = paste(i, "Shannon diversity")) +
    theme_bw() +
    theme(axis.text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.direction = "horizontal", legend.position = "top") + 
    stat_compare_means(comparisons = list(c("MSI-H_Blood Derived Normal", "MSI-H_Solid Tissue Normal"),
                                          c("MSI-H_Blood Derived Normal", "MSI-H_Primary Tumor"),
                                          c("MSI-H_Solid Tissue Normal", "MSI-H_Primary Tumor"),
                                          c("MSS_Primary Tumor", "MSS_Blood Derived Normal"),
                                          c("MSS_Primary Tumor", "MSS_Solid Tissue Normal"),
                                          c("MSS_Solid Tissue Normal", "MSS_Blood Derived Normal"),
                                          c("MSS-hiCIRC_Blood Derived Normal", "MSS-hiCIRC_Primary Tumor"),
                                          c("MSS-hiCIRC_Blood Derived Normal", "MSS-hiCIRC_Solid Tissue Normal"),
                                          c("MSS-hiCIRC_Primary Tumor", "MSS-hiCIRC_Solid Tissue Normal")),
                       label = "p.signif", method = "wilcox.test")
  filen <- paste0(i, ".pdf")
  ggplot2:: ggsave(filen, plot = temp_plot, device = "pdf",
                   path = "./Figures/PathSeq/Diversity/MMR_Samp_Type",
                   height = 6, width = 6)
  }

div <- droplevels(subset(pathseq4b, type == "species"))
library(vegan)
thesecols <- c("Sample.ID", "name", "reads")
div1 <- div[, colnames(div) %in% thesecols]
div2 <- spread(div1, key = "name", value = "reads") %>%
  column_to_rownames(., var = "Sample.ID") %>%
  diversity(., index = "shannon") %>% as.data.frame() %>% rownames_to_column(., var = "Sample.ID")

colnames(div2)[colnames(div2) == "."] <- "shannon"
div4 <- merge(div3, samp_pat_sub, by = "Sample.ID")

pdf("./Figures/PathSeq/Diversity/Root.pdf", height = 6, width = 6)
ggplot(div4, aes(x = Subtype, y = shannon)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(aes(Subtype, fill = Subtype),
              scale = "width", alpha = 0.8) +
  scale_fill_manual(values = cbcols) +
  labs(x = "Subtype", y = "root Shannon index") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top") + 
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test")
dev.off()

# Compare by Sample.Type, and MMR status
samp_pat_sub <- merge(samp_pat, pat_sub, by = "Patient.ID")
pathseq4a <- droplevels(subset(pathseq_tum, kingdom != "root"))
pathseq4b <- pathseq4a[pathseq4a$Patient.ID %in% samp_pat_sub$Patient.ID, ]

# PCA of data
set.seed(1)

prin <- droplevels(subset(pathseq6, type == "species"))
thesecols <- c("Sample.ID", "name", "reads")
prin1 <- prin[, colnames(prin) %in% thesecols]
prin1a <- prin1[prin1$Sample.ID %in% samp_pat_sub$Sample.ID, ]
nlevels(prin1a$Sample.ID)

prin2 <- spread(prin1, key = "name", value = "reads") %>% column_to_rownames(., var = "Sample.ID")
prin_comp <- prcomp(prin2, scale. = T)

library(ggbiplot)
pdf("./Figures/PathSeq/.pdf", height = 6, width = 6)
ggbiplot(prin_comp, obs.scale = 1, var.scale = 1, 
         #groups = Pcluster,
         circle = T, var.axes = F) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top") 
dev.off()

## Remove the outlier
rotations <- prin_comp$x %>% as.data.frame() %>% rownames_to_column(., var = "Sample.ID")
outlier <- droplevels(subset(rotations, PC2 < -200 | PC1 > 200))$Sample.ID
prin3 <- prin2[!('%in%'(rownames(prin2), outlier)), ]
prin4 <- prin3[, colSums(prin3 != 0) > 0]
prin_comp <- prcomp(prin4, scale. = T)

samp_pat_sub1 <- samp_pat_sub[!('%in%'(samp_pat_sub$Sample.ID, outlier)), ]


rotations <- prin_comp$x %>% as.data.frame() %>% rownames_to_column(., var = "Sample.ID")
rotations$Sample.ID <- as.factor(rotations$Sample.ID)

samp_pat_sub2 <- samp_pat_sub1['%in%'(samp_pat_sub1$Sample.ID, rotations$Sample.ID), ]
samp_pat_sub3 <- merge(samp_pat_sub2, GDC_convert[, c("Sample.ID", "Sample.Type")], by = "Sample.ID")
samp_pat_sub3$Unique <- paste(samp_pat_sub3$Subtype, samp_pat_sub3$Sample.Type, sep = "_")


ggbiplot(prin_comp, obs.scale = 1, var.scale = 1, 
        groups = samp_pat_sub3$Unique,
         circle = T, var.axes = F) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top") 
dev.off()


## PCA doesn't find it...

# tSNE
## Perform Phenograph and kmeans
library(Rphenograph)
set.seed(1)
a <- cbind(prin4, cluster = factor(Rphenograph(prin4)[[2]]$membership))

clusterP <- a[, "cluster"]

# RTsne
set.seed(1)
library(Rtsne)
tsne_out <- Rtsne(as.matrix(prin4))
tsne_dimensions <- as.data.frame(tsne_out$Y)
colnames(tsne_dimensions) <- c("Dim1", "Dim2")


## tSNE plot - looks the same as PCA just on a different axis
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# pdf("./Figures/Clustering/CIRC_PhenoGraph_tsNE.pdf")
ggplot(tsne_dimensions, aes(x = Dim1, y = Dim2, colour = samp_pat_sub3$Unique)) +
  geom_point(size = 4, alpha = 0.8, pch = 20) +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top")
dev.off()

# tSNE also doesn't find distinct clusters of patients based on the MSS status

# Take into account taxonomy ----
# data is rownames(species), Genus, Family, Order, Superorder. Subclass

data(dune)
data(dune.taxon)
taxdis <- taxa2dist(dune.taxon, varstep=TRUE)


mod <- taxondive(dune, taxdis)


# ## Double-check
# output <- data.frame(stringsAsFactors = F)
# c <- 1
# for(i in levels(pathseq3$Patient.ID)){
#   work <- droplevels(subset(pathseq3, Patient.ID == i))
#   these <- nlevels(work$name)
#   output[c, "Patient.ID"] <- i
#   output[c, "num_of_tax"] <- these
#   c <- c + 1
# }

# Merge with the pat_sub file
library(reshape2)
droplevels(subset(pathseq6, Subtype == "MSS"))$Patient.ID %>% nlevels()
droplevels(subset(pathseq6, Subtype == "MSS-hiCIRC"))$Patient.ID %>% nlevels()
droplevels(subset(pathseq6, Subtype == "MSI-H"))$Patient.ID %>% nlevels()

empty_ <- data.frame("name" = character(),
                     "num_lev" = double(),
                     stringsAsFactors = F)
c <- 1
for(i in levels(pathseq6$name)){
  print(i)
  working <- droplevels(subset(pathseq6, name == i))
  empty_[c, "name"] <- i
  empty_[c, "num_lev"] <- nlevels(working$Patient.ID)
  c <- c + 1
}

# Combine stats into a list, not working because some bugs are 0
conversion <- samp_pat_sub3
pathseq6 <- merge(pathseq4, conversion, by = "Sample.ID")
pathseq6a <- droplevels(subset(pathseq6, Sample.Type == "Primary Tumor"))

maximum_nums <- data.frame(stringsAsFactors = F)
c <- 1
for(i in levels(pathseq6a$name)){
  print(i)
  work <- droplevels(subset(pathseq6a, name == i))
  max_num <- max(work$score)
  maximum_nums[c, "name"] <- i
  maximum_nums[c, "max"] <- max_num
  c <- c + 1
}

remove_these <- droplevels(subset(maximum_nums, max == 0))$name
pathseq5 <- pathseq6a[!('%in%'(pathseq6a$name, remove_these)), ] %>% droplevels()

# save.image("./R_Data/PathSeq.RData")
load("./R_Data/PathSeq.RData")
head(pathseq)
# check bacteroides/firmicutes ratio, prediction increased bacteroides in hiCIRC patients.

library(ggpubr)
stat_list <- list()
c <- 1
for(i in levels(pathseq5$name)){
  name <- basename(i)
  cat("Processing", i, "\n")
  workingon <- droplevels(subset(pathseq5, name == i))
  # assign your ggplot call to the i"th position in the list
  x <- compare_means(score ~ Subtype, data = workingon, method = "wilcox.test")
  y <- as.data.frame(x)
  stat_list[[i]]  <- cbind(i, y)
  c <- c + 1}


# Bind and remove row names
z <- do.call(rbind, stat_list)
rownames(z) <- c()

sigs <- droplevels(subset(z, p.signif != "ns" & p.signif != "NS" & p.signif != "*"))

# Write out the statistics
write.csv("Output/significance_pathseq.csv", x = sigs, row.names = F)

these <- data.frame(stringsAsFactors = F)
c <- 1
for(j in levels(sigs$i)){
  work <- droplevels(subset(sigs, i = j))
  these[c, "patho"] <- j
  these[c, "num_rows"] <- nrow(work)
  c <- c + 1
}


head(sigs)
MSS_ <- droplevels(subset(sigs, group2 == "MSS-hiCIRC" & group1 == "MSS" | group1 == "MSS-hiCIRC" & group2 == "MSS"))

hiCIRC_H <- droplevels(subset(sigs, group2 == "MSS-hiCIRC" & group1 == "MSI-H" | group1 == "MSS-hiCIRC" & group2 == "MSI-H"))

bugs <- levels(MSS_$i)[levels(MSS_$i) %in% levels(hiCIRC_H$i)]
# try2 <- levels(hiCIRC_H$i)[levels(hiCIRC_H$i) %in% levels(MSS_$i)]

hiCIRC_sp <- z[z$i %in% bugs, ]


droplevels(subset(z, i == "Janthinobacterium_sp._B9-8"))
droplevels(subset(z, i == "Lactobacillus_hominis"))
droplevels(subset(z, i == "Lactobacillus_hominis_DSM_23910_=_CRBIP_24.179"))
droplevels(subset(z, i == "Methylotenera_sp._N17"))
droplevels(subset(z, i == "Nocardioidaceae"))
droplevels(subset(z, i == "Pseudogulbenkiania_sp._MAI-1"))
droplevels(subset(z, i == "Staphylococcus_sp._HMSC075F12"))
droplevels(subset(z, i == "Staphylococcus_sp._HMSC076B11"))
droplevels(subset(z, i == "Staphylococcus_sp._HMSC076E07"))
droplevels(subset(z, i == "Xanthomonas_citri_pv._glycines"))
droplevels(subset(z, i == "Xanthomonas_citri_pv._glycines_str._12-2"))
droplevels(subset(z, i == "Zymobacter_group"))


View(MSS_)

try <- pathseq5[pathseq5$name %in% bugs, ] %>% droplevels()
levels(try$Sample.Type)



c <- 1
for(i in levels(try$name)){
  work <- droplevels(subset(try, name == i))
  work$Rank <- rank(work$score)
  temp_plot <- ggplot(work, aes(x = Subtype, y = Rank)) +
    geom_boxplot(alpha = 0.5, width = 0.2) + 
    geom_violin(aes(Subtype, fill = Subtype),
                scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) +
    labs(x = "Subtype", y = "Normalised bact") +
    theme_bw() +
    theme(axis.text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.direction = "horizontal", legend.position = "top") + 
    stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif", method = "wilcox.test")
  filen <- paste0("./Figures/PathSeq/", i, ".pdf")
  ggsave(filen, temp_plot)
  c <- c + 1
}


# Next 'level'
head(pathseq3)
levels(pathseq3$type)

droplevels(subset(pathseq3, type == "superkingdom"))

empty_ <- data.frame("name" = character(),
                     "num_lev" = double(),
                     stringsAsFactors = F)
c <- 1
for(i in levels(pathseq4$name)){
  print(i)
  working <- droplevels(subset(pathseq4, name == i))
  empty_[c, "name"] <- i
  empty_[c, "num_lev"] <- nlevels(working$Patient.ID)
  c <- c + 1
}

empty_[empty_$num_lev != 529, ]


try <- pathseq4[grepl("Clostridium", pathseq4$name, ignore.case = T), ] %>% droplevels()

Ecol <- droplevels(subset(pathseq4, name == "Buty"))


# Coprobacillus_sp._D6, Bifidobacterium_animalis_subsp._lactis (not exact to one in paper), Clostridium_sp._ASF356, Clostridium_saccharobutylicum (not exact)
# Eubacterium_sp._3_1_31 (up MSI-h), Erysipelotrichaceae_bacterium_21_3 (up in MSI-H, slightly), Firmicutes_bacterium_ASF500,
# Clostridium_hathewayi 12489931 (not in), Ruminococcus_gnavus_AGR2154 (MSI-H), Ruminococcus bromii (not in) or Ruminococcus obeum
# Subdoligranulum_sp._4_3_54A2FAA (highly expressed in MSI-H), Bifidobacterium_breve (not eact), Clostridium symbiosium (not in)
# Bacteroides_dorei (not exact, but up in MSI-H)



# Combine stats into a list, not working because some bugs are 0
maximum_nums <- data.frame(stringsAsFactors = F)
c <- 1
for(i in levels(pathseq4$name)){
  print(i)
  work <- droplevels(subset(pathseq4, name == i))
  max_num <- max(work$score)
  maximum_nums[c, "name"] <- i
  maximum_nums[c, "max"] <- max_num
  c <- c + 1
}

remove_these <- droplevels(subset(maximum_nums, max == 0))$name
pathseq5 <- pathseq4[!('%in%'(pathseq4$name, remove_these)), ] %>% droplevels()

# save.image("./R_Data/PathSeq.RData")
load("./R_Data/PathSeq.RData")

library(ggpubr)
stat_list <- list()
c <- 1
for(i in levels(pathseq5$name)){
  name <- basename(i)
  cat("Processing", i, "\n")
  workingon <- droplevels(subset(pathseq5, name == i))
  # assign your ggplot call to the i"th position in the list
  x <- compare_means(score ~ Subtype, data = workingon, method = "wilcox.test")
  y <- as.data.frame(x)
  stat_list[[i]]  <- cbind(i, y)
  c <- c + 1}


# Bind and remove row names
z <- do.call(rbind, stat_list)
rownames(z) <- c()

sigs <- droplevels(subset(z, p.signif != "ns" & p.signif != "NS" & p.signif != "*"))

# Write out the statistics
write.csv("Output/significance_pathseq.csv", x = sigs, row.names = F)

these <- data.frame(stringsAsFactors = F)
c <- 1
for(j in levels(sigs$i)){
  work <- droplevels(subset(sigs, i = j))
  these[c, "patho"] <- j
  these[c, "num_rows"] <- nrow(work)
  c <- c + 1
}


# Normals
head(GDC_convert_n)
GDC_convert_n_t <- droplevels(subset(GDC_convert_n, Sample.Type == "Solid Tissue Normal"))

lists1 <- lists[names(lists) %in% GDC_convert_n_t$File.Name]

lists1a <- lists1[sapply(lists1, function(x) dim(x)[1]) > 0]
lists2 <- lapply(names(lists1a), 
                 function(n, x){
                   x[[n]]$File.Name <- n
                   return (x[[n]])},
                 lists1)


multi_join <- function(list_of_loaded_data, join_func, ...){
  require("dplyr")
  output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
  return(output)
}

pathseq <- multi_join(lists2, full_join, by = c("tax_id", "taxonomy", "type", "name",
                                                "kingdom", "score", "score_normalized",
                                                "reads", "unambiguous", "reference_length", "File.Name"))

pathseq1 <- merge(pathseq, GDC_convert, by = "File.Name") %>% droplevels()
pathseq1$Patient.ID <- as.factor(pathseq1$Patient.ID)
pathseq1$kingdom <- as.factor(pathseq1$kingdom)
pathseq1$name <- as.factor(pathseq1$name)
pathseq1$taxonomy <- as.factor(pathseq1$taxonomy)
pathseq1$Patient.ID <- as.factor(pathseq1$Patient.ID)
pathseq1$type <- as.factor(pathseq1$type)


pathseq2 <- pathseq1[, c("Patient.ID", "name", "score",
                         "score_normalized", "reads", "kingdom", "type", "unambiguous", "reference_length")]

# Fill in the missing data for microbes not found in those samples
microbe_taxa <- pathseq1[, c("kingdom", "type", "name", "reference_length")]
microbe_taxa <- microbe_taxa[!duplicated(microbe_taxa), ]

missing_microbes <- list()
for(i in levels(pathseq2$Patient.ID)){
  work <- droplevels(subset(pathseq2, Patient.ID == i))
  pat_microbes <- as.character(levels(work$name))
  all_microbes <- as.character(levels(pathseq2$name))
  see <- all_microbes[!(all_microbes %in% pat_microbes)]
  missing <- ifelse((see != 0), c(see), "None")
  missing_microbes[[i]] <- missing
}

missing_microbes1 <- lapply(missing_microbes, as.data.frame, stringsAsFactors = F)
missing_microbes2 <- lapply(missing_microbes1, setNames, c("name"))

missing_microbes3 <- lapply(names(missing_microbes2), 
                            function(n, x){
                              x[[n]]$Patient.ID <- n
                              return (x[[n]])},
                            missing_microbes2)

missing_microbes3 <- mapply(cbind, missing_microbes3, "score" = 0, SIMPLIFY = F)
missing_microbes3 <- mapply(cbind, missing_microbes3, "score_normalized" = 0, SIMPLIFY = F)
missing_microbes3 <- mapply(cbind, missing_microbes3, "reads" = 0, SIMPLIFY = F)
missing_microbes3 <- mapply(cbind, missing_microbes3, "unambiguous" = 0, SIMPLIFY = F)


missing_data <- multi_join(missing_microbes3, full_join, by = c("Patient.ID", "name", "score", "score_normalized",
                                                                "reads", "unambiguous"))
missing_data <- merge(missing_data, microbe_taxa, by = "name")
pathseq3 <- rbind(pathseq2, missing_data)

# ## Double-check
# output <- data.frame(stringsAsFactors = F)
# c <- 1
# for(i in levels(pathseq3$Patient.ID)){
#   work <- droplevels(subset(pathseq3, Patient.ID == i))
#   these <- nlevels(work$name)
#   output[c, "Patient.ID"] <- i
#   output[c, "num_of_tax"] <- these
#   c <- c + 1
# }

# Merge with the pat_sub file
pat_sub <- read.csv("./Output/Patient_Subtypes.csv")
pathseq4 <- merge(pathseq3, pat_sub, by = "Patient.ID")
pathseq4 <-  pathseq4[, c("Patient.ID", "kingdom", "type", "name",
                          "reference_length", "score", "score_normalized",
                          "reads", "unambiguous", "CIRC_Genes", "Subtype")] %>% droplevels()


library(reshape2)

droplevels(subset(pathseq4, Subtype == "MSS"))$Patient.ID %>% nlevels()
droplevels(subset(pathseq4, Subtype == "MSS-hiCIRC"))$Patient.ID %>% nlevels()
droplevels(subset(pathseq4, Subtype == "MSI-H"))$Patient.ID %>% nlevels()

empty_ <- data.frame("name" = character(),
                     "num_lev" = double(),
                     stringsAsFactors = F)
c <- 1
for(i in levels(pathseq4$name)){
  print(i)
  working <- droplevels(subset(pathseq4, name == i))
  empty_[c, "name"] <- i
  empty_[c, "num_lev"] <- nlevels(working$Patient.ID)
  c <- c + 1
}

empty_[empty_$num_lev != 529, ]

# Combine stats into a list, not working because some bugs are 0
maximum_nums <- data.frame(stringsAsFactors = F)
c <- 1
for(i in levels(pathseq4$name)){
  print(i)
  work <- droplevels(subset(pathseq4, name == i))
  max_num <- max(work$score)
  maximum_nums[c, "name"] <- i
  maximum_nums[c, "max"] <- max_num
  c <- c + 1
}

remove_these <- droplevels(subset(maximum_nums, max == 0))$name
pathseq5 <- pathseq4[!('%in%'(pathseq4$name, remove_these)), ] %>% droplevels()

# save.image("./R_Data/PathSeq.RData")
load("./R_Data/PathSeq.RData")

library(ggpubr)
stat_list <- list()
c <- 1
for(i in levels(pathseq5$name)){
  name <- basename(i)
  cat("Processing", i, "\n")
  workingon <- droplevels(subset(pathseq5, name == i))
  # assign your ggplot call to the i"th position in the list
  x <- compare_means(score ~ Subtype, data = workingon, method = "wilcox.test")
  y <- as.data.frame(x)
  stat_list[[i]]  <- cbind(i, y)
  c <- c + 1}


# Bind and remove row names
z1 <- do.call(rbind, stat_list)
rownames(z1) <- c()

sigs <- droplevels(subset(z1, p.signif != "ns" & p.signif != "NS" & p.signif != "*"))

# Write out the statistics
write.csv("Output/significance_pathseq.csv", x = sigs, row.names = F)

these <- data.frame(stringsAsFactors = F)
c <- 1
for(j in levels(sigs$i)){
  work <- droplevels(subset(sigs, i = j))
  these[c, "patho"] <- j
  these[c, "num_rows"] <- nrow(work)
  c <- c + 1
}


head(sigs)
MSS_ <- droplevels(subset(sigs, group2 == "MSS-hiCIRC" & group1 == "MSS" | group1 == "MSS-hiCIRC" & group2 == "MSS"))

hiCIRC_H <- droplevels(subset(sigs, group2 == "MSS-hiCIRC" & group1 == "MSI-H" | group1 == "MSS-hiCIRC" & group2 == "MSI-H"))

bugs1 <- levels(MSS_$i)
bugs2 <- levels(hiCIRC_H)

try <- pathseq5[pathseq5$name %in% bugs1, ] %>% droplevels()
try1 <- pathseq5[pathseq5$name %in% bugs2, ] %>% droplevels()

c <- 1
for(i in levels(try$name)){
  work <- droplevels(subset(try, name == i))
  work$Rank <- rank(work$score)
  temp_plot <- ggplot(work, aes(x = Subtype, y = Rank)) +
    geom_boxplot(alpha = 0.5, width = 0.2) + 
    geom_violin(aes(Subtype, fill = Subtype),
                scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) +
    labs(x = "Subtype", y = "Normalised bact") +
    theme_bw() +
    theme(axis.text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.direction = "horizontal", legend.position = "top") + 
    stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif", method = "wilcox.test")
  filen <- paste0("./Figures/PathSeq/Norm/", i, ".pdf")
  ggsave(filen, temp_plot)
  c <- c + 1
  }

# Blood
levels(GDC_convert_n$Sample.Type)
GDC_convert_n_b <- droplevels(subset(GDC_convert_n, Sample.Type == "Blood Derived Normal"))

lists1 <- lists[names(lists) %in% GDC_convert_n_b$File.Name]

lists1a <- lists1[sapply(lists1, function(x) dim(x)[1]) > 0]
lists2 <- lapply(names(lists1a), 
                 function(n, x){
                   x[[n]]$File.Name <- n
                   return (x[[n]])},
                 lists1)


multi_join <- function(list_of_loaded_data, join_func, ...){
  require("dplyr")
  output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
  return(output)
}

pathseq <- multi_join(lists2, full_join, by = c("tax_id", "taxonomy", "type", "name",
                                                "kingdom", "score", "score_normalized",
                                                "reads", "unambiguous", "reference_length", "File.Name"))

pathseq1 <- merge(pathseq, GDC_convert, by = "File.Name") %>% droplevels()
pathseq1$Patient.ID <- as.factor(pathseq1$Patient.ID)
pathseq1$kingdom <- as.factor(pathseq1$kingdom)
pathseq1$name <- as.factor(pathseq1$name)
pathseq1$taxonomy <- as.factor(pathseq1$taxonomy)
pathseq1$Patient.ID <- as.factor(pathseq1$Patient.ID)
pathseq1$type <- as.factor(pathseq1$type)


pathseq2 <- pathseq1[, c("Patient.ID", "name", "score",
                         "score_normalized", "reads", "kingdom", "type", "unambiguous", "reference_length")]

# Fill in the missing data for microbes not found in those samples
microbe_taxa <- pathseq1[, c("kingdom", "type", "name", "reference_length")]
microbe_taxa <- microbe_taxa[!duplicated(microbe_taxa), ]

missing_microbes <- list()
for(i in levels(pathseq2$Patient.ID)){
  work <- droplevels(subset(pathseq2, Patient.ID == i))
  pat_microbes <- as.character(levels(work$name))
  all_microbes <- as.character(levels(pathseq2$name))
  see <- all_microbes[!(all_microbes %in% pat_microbes)]
  missing <- ifelse((see != 0), c(see), "None")
  missing_microbes[[i]] <- missing
}

missing_microbes1 <- lapply(missing_microbes, as.data.frame, stringsAsFactors = F)
missing_microbes2 <- lapply(missing_microbes1, setNames, c("name"))

missing_microbes3 <- lapply(names(missing_microbes2), 
                            function(n, x){
                              x[[n]]$Patient.ID <- n
                              return (x[[n]])},
                            missing_microbes2)

missing_microbes3 <- mapply(cbind, missing_microbes3, "score" = 0, SIMPLIFY = F)
missing_microbes3 <- mapply(cbind, missing_microbes3, "score_normalized" = 0, SIMPLIFY = F)
missing_microbes3 <- mapply(cbind, missing_microbes3, "reads" = 0, SIMPLIFY = F)
missing_microbes3 <- mapply(cbind, missing_microbes3, "unambiguous" = 0, SIMPLIFY = F)


missing_data <- multi_join(missing_microbes3, full_join, by = c("Patient.ID", "name", "score", "score_normalized",
                                                                "reads", "unambiguous"))
missing_data <- merge(missing_data, microbe_taxa, by = "name")
pathseq3 <- rbind(pathseq2, missing_data)

# ## Double-check
# output <- data.frame(stringsAsFactors = F)
# c <- 1
# for(i in levels(pathseq3$Patient.ID)){
#   work <- droplevels(subset(pathseq3, Patient.ID == i))
#   these <- nlevels(work$name)
#   output[c, "Patient.ID"] <- i
#   output[c, "num_of_tax"] <- these
#   c <- c + 1
# }

# Merge with the pat_sub file
pat_sub <- read.csv("./Output/Patient_Subtypes.csv")
pathseq4 <- merge(pathseq3, pat_sub, by = "Patient.ID")
pathseq4 <-  pathseq4[, c("Patient.ID", "kingdom", "type", "name",
                          "reference_length", "score", "score_normalized",
                          "reads", "unambiguous", "CIRC_Genes", "Subtype")] %>% droplevels()


library(reshape2)

droplevels(subset(pathseq4, Subtype == "MSS"))$Patient.ID %>% nlevels()
droplevels(subset(pathseq4, Subtype == "MSS-hiCIRC"))$Patient.ID %>% nlevels()
droplevels(subset(pathseq4, Subtype == "MSI-H"))$Patient.ID %>% nlevels()

empty_ <- data.frame("name" = character(),
                     "num_lev" = double(),
                     stringsAsFactors = F)
c <- 1
for(i in levels(pathseq4$name)){
  print(i)
  working <- droplevels(subset(pathseq4, name == i))
  empty_[c, "name"] <- i
  empty_[c, "num_lev"] <- nlevels(working$Patient.ID)
  c <- c + 1
}

empty_[empty_$num_lev != 504, ]

# Combine stats into a list, not working because some bugs are 0
maximum_nums <- data.frame(stringsAsFactors = F)
c <- 1
for(i in levels(pathseq4$name)){
  print(i)
  work <- droplevels(subset(pathseq4, name == i))
  max_num <- max(work$score)
  maximum_nums[c, "name"] <- i
  maximum_nums[c, "max"] <- max_num
  c <- c + 1
}

remove_these <- droplevels(subset(maximum_nums, max == 0))$name
pathseq5 <- pathseq4[!('%in%'(pathseq4$name, remove_these)), ] %>% droplevels()

library(ggpubr)
stat_list <- list()
c <- 1
for(i in levels(pathseq5$name)){
  name <- basename(i)
  cat("Processing", i, "\n")
  workingon <- droplevels(subset(pathseq5, name == i))
  # assign your ggplot call to the i"th position in the list
  x <- compare_means(score ~ Subtype, data = workingon, method = "wilcox.test")
  y <- as.data.frame(x)
  stat_list[[i]]  <- cbind(i, y)
  c <- c + 1}


# Bind and remove row names
z2 <- do.call(rbind, stat_list)
rownames(z2) <- c()

sigs <- droplevels(subset(z2, p.signif != "ns" & p.signif != "NS" & p.signif != "*"))

# Write out the statistics
# write.csv("Output/significance_pathseq.csv", x = sigs, row.names = F)

these <- data.frame(stringsAsFactors = F)
c <- 1
for(j in levels(sigs$i)){
  work <- droplevels(subset(sigs, i = j))
  these[c, "patho"] <- j
  these[c, "num_rows"] <- nrow(work)
  c <- c + 1
}


head(sigs)
MSS_ <- droplevels(subset(sigs, group2 == "MSS-hiCIRC" & group1 == "MSS" | group1 == "MSS-hiCIRC" & group2 == "MSS"))

hiCIRC_H <- droplevels(subset(sigs, group2 == "MSS-hiCIRC" & group1 == "MSI-H" | group1 == "MSS-hiCIRC" & group2 == "MSI-H"))

bugs1 <- levels(MSS_$i)
bugs2 <- levels(hiCIRC_H)

try <- pathseq5[pathseq5$name %in% bugs1, ] %>% droplevels()
try1 <- pathseq5[pathseq5$name %in% bugs2, ] %>% droplevels()
levels(try$name)
c <- 1
for(i in levels(try$name)){
  work <- droplevels(subset(try, name == i))
  work$Rank <- rank(work$score)
  temp_plot <- ggplot(work, aes(x = Subtype, y = Rank)) +
    geom_boxplot(alpha = 0.5, width = 0.2) + 
    geom_violin(aes(Subtype, fill = Subtype),
                scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) +
    labs(x = "Subtype", y = "Normalised bact") +
    theme_bw() +
    theme(axis.text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.direction = "horizontal", legend.position = "top") + 
    stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif", method = "wilcox.test")
  filen <- paste0("./Figures/PathSeq/Norm/Blood/", i, ".pdf")
  ggsave(filen, temp_plot)
  c <- c + 1
}

