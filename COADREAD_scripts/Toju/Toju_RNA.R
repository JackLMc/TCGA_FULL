## Packages
required <- c("tidyverse", "ggpubr", "ggbiplot", "devtools", "gplots", "UsefulFunctions")
for (lib in required){
  if (!require(lib, character.only = T)){
    install.packages(lib)
    suppressMessages(library(lib, character.only = T, quietly = T))
  }
}

# Comparisons and Colours
my_comparisons <- list(c("VD1.CD27HI", "VD1.CD27LO"), c("VD1.CD27HI", "CD8.EMRA"), c("VD1.CD27HI", "CD8.Naive"), c("VD1.CD27HI", "VD2"),
                       c("VD1.CD27LO", "CD8.EMRA"), c("VD1.CD27LO", "CD8.Naive"), c("VD1.CD27LO", "VD2"),
                       c("CD8.EMRA", "CD8.Naive"), c("CD8.EMRA", "VD2"), c("CD8.Naive", "VD2"))

# Colours
cbcols <- c("VD1.CD27LO" = "#999999", "CD8.EMRA" = "#56B4E9",
            "CD8.Naive" = "#E69F00", "VD1.CD27HI" = "#009E73",
            "VD2" = "#CC79A7")

# source("https://bioconductor.org/biocLite.R")
# biocLite("Rsubread", dependencies = T)
library(Rsubread)

# Run after having run the terminal protocol
# bam.files <- get_sorbam("/Volumes/ResearchData/Willcox Group/Jack/GD_RNA_Comb/BAM/")
bam.files <- list.files("/Volumes/JackMcMurray/Toju_data/bam/", pattern = ".bam$", full.names = T)

# patientID <- gsub(".bam", "", bam.files)
# fastq.files <- list.files("/Volumes/JackMcMurray/Toju_data/fastq/", pattern = ".fastq.gz$")
# sample.lsit <- read.delim("/Volumes/JackMcMurray/Toju_data/fastq/sample.list.txt", header = F)
# 
# 
# edit.fastq <- gsub("_S[[:digit:]]{1,}_R1_001.fastq.gz$", "", fastq.files)
# edit.fastq <- gsub("_[[:digit:]]{1}$", "", edit.fastq)
# fastq.files[duplicated(edit.fastq)]
# fastq.files[grep("S377079", fastq.files)]

## Summary of the proportion of reads that are mapped
# props <- propmapped(files = bam.files)
# props

## Counting
### Contains inbuilt annotation for hg19 genome assembly
fc <- featureCounts(bam.files, annot.inbuilt = "hg38")
Counts <- as.data.frame(fc$counts)

# Gain genelengths
## Create a dataframe to allow merging
Counts1 <- tibble:: rownames_to_column(Counts, var = "GeneID")
Counts1$GeneID <- as.factor(Counts1$GeneID)
annotation <- as.data.frame(fc$annotation)
ann <- annotation[ , which(names(annotation) %in% c("GeneID", "Length"))]
ann$GeneID <- as.factor(ann$GeneID)

Counts2 <- droplevels(merge(Counts1, ann, by = "GeneID"))

# Gather the Pats
head(Counts2)

Counts3 <- Counts2 %>% gather(contains(".bam"), key = "ID", value = "Count")
Counts3$ID <- as.factor(Counts3$ID)

head(Counts3)

# loop to make individual dataframes in a list
listofsamp <- list()

c <- 1
for(i in levels(Counts3$ID)){
  print(i)
  work <- droplevels(subset(Counts3, ID == i))
  listofsamp[[i]] <- work
  c <- c + 1
}

# Get rid of "ID" in the list
LoS <- lapply(listofsamp, function(x) x[!(names(x) %in% c("ID"))])

# write Csvs for each member of the list
setwd("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/6_TCGA_Full/TCGA_Full/Data/Toju/Counts")
for(i in names(LoS)){
  write.table(LoS[[i]], paste0(i,".txt"), sep = "\t", row.names = F, quote = F)
}
save.image(file = "../Toju_counts.RData")

setwd("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/6_TCGA_Full/TCGA_Full/")




### CQN Normalisation
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
filelist1 <- list.files("./Data/Toju/Counts", pattern = ".txt$", full.names = T)
filelist = unlist(filelist1)

# Read in files and combine
lists <- lapply(filelist, read.delim, header = T)
names(lists) <- gsub("./Data/Toju/Counts/", "", filelist)
for(i in names(lists)){
  lists[[i]]$Patient.ID <- gsub(".bam.txt", "", i)}



multi_join <- function(list_of_loaded_data, join_func, ...){
  require("dplyr")
  output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
  return(output)
}

combined_df <- multi_join(lists, full_join)
combined_df$Patient.ID <- as.factor(combined_df$Patient.ID)


levels(combined_df$Patient.ID)
lists[["S377079.1.bam.txt"]]


combined_df <- cdf


# Check that each patient has one file id
combined_df$Pat_short <- gsub("\\.S*[[:digit:]]{1,}$", "", combined_df$Patient.ID) %>% as.factor()


try <- droplevels(subset(combined_df, GeneID == "1")) %>% droplevels()
nlevels(try$Patient.ID)
nlevels(try$Pat_short)

these_pats <- try[duplicated(try$Pat_short), ]$Pat_short %>% droplevels() %>% levels()

multiples <- combined_df[combined_df$Pat_short %in% these_pats, ] %>% droplevels()
multiples$Pat_short <- as.factor(multiples$Pat_short)

# Decide how to subset the patients which have two runs - smallest library size = removed
decide <- data.frame()
c <- 1
for(i in levels(multiples$Patient.ID)){
  print(i)
  work <- droplevels(subset(multiples, Patient.ID == i))
  
  library_size <- sum(work$Count)
  
  decide[c, "Pat_short"] <- levels(work$Pat_short)
  decide[c, "Patient.ID"] <- i
  decide[c, "Library"] <- library_size
  c <- c + 1}


decide <- droplevels(decide)


remove_these <- data.frame()
decide$Patient.ID <- as.factor(decide$Patient.ID)
decide$Pat_short <- as.factor(decide$Pat_short)

# levels(decide$Patient.ID)

c <- 1
for(i in levels(decide$Pat_short)){
  print(i)
  work <- droplevels(subset(decide, Pat_short == i))
  samp_to_remove <- with(work, Patient.ID[Library == min(Library)])
  tjis <- as.character(samp_to_remove)
  remove_these[c, "Patient.ID"] <- tjis
  remove_these[c, "Pat_short"] <- i
  c <- c + 1
}

temp_df1 <- combined_df[!'%in%'(combined_df$Patient.ID, remove_these), ] %>% 
  droplevels()

nlevels(temp_df1$Pat_short)

#### temp_df1 = each patient has a unique Patient.ID, can remove the old Patient.ID column, and rename Pat_short to Patient.ID
temp_df1$Patient.ID <- temp_df1$Pat_short
temp_df2 <- temp_df1[, !'%in%'(colnames(temp_df1), "Pat_short")]


# Start the reshaping....
library(reshape2)

Counts <- dcast(temp_df2, GeneID ~ Patient.ID, sum, value.var = "Count") # This does nothing, but get it in the correct format

# BiocManager::install("biomaRt")
library(biomaRt)
ensembl_DB <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

Gene_Map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id", values = Counts_cleaned$Gene, ensembl_DB) # Only finds 56,543 of the 60,483 genes


rm(list = setdiff(ls(), c("Gene_Map", "Counts_cleaned", "ensembl_DB"))) # Clean environment
# save.image("./R_Data/Counts_clean1.RData")

load("./R_Data/Counts_clean1.RData")
library(UsefulFunctions)
library(tidyverse)
library(ggpubr)
library(reshape2)

colnames(Counts_cleaned)[colnames(Counts_cleaned) == "Gene"] <- "ensembl_gene_id"

Counts_merged <- merge(Counts_cleaned, Gene_Map, by = "ensembl_gene_id")  
Counts_SYMBOLS <- Counts_merged[, !'%in%'(colnames(Counts_merged), "ensembl_gene_id")]

Counts_long <- Counts_SYMBOLS %>% 
  gather(-hgnc_symbol, value = "Count", key = "Patient.ID")

# The below function adds the counts together to account for the many:1 mapping of ensembl_gene_id to hgnc_symbols
# hgnc_symbols are much easier to work with for my analysis - Neeraj determined all of CIRC in this
# Big assumption, however, is that I assume the various loci are producing proteins with the same function
Counts_totalled <- dcast(Counts_long, hgnc_symbol ~ Patient.ID, sum, value.var = "Count")

Counts <- Counts_totalled[!is.na(Counts_totalled$hgnc_symbol), ]
Counts$hgnc_symbol <- as.factor(Counts$hgnc_symbol)

Counts <- droplevels(subset(Counts, hgnc_symbol != "")) # Remove genes which don't have a hgnc_symbol
row.names(Counts) <- NULL
Counts <- column_to_rownames(Counts, var = "hgnc_symbol")

# library_sizes <- colSums(Counts[!'%in%'(colnames(Counts), "hgnc_symbol")]) %>% as.data.frame() %>%
#   rownames_to_column(., var = "Patient.ID")
# colnames(library_sizes) <- c("Patient.ID", "Library_size")
# CPM_scaling <- min(library_sizes$Library_size)/1000000 # Calcualte the cpm scaling value based on the smallest library size
# cpm_scale <- 2/CPM_scaling # minimum 1 count
# 
# library(edgeR)
# cpms <- cpm(Counts)
# 
# keep.exprs <- rowSums(cpms>cpm_scale)>=54
# 
# 
# x <- cpms[keep.exprs,]
# dim(x)
# CIRC_IG <- read.csv("./Exploratory_Data/Genesets/CIRC.csv")
# CIRC_IG$SYMBOL <- as.factor(CIRC_IG$SYMBOL)
# 
# CIRC <- CIRC_IG[CIRC_IG$CIRC,]$SYMBOL
# 
# x[rownames(x) %in% CIRC, ] %>% dim()

Genes_to_keep <- Counts[rowSums(Counts) != 0, ] %>% rownames() # Keep all genes which are expressed in at least 1 patient


CIRC_IG <- read.csv("./Exploratory_Data/Genesets/CIRC.csv") # Check that the CIRC genes are in this data
CIRC_IG$SYMBOL <- as.factor(CIRC_IG$SYMBOL)
CIRC <- CIRC_IG[CIRC_IG$CIRC, ]$SYMBOL %>% droplevels() %>% levels()

Counts_above_0 <- Counts[rownames(Counts) %in% Genes_to_keep, ]
# Counts_above_0[rownames(Counts_above_0) %in% CIRC, ] %>% dim()

rm(list = setdiff(ls(), c("Gene_Map", "Counts_above_0", "ensembl_DB"))) # Clean environment
save.image("./R_Data/Counts_clean2.RData")

# Normalisation
# BiocManager::install("cqn")
load("./R_Data/Counts_clean2.RData")
library(cqn)
Gene_Map1 <- getBM(attributes = c("hgnc_symbol", "percentage_gene_gc_content", "start_position", "end_position"),
                   filters = "hgnc_symbol", values = rownames(Counts_above_0), ensembl_DB)
libraries <- colSums(Counts_above_0[!'%in%'(colnames(Counts_above_0), "hgnc_symbol")])

duplicate_info <- Gene_Map1$hgnc_symbol[duplicated(Gene_Map1$hgnc_symbol)] %>% unique() # These genes have multiple numbers for GC or Start or end positions
length(duplicate_info)

Gene_Map1$hgnc_symbol <- as.factor(Gene_Map1$hgnc_symbol)
Gene_Map2 <- data.frame()
c <- 1
for(i in levels(Gene_Map1$hgnc_symbol)){
  print(i)
  working <- droplevels(subset(Gene_Map1, hgnc_symbol == i))
  Average_GC <- mean(working$percentage_gene_gc_content)
  Average_start <- mean(working$start_position)
  Average_end <- mean(working$end_position)
  Gene_Map2[c, "hgnc_symbol"] <- i
  Gene_Map2[c, "Average_percentage_gc_content"] <- Average_GC
  Gene_Map2[c, "Average_start_position"] <- Average_start
  Gene_Map2[c, "Average_end_position"] <- Average_end
  c <- c + 1
}

Gene_Map2$length <- Gene_Map2$Average_end_position - Gene_Map2$Average_start_position
row.names(Gene_Map2) <- NULL

Gene_Map3 <- column_to_rownames(Gene_Map2, var = "hgnc_symbol")
stopifnot(all(rownames(Counts_above_0) == rownames(Gene_Map3)))
stopifnot(colnames(Counts_above_0) == names(libraries))

# Normalisation
cqn_Counts <- cqn(Counts_above_0, lengths = Gene_Map3$length, x = Gene_Map3$Average_percentage_gc_content,
                  sizeFactors = libraries, verbose = T)


## Look at systematic effects, n argument refers to the effect, 1 is always the covariate specified by the x argument above (gc content)
## whereas 2 is lengths
par(mfrow = c(1, 2))
cqnplot(cqn_Counts, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7))
cqnplot(cqn_Counts, n = 2, xlab = "length", lty = 1, ylim = c(1,7))


## Normalised expression values
Counts_cqn <- cqn_Counts$y + cqn_Counts$offset

rm(list = setdiff(ls(), c("Gene_Map3", "ensembl_DB", "Counts_cqn", "cqn_Counts", "libraries"))) # Clean environment


head(Counts_cqn) # These values are on a log2 scale.
write.table("./Output/Patient_list.txt", x = as.factor(colnames(Counts_cqn)), row.names = F)
save.image("./R_Data/Counts_clean.RData")

