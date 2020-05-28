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
