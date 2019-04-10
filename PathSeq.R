# PathSeq analysis
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
thousand.folders <- list.dirs(path = "/Volumes/2018/beggsa-tcgacolorectal/download_rest/bacterial_project/results", full.names = T)
filelist1 <- sapply(thousand.folders, function(x){
  list.files(x, pattern = ".pathseq.txt$", full.names = T)})
filelist = unlist(filelist1)



# Read in files and combine
lists <- lapply(filelist, read.delim, header = T)
listsDF <- lists
lists <- listsDF

names(lists) <- filelist
names(lists) <- gsub("/Volumes/2018/beggsa-tcgacolorectal/download_rest/bacterial_project/results/", "", names(lists))
names(lists) <- gsub("_gdc_realn.pathseq.txt", "", names(lists))
names(lists) <- gsub("_hg19_Illumina", "", names(lists))


names(lists)
# Colnames = c(tax_id, taxonomy, type, name,
# kingdom, score, score_normalized, reads, unambiguous, reference_length)
GDC_convert <- read.delim("./Data/GDC_large_mapping_TCGA.txt")
GDC_convert$Patient.ID <- gsub("-", ".", GDC_convert$cases.0.submitter_id)

colnames(GDC_convert)[colnames(GDC_convert) == "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"] <- "file.ID"
head(GDC_convert)

COADREAD_GDC <- droplevels(subset(GDC_convert, cases.0.project.project_id == "TCGA-COAD" | cases.0.project.project_id == "TCGA-READ"))

