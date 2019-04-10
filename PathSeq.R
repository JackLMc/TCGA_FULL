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
lists <- lapply(filelist, read.delim, header = T, stringsAsFactors = F)
listsDF <- lists
lists <- listsDF
names(lists) <- filelist
names(lists) <- gsub("/Volumes/2018/beggsa-tcgacolorectal/download_rest/bacterial_project/results/", "", names(lists))
names(lists) <- gsub("pathseq.txt", "bam", names(lists))

# Colnames = c(tax_id, taxonomy, type, name,
# kingdom, score, score_normalized, reads, unambiguous, reference_length)
GDC_convert <- read.delim("./Data/PathSeq/gdc_sample_sheet.2019-04-10.tsv")
GDC_convert$Patient.ID <- gsub("-", ".", GDC_convert$Case.ID)

# colnames(GDC_convert)[colnames(GDC_convert) == "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"] <- "file.ID"
GDC_convert_t <- droplevels(subset(GDC_convert, Sample.Type != "Blood Derived Normal" &
                                Sample.Type != "Solid Tissue Normal"))
GDC_convert_n <- droplevels(subset(GDC_convert, Sample.Type == "Blood Derived Normal" |
                                     Sample.Type == "Solid Tissue Normal"))


lists1 <- lists[names(lists) %in% GDC_convert_n$File.Name]
lists2 <- lapply(names(lists1), 
                  function(n, x){
                    x[[n]]$File.Name <- n
                    return (x[[n]])
                  }, lists1)


multi_join <- function(list_of_loaded_data, join_func, ...){
  require("dplyr")
  output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
  return(output)
}

pathseq <- multi_join(lists2, full_join, by = c("tax_id", "taxonomy", "type", "name",
                                                    "kingdom", "score", "score_normalized",
                                                    "reads", "unambiguous", "reference_length", "File.Name"))

pathseq1 <- merge(pathseq, GDC_convert, by = "File.Name") %>% droplevels()
pathseq1$kingdom <- as.factor(pathseq1$kingdom)
pathseq1$name <- as.factor(pathseq1$name)
pathseq1$taxonomy <- as.factor(pathseq1$taxonomy)
pathseq1$Patient.ID <- as.factor(pathseq1$Patient.ID)
pathseq1$type <- as.factor(pathseq1$type)


pathseq2 <- pathseq1[, c("Patient.ID", "kingdom", "type", "name", "score",
                         "score_normalized", "reads", "unambiguous", "reference_length")]

levels(pathseq2$kingdom)
head(pathseq2)

droplevels(subset(pathseq2, type == "no_rank"))



output <- data.frame(stringsAsFactors = F)
c <- 1
for(i in levels(pathseq1$Patient.ID)){
  work <- droplevels(subset(pathseq1, Patient.ID == i))
  these <- nlevels(work$taxonomy)
  output[c, "Patient.ID"] <- i
  output[c, "num_of_tax"] <- these
  c <- c + 1
}

head(output)


pat_sub <- read.csv("./Output/Patient_Subtypes.csv")
