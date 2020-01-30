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
thousand.folders <- list.dirs(path = "./Data/Counts", full.names = T)
filelist1 <- sapply(thousand.folders[-1], function(x){
  list.files(x, pattern = "htseq.counts.gz$", full.names = T)})
filelist = unlist(filelist1)

     # Read in files and combine
lists <- lapply(filelist, read.delim, header = F)
listsDF <- lists
lists <- listsDF

for (i in names(lists)){
  p <- gsub("./Data/Counts/", "", i)
  colnames(lists[[i]]) <- c("Gene", p)
}

names(lists) <- gsub("./Data/Counts/", "", names(lists))

converter <- read.delim("./Data/Important/gdc_sample_sheet_counts.tsv")
converter$Patient.ID <- gsub("-", ".", converter$Case.ID)
converter$File.Name <- gsub(".htseq.counts.gz", "", converter$File.Name)


lists1 <- lists[names(lists) %in% converter$File.ID]



multi_join <- function(list_of_loaded_data, join_func, ...){
  require("dplyr")
  output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
  return(output)
}

combined_df <- multi_join(lists1, full_join)

temp_df <- combined_df %>% gather(key = "File.ID", value = "Count", -Gene)
converter2 <- droplevels(subset(converter, Sample.Type == "Primary Tumor"))

temp_df1 <- merge(converter2[, c("Patient.ID", "File.ID")], temp_df, by = "File.ID")
library(reshape2)

Counts <- dcast(temp_df1, Gene ~ Patient.ID, sum, value.var = "Count") # This does nothing, but get in correct format.

# library_sizes <- colSums(Counts[!'%in%'(colnames(Counts), "Gene")]) %>% as.data.frame() %>%
#   rownames_to_column(., var = "Patient.ID")
# colnames(library_sizes) <- c("Patient.ID", "Library_size")
# write.csv("./Output/Library_Sizes.csv", x = library_sizes, row.names = F)

Counts_cleaned <- droplevels(Counts[!'%in%'(Counts$Gene, 
                                            c("__alignment_not_unique", "__no_feature",
                                              "__not_aligned", "__too_low_aQual", "__ambiguous")), ])

# Remove the version number from the ensembl Gene IDs (rownames)
Counts_cleaned$Gene <- gsub("\\..*", "", Counts_cleaned$Gene)

# BiocManager::install("biomaRt")
library(biomaRt)
ensembl_DB <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

Gene_Map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
      filters = "ensembl_gene_id", values = Counts_cleaned$Gene, mart = ensembl_DB)


rm(list = setdiff(ls(), c("Gene_Map", "Counts_cleaned"))) # Clean environment
save.image("./R_Data/Counts_clean.RData")

load("./R_Data/Counts_clean.RData")
library(UsefulFunctions)
library(tidyverse)
library(ggpubr)

colnames(Counts_cleaned)[colnames(Counts_cleaned) == "Gene"] <- "ensembl_gene_id"
Counts_SYMBOLS <- merge(Counts_cleaned, Gene_Map, by = "ensembl_gene_id") %>% 
  dplyr:: select(., -ensembl_gene_id) %>% gather(., key = contains("TCGA"), value = "Count")


# The below function adds the counts together to account for the many:1 mapping of ensembl_gene_id to hgnc_symbols
# hgnc_symbols are much easier to work with for my analysis - Neeraj determined all of CIRC in this
# Big assumption, however, is that I assume the various loci are producing proteins with the same function
Counts_totalled <- dcast(Counts_SYMBOLS, hgnc_symbol ~ Patient.ID, sum, value.var = "Count")

Counts <- Counts_totalled[!is.na(Counts_totalled$hgnc_symbol), ]
rownames(Counts) <- NULL

## Should these be used to talk about library sizes?
# Patient_list <- colnames(FPKMs)[!'%in%'(colnames(FPKMs), "Gene")]
# write.table(Patient_list, "./Output/Patient_list.txt", row.names = F)



BiocManager::install("cqn")

FPKM2[!'%in%'(colnames(FPKM2), c("SYMBOL"))] <- log2(FPKM2[!'%in%'(colnames(FPKM2), c("SYMBOL"))] + 1)
library(cqn)

?cqn