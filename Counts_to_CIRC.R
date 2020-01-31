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

Counts_cleaned <- droplevels(Counts[!'%in%'(Counts$Gene, 
                                            c("__alignment_not_unique", "__no_feature",
                                              "__not_aligned", "__too_low_aQual", "__ambiguous")), ])

# Remove the version number from the ensembl Gene IDs (rownames)
Counts_cleaned$Gene <- gsub("\\..*", "", Counts_cleaned$Gene)

# BiocManager::install("biomaRt")
library(biomaRt)
ensembl_DB <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

Gene_Map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
      filters = "ensembl_gene_id", values = Counts_cleaned$Gene, ensembl_DB) # Only finds 56,543 of the 60,483 genes


rm(list = setdiff(ls(), c("Gene_Map", "Counts_cleaned", "ensembl_DB"))) # Clean environment
save.image("./R_Data/Counts_clean.RData")

load("./R_Data/Counts_clean.RData")
library(UsefulFunctions)
library(tidyverse)
library(ggpubr)

colnames(Counts_cleaned)[colnames(Counts_cleaned) == "Gene"] <- "ensembl_gene_id"

Counts_merged <- merge(Counts_cleaned, Gene_Map, by = "ensembl_gene_id")  
Counts_SYMBOLS <- Counts_merged[, !'%in%'(colnames(Counts_merged), "ensembl_gene_id")]

Counts_long <- Counts_SYMBOLS %>% gather(-hgnc_symbol, value = "Count", key = "Patient.ID")

# The below function adds the counts together to account for the many:1 mapping of ensembl_gene_id to hgnc_symbols
# hgnc_symbols are much easier to work with for my analysis - Neeraj determined all of CIRC in this
# Big assumption, however, is that I assume the various loci are producing proteins with the same function
Counts_totalled <- dcast(Counts_long, hgnc_symbol ~ Patient.ID, sum, value.var = "Count")

Counts <- Counts_totalled[!is.na(Counts_totalled$hgnc_symbol), ]
Counts$hgnc_symbol <- as.factor(Counts$hgnc_symbol)

Counts <- droplevels(subset(Counts, hgnc_symbol != "")) # Remove genes which don't have a hgnc_symbol
row.names(Counts) <- NULL
Counts <- column_to_rownames(Counts, var = "hgnc_symbol")

library_sizes <- colSums(Counts[!'%in%'(colnames(Counts), "hgnc_symbol")]) %>% as.data.frame() %>%
  rownames_to_column(., var = "Patient.ID")
colnames(library_sizes) <- c("Patient.ID", "Library_size")
CPM_scaling <- min(library_sizes$Library_size)/1000000 # Calcualte the cpm scaling value based on the smallest library size
cpm_scale <- 6.5/CPM_scaling # 6.5 count = minimum I believe (this is from edgeR user manual)

cpms <- cpm(Counts)
average <- rowMeans(cpms)

Genes_to_keep <- average[average >= cpm_scale] %>% names() # Only keep the genes which have an average cpm above the scaling value
# This is how it was done in the GSVA paper

Counts_above_6.5 <- Counts[rownames(Counts) %in% Genes_to_keep, ]

rm(list = setdiff(ls(), c("Gene_Map", "Counts_above_6.5", "ensembl_DB"))) # Clean environment

# Normalisation
# BiocManager::install("cqn")
library(cqn)
Gene_Map1 <- getBM(attributes = c("hgnc_symbol", "percentage_gene_gc_content", "start_position", "end_position"),
                   filters = "hgnc_symbol", values = rownames(Counts_above_6.5), ensembl_DB)
libraries <- colSums(Counts_above_6.5[!'%in%'(colnames(Counts_above_6.5), "hgnc_symbol")])

duplicate_info <- Gene_Map1$hgnc_symbol[duplicated(Gene_Map1$hgnc_symbol)] # These genes have multiple numbers for GC or Start or end positions
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
stopifnot(all(rownames(Counts_above_6.5) == rownames(Gene_Map3)))
stopifnot(colnames(Counts_above_6.5) == names(libraries))

# Normalisation
cqn_Counts <- cqn(Counts_above_6.5, lengths = Gene_Map3$length, x = Gene_Map3$Average_percentage_gc_content,
                  sizeFactors = libraries, verbose = T)


## Look at systematic effects, n argument refers to the effect, 1 is always the covariate specified by the x argument above (gc content)
## whereas 2 is lengths
par(mfrow = c(1, 2))
cqnplot(cqn_Counts, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7))
cqnplot(cqn_Counts, n = 2, xlab = "length", lty = 1, ylim = c(1,7))


## Normalised expression values
Counts_cqn <- cqn_Counts$y + cqn_Counts$offset

rm(list = setdiff(ls(), c("Gene_Map3", "Counts_above_6.5", "ensembl_DB", "Counts_cqn", "cqn_Counts", "libraries"))) # Clean environment


head(Counts_cqn) # These values are on a log2 scale.

save.image("./R_Data/Counts_clean.RData")
