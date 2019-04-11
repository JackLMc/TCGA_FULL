library(data.table)
library(tidyverse)
library(UsefulFunctions)

these_pats <- read.delim("./Data/Important/gdc_sample_sheet_FPKM.tsv")

these <- these_pats$Sample.ID
these <- as.character(these)
columns_I_want <- append(c("Composite Element REF"), as.character(these))

write.table(columns_I_want, file = "./Data/Methylation/these_cols.txt", row.names = F, col.names = F)


library(data.table)
methyl <- fread("./Data/Methylation/jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.tsv")

# Need to edit the column names of methyl to match the columns_I_want
## need to rip off the last 12 stuff (e.g. -11D-A41L-05)
these_cols <- colnames(methyl)[colnames(methyl) %in% columns_I_want]


# Need to subset the dataframe by these_cols


pat_sub <- read.csv("./Output/Patient_Subtypes.csv")

