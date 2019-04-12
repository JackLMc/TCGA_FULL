library(tidyverse)
library(UsefulFunctions)

these_pats <- read.delim("./Data/Important/gdc_sample_sheet_FPKM.tsv")$Sample.ID
these <- as.character(these_pats)
columns_I_want <- append(c("Composite Element REF"), as.character(these))

write.table(columns_I_want, file = "./Data/Methylation/these_cols.txt", row.names = F, col.names = F)



methyl <- data.table:: fread("./Data/Methylation/jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.tsv")
methyl <- as.data.frame(methyl)

# Need to edit the column names of methyl to match the columns_I_want
## need to rip off the last 12 stuff (e.g. -11D-A41L-05)
splist <- colnames(methyl)
splist <- as.character(splist)
lslist <- vector(mode = "character", length = length(splist))
samples <- regexpr("^TCGA", splist)
spselect <- samples != -1
lslist[!spselect] <- splist[!spselect]
sp_pos_fix <- regexpr("TCGA[[:punct:]]{1}[[:alnum:]]{2}[[:punct:]]{1}[[:alnum:]]{4}[[:punct:]]{1}[[:alnum:]]{3}",
                      splist[spselect], perl = T)
lslist[spselect] <- substr(splist[spselect],
                           sp_pos_fix,sp_pos_fix+attributes(sp_pos_fix)[[1]]-1)
colnames(methyl) <- lslist
methyl_data <- methyl[, colnames(methyl) %in% columns_I_want]

# Write this out... Stop doing it again.
write.csv(methyl_data, file = "./Data/Methylation/COADREAD_methyl.csv", row.names = F)

methyl_data <- read.csv("./Data/Methylation/COADREAD_methyl.csv")

gathered_meth <- methyl_data %>% gather(contains("TCGA"), key = "Sample.ID", value = "beta_val")
gathered_meth$Patient.ID <- samptopat(gathered_meth$Sample.ID)
gathered_meth$Patient.ID <- gsub("-", ".", gathered_meth$Patient.ID)

pat_sub <- read.csv("./Output/Patient_Subtypes.csv")
gat_meth_clin <- merge(gathered_meth, pat_sub[, c("Patient.ID", "Subtype")], by = "Patient.ID")

# Annotation of probes
# BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
# BiocManager::install("minfi")
# BiocManager::install("readr")

library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
annotation.table <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19) %>% as.data.frame()
annotation.table <- as.data.frame(annotation.table)

important_annotations <- annotation.table[, c("Name", "chr", "Type", "Islands_Name", "UCSC_RefGene_Name",
                                              "UCSC_RefGene_Group")]
colnames(gat_meth_clin)[colnames(gat_meth_clin) == "Composite.Element.REF"] <- "Name"

might_be_big <- merge(gat_meth_clin, important_annotations, by = "Name")


head(might_be_big)
