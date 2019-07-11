library(UsefulFunctions)
library(tidyverse)
library(ggpubr)

# # This is to remake the lost data...
# load("./R_Data/PathSeq.RData")
# write.table(pathseq, file = "./Data/PathSeq/Raw_Data_20190711/pathseq.txt", sep = "\t", row.names = F)


# Rip out the tissue source codes from the Pathseq patient data
pathseq <- read.delim("./Data/PathSeq/Raw_Data_20190711/pathseq.txt", stringsAsFactors = F)
GDC_convert <- read.delim("./Data/GDC_large_mapping_TCGA.txt")[, c("file_name", "cases.0.submitter_id", "cases.0.samples.0.sample_type",
                                                                   "data_category",
                                                                   "cases.0.samples.0.submitter_id",
                                                                   "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id")]

colnames(GDC_convert) <- c("File.Name", "Patient.ID", "Sample.Type", "Data.Type", "TCGA_Submitter_ID", "TCGA.ID")
GDC_convert$Patient.ID <- gsub("-", ".", GDC_convert$Patient.ID)
GDC_convert <- GDC_convert[GDC_convert$File.Name %in% pathseq$File.Name, ]


pathseq1 <- merge(pathseq, GDC_convert, by = "File.Name") %>% droplevels()
pathseq1$Patient.ID <- as.factor(pathseq1$Patient.ID)
pathseq1$File.Name <- as.factor(pathseq1$File.Name)


## Ripping out data
TCGA_ID <- pathseq1$TCGA.ID[!duplicated(pathseq1$TCGA.ID)] %>% as.character()
Dup <- cbind(TCGA_ID, TCGA_ID) %>% as.data.frame()
colnames(Dup)[2] <- "TCGA_TO_SEP"
gutting_code <- separate(Dup, TCGA_TO_SEP, into = c("Project", "TSS.Code", "Participant",
                                    "Sample.Type", "Portion", "Plate", "ACC.Code"), sep = "-")

# pos <- regexpr(paste0("(TCGA-[[:alnum:]]{2})"), TCGA_ID, perl = T)
# samplestring_unfixed <- substr(TCGA_ID, pos, pos + attributes(pos)[[1]]-1)
# Tissue_Source_Site <- gsub("TCGA-", "", samplestring_unfixed)

We_care <- gutting_code[, c("TCGA_ID", "TSS.Code", "ACC.Code")]

ACC <- read.csv("./Data/Analysis_Centre_Codes.csv")
colnames(ACC)[colnames(ACC) == "Code"] <- "ACC.Code"
colnames(ACC)[colnames(ACC) != "ACC.Code"] <- paste("Analysis", colnames(ACC)[colnames(ACC) != "ACC.Code"], sep = ".")
ACC$ACC.Code <- ifelse((ACC$ACC.Code < 10), paste0("0", as.character(ACC$ACC.Code)), as.character(ACC$ACC.Code))

TSS <- read.csv("./Data/Tissue_Source_Site_Codes.csv")
colnames(TSS)[colnames(TSS) == "Source.Site"] <- "Tissue.Source.Site"

TSS_add <- merge(We_care, TSS, by = "TSS.Code")

Centres <- merge(TSS_add, ACC, by = "ACC.Code") %>% droplevels()
Centres$Patient.ID <- samptopat(Centres$TCGA_ID)
Centres$Patient.ID <- gsub("-", ".", Centres$Patient.ID)
pat_sub <- read.csv("./Output/Patient_Subtypes.csv")

Centre_sub <- merge(Centres, pat_sub[, c("Patient.ID", "Subtype")], by = "Patient.ID") %>% droplevels()

CS <- Centre_sub[, !('%in%'(colnames(Centre_sub), "TCGA_ID"))]
CS1 <- CS[!duplicated(CS), ] 

dcast(CS1, Subtype ~ Source.Site, length)
dcast(CS1, Subtype ~ Analysis.Short.Name, length)
dcast(CS1, Subtype ~ BCR, length)


