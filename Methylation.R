library(tidyverse)
library(UsefulFunctions)
library(ggpubr)

my_comparisons <- list(c("MSS-hiCIRC", "MSI-H"),
                       c("MSS-hiCIRC", "MSS"),
                       c("MSI-H", "MSS"))
cbcols <- c("MSS-hiCIRC" = "#999999",
            "MSI-H" = "#56B4E9",
            "MSS" = "#009E73",
            "MSI-L" = "#E69F00")


COADREAD_pats <- read.delim("./Data/Important/gdc_sample_sheet_FPKM.tsv")$Sample.ID
COADREAD_pats1 <- as.character(COADREAD_pats)
columns_I_want <- append(c("Composite Element REF"), as.character(COADREAD_pats1))

write.table(columns_I_want, file = "./Data/Methylation/COADREAD_pats.txt", row.names = F, col.names = F)


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

important_annotations <- annotation.table[, c("Name", "chr", "Type", "Islands_Name", "UCSC_RefGene_Name",
                                              "UCSC_RefGene_Group")]
colnames(gat_meth_clin)[colnames(gat_meth_clin) == "Composite.Element.REF"] <- "Name"

important_annotations <- droplevels(subset(important_annotations, UCSC_RefGene_Name != ""))
important_annotations$UCSC_RefGene_Name <- as.factor(important_annotations$UCSC_RefGene_Name)
IA <- important_annotations[important_annotations$Name %in% gat_meth_clin$Name, ]
rownames(IA) <- NULL

# rm(important_annotations)
# rm(annotation.table)
# rm(IlluminaHumanMethylation450kanno.ilmn12.hg19)

write.csv(IA, file = "./Data/Methylation/Probe_Annotations.csv", row.names = F)
IA <- read.csv("./Data/Methylation/Probe_Annotations.csv")


### Just Grep for the genes of interest out of the RefGene_Name column
CIITA <- droplevels(IA[grepl("CIITA", IA$UCSC_RefGene_Name), ])

working_on <- merge(gat_meth_clin, CIITA, by = "Name") %>% droplevels()

for(i in levels(working_on$Name)){
  print(paste("Working on the ", i, " of CIITA"))
  work <- droplevels(subset(working_on, Name == i))
  temp_plot <- ggplot(work, aes(x = Subtype, y = beta_val)) +
    geom_boxplot(alpha = 0.5, width = 0.2) + 
    geom_violin(aes(Subtype, fill = Subtype),
                scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) +
    labs(x = "Subtype", y = paste("CIITA ", i, " Beta Value")) +
    theme_bw() +
    theme(axis.text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.direction = "horizontal", legend.position = "top") + 
    stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif", method = "wilcox.test")
  filen <- paste0("CIITA", i, ".pdf")
  ggplot2:: ggsave(filen, plot = temp_plot, device = "pdf",
                   path = "./Figures/Methylation",
                   height = 6, width = 6)}


### TLR4
TLR4 <- droplevels(IA[grepl("TLR4", IA$UCSC_RefGene_Name), ])
working_on <- merge(gat_meth_clin, TLR4, by = "Name") %>% droplevels()

for(i in levels(working_on$Name)){
  print(paste("Working on the ", i, "probe of TLR4"))
  work <- droplevels(subset(working_on, Name == i))
  temp_plot <- ggplot(work, aes(x = Subtype, y = beta_val)) +
    geom_boxplot(alpha = 0.5, width = 0.2) + 
    geom_violin(aes(Subtype, fill = Subtype),
                scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) +
    labs(x = "Subtype", y = paste("TLR4 ", i, " Beta Value")) +
    theme_bw() +
    theme(axis.text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.direction = "horizontal", legend.position = "top") + 
    stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif", method = "wilcox.test")
  filen <- paste0("TLR4_", i, ".pdf")
  ggplot2:: ggsave(filen, plot = temp_plot, device = "pdf",
                   path = "./Figures/Methylation",
                   height = 6, width = 6)}

