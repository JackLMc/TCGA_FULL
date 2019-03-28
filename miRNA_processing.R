library(UsefulFunctions)
library(tidyverse)
library(ggpubr)

my_comparisons <- list(c("MSS-hiCIRC", "MSI-H"),
                       c("MSS-hiCIRC", "MSS"),
                       c("MSI-H", "MSS"))
cbcols <- c("MSS-hiCIRC" = "#999999",
            "MSI-H" = "#56B4E9",
            "MSS" = "#009E73")

miRNA <- read.csv("./Data/miRNA/miRNA.csv")
miRNA1 <-miRNA %>% gather(contains("TCGA"), key = "Patient.ID", value = "FPKM")
colnames(miRNA1)[colnames(miRNA1) == "Patient.ID"] <- "Sample.ID"

pat_sub <- read.csv("./Output/Patient_Subtypes.csv")
 
miRNA2 <- merge(miRNA1[, c("Patient.ID", "Genes", "FPKM")], pat_sub[, c("Patient.ID", "Subtype")], by = "Patient.ID")


miRNA21s <- c("hsa-miR-21-3p", "hsa-miR-21-5p")
GOI <- droplevels(subset(miRNA2, Genes == "hsa-miR-21-5p"))
GOI$Rank <- rank(GOI$FPKM)
ggplot(GOI, aes(x = Subtype, y = Rank)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(aes(Subtype, fill = Subtype),
              scale = "width", alpha = 0.8) +
  scale_fill_manual(values = cbcols) +
  labs(x = "MSI Status", y = paste("Rank transformed", "FPKM", sep = " ")) +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top") + 
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test")





miRNA1$Patient.ID <- samptopat(miRNA1$Sample.ID)
miRNA1$Patient.ID <- gsub("-", ".", miRNA1$Patient.ID)

levels(miRNA1$Genes)[grepl("21", levels(miRNA1$Genes))]

mir_Int <- read.csv("./Exploratory_Data/miRNA_interactions.csv")
CIITA <- droplevels(subset(mir_Int, Target.Gene == "CIITA"))




try <- miRNA1[miRNA1$Gene %in% CIITA$miRNA,]

head(try)

pat_sub <- read.csv("Output/Patient_Subtypes.csv")

these <- merge(pat_sub, try[, c("Patient.ID", "Genes", "FPKM")]) %>% droplevels()

head(these)
levels(these$Genes)

for(i in levels(these$Genes)){
  print(i)
  work <- droplevels(subset(these, Genes == i))
  work$Rank <- rank(work$FPKM)
  temp_plot <- ggplot(work, aes(x = Subtype, y = Rank)) +
    geom_boxplot(alpha = 0.5, width = 0.2) + 
    geom_violin(aes(Subtype, fill = Subtype),
                scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) +
    labs(x = "MSI Status", y = paste("Rank transformed", i, "FPKM", sep = " ")) +
    theme_bw() +
    theme(axis.text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.direction = "horizontal", legend.position = "top") + 
    stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif", method = "wilcox.test")
  filen <- paste0(i, ".pdf")
  ggplot2:: ggsave(filen, plot = temp_plot, device = "pdf",
                   path = "./Figures/miRNAs",
                   height = 6, width = 6)
}

head(mir_Int)
genesofI <- droplevels(subset(mir_Int, miRNA == "hsa-miR-335-5p"& Support.Type != "Functional MTI (Weak)"))


levels(genesofI$Target.Gene)


droplevels(subset(mir_Int, Target.Gene == "HLA-DRB5" ))
           