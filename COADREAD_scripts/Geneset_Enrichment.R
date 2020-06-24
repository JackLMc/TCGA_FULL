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

# source("Clinical.R") # Run to gain the clinical dataframe that's in Output (Clin_540)

load("./R_Data/Counts_clean.RData")

pat_sub <- read.csv("./Output/Patient_Subtypes_09_03.csv")
library(reshape2)
dcast(pat_sub, Subtype ~., length)

Ccqn <- Counts_cqn

# MSS_pats <- droplevels(subset(pat_sub, Subtype != "MSI-H"))
# Ccqn <- Counts_cqn[, colnames(Counts_cqn) %in% MSS_pats$Patient.ID]


# GENESET Interrogation ----
# BiocManager::install("GSVA")
library(GSVA)

## Correlate with the CIRC
SigGen <- read.csv("./Exploratory_Data/Genesets/Signature_Genesets.csv")
SigGen <- factorthese(SigGen, c("Signature", "HUGO.symbols"))

head(SigGen)

SigGen_List <- list()
c <- 1
for(i in levels(SigGen$Signature)){
  print(i)
  work <- droplevels(subset(SigGen, Signature == i))
  Genes <- levels(work$HUGO.symbols)
  SigGen_List[[i]] <- Genes
  c <- c + 1
}

Enrichment_SigGen <- gsva(Ccqn, SigGen_List)
Enrichment_SigGen1 <- Enrichment_SigGen %>% as.data.frame() %>%
  rownames_to_column(., var = "Geneset") %>%
  gather(contains("TCGA"), key = "Patient.ID", value = "Enrich") %>%
  spread(., key = "Geneset", value = "Enrich")

Enrichment_SigGen1$Patient.ID <- as.factor(Enrichment_SigGen1$Patient.ID)
Enrich <- merge(pat_sub, Enrichment_SigGen1, by = "Patient.ID")

Enrich1 <- Enrich %>% gather(key = "Parameter", value = "Enrichment", -CIRC_Genes, -Patient.ID, -Subtype)
Enrich1$Parameter <- as.factor(Enrich1$Parameter)

## Plot pearson correlation
for(i in levels(Enrich1$Parameter)){
  print(i)
  work <- droplevels(subset(Enrich1, Parameter == i))
  temp_plot <- ggplot(work, aes(y = Enrichment, x = CIRC_Genes))+
    geom_point(alpha = 0.8, size = 4, colour = "slategray") +
    labs(x = "CIRC enrichment score", y = paste(i, "enrichment score")) +
    theme_bw() +
    # geom_text(aes(x = -0.3, y = .75, label = lm_eqn(lm(CIRC_Genes ~ Enrichment, work))), parse = T) +
    # scale_color_manual(values = cbcols) +
    geom_smooth(method = "lm", se = F) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.position = "none") +
    stat_cor()
  filen <- paste0(i, ".pdf")
  ggplot2:: ggsave(filen, plot = temp_plot, device = "pdf", path = "./Figures/Gene_Sets/Pearson",
                   height = 6, width = 6)
}

# Cell Types (immunome and Castro [Th17])
CTGenesets <- read.csv("./Exploratory_Data/Genesets/Investigative_Genesets.csv", stringsAsFactors = T)
Genesets <- deduplicate(CTGenesets)
geneset_list <- list()
for(i in levels(Genesets$Cell.population)){
  print(i)
  work <- droplevels(subset(Genesets, Cell.population == i))
  genes <- levels(work$HUGO.symbols)
  geneset_list[[i]] <- genes
}

Enrichments <- gsva(Ccqn, geneset_list) %>% as.data.frame() %>%
  rownames_to_column(., "Geneset") %>% gather(contains("TCGA"), key = "Patient.ID", value = "Enrichment") %>%
  merge(., pat_sub, by = "Patient.ID")

Enrichments$Geneset <- as.factor(Enrichments$Geneset)

## Pearson
No_MSI <- droplevels(subset(Enrichments, Subtype != "MSI-H")) # Remove these patients from pearson correlations
for(i in levels(No_MSI$Geneset)){
  print(i)
  work <- droplevels(subset(No_MSI, Geneset == i))
  temp_plot <- ggplot(work, aes(y = Enrichment, x = CIRC_Genes))+
    geom_point(alpha = 0.8, size = 4, colour = "slategray") +
    labs(x = "CIRC enrichment score", y = paste(i, "enrichment score")) +
    theme_bw() +
    # geom_text(aes(x = -0.3, y = .75, label = lm_eqn(lm(CIRC_Genes ~ Enrichment, work))), parse = T) +
    # scale_color_manual(values = cbcols) +
    geom_smooth(method = "lm", se = F) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.position = "none") +
    stat_cor()
  filen <- paste0(i, ".pdf")
  ggplot2:: ggsave(filen, plot = temp_plot, device = "pdf", path = "./Figures/Gene_Sets/Pearson",
                   height = 6, width = 6)
}

## Enrichment for Subtype
for(i in levels(Enrichments$Geneset)){
  print(i)
  work <- droplevels(subset(Enrichments, Geneset == i))
  temp_plot <- ggplot(work, aes(x = Subtype, y = Enrichment)) +
    geom_boxplot(alpha = 0.5, width = 0.2) + 
    geom_violin(aes(Subtype, fill = Subtype),
                scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) +
    labs(x = "MSI Status", y = paste(i, "enrichment score")) +
    theme_bw() +
    theme(axis.text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.direction = "horizontal", legend.position = "top") + 
    stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif", method = "wilcox.test")
  filen <- paste0(i, ".pdf")
  ggplot2:: ggsave(filen, plot = temp_plot, device = "pdf",
                   path = "./Figures/Gene_Sets/Enrichment",
                   height = 6, width = 6)}

# GO TERMS
## Reactive Oxygen Species
ROS <- read.csv("./Exploratory_Data/Genesets/GO_term_summary_20190320_151206.csv")

ROS_list <- list()
c <- 1
for(i in levels(ROS$Annotated.Term)){
  print(i)
  work <- droplevels(subset(ROS, Annotated.Term == i))
  Genes <- toupper(levels(work$Symbol))
  ROS_list[[i]] <- Genes
  c <- c + 1
}

Enrichment_book <- gsva(Ccqn, ROS_list)
Enrichment_book1 <- Enrichment_book %>% as.data.frame() %>%
  rownames_to_column(., var = "Geneset") %>%
  gather(contains("TCGA"), key = "Patient.ID", value = "Enrich") %>%
  spread(., key = "Geneset", value = "Enrich")


Enrichment_book1$Patient.ID <- as.factor(Enrichment_book1$Patient.ID)
Enrich <- merge(pat_sub, Enrichment_book1, by = "Patient.ID")

Enrich1 <- Enrich %>% gather(key = "Parameter", value = "Enrichment", -CIRC_Genes, -Patient.ID, -Subtype)
Enrich1$Parameter <- as.factor(Enrich1$Parameter)

for(i in levels(Enrich1$Parameter)){
  print(i)
  work <- droplevels(subset(Enrich1, Parameter == i))
  temp_plot <- ggplot(work, aes(x = Subtype, y = Enrichment)) +
    geom_boxplot(alpha = 0.5, width = 0.2) + 
    geom_violin(aes(Subtype, fill = Subtype),
                scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) +
    labs(x = "Subtype", y = paste(i, "enrichment score")) +
    theme_bw() +
    theme(axis.text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.direction = "horizontal", legend.position = "top") + 
    stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif", method = "wilcox.test")
  filen <- paste0(i, ".pdf")
  ggplot2:: ggsave(filen, plot = temp_plot, device = "pdf",
                   path = "./Figures/Gene_Sets/Enrichment/ROS",
                   height = 6, width = 6)}

# Fatty acid metabolism
FAM  <- read.csv("./Exploratory_Data/Genesets/GO_term_summary_20190603_065405.csv")

FAM_list <- list()
c <- 1
for(i in levels(FAM$Annotated.Term)){
  print(i)
  work <- droplevels(subset(FAM, Annotated.Term == i))
  Genes <- toupper(levels(work$Symbol))
  FAM_list[[i]] <- Genes
  c <- c + 1
}

library(GSVA)
Enrichment_book <- gsva(Ccqn, FAM_list)
Enrichment_book1 <- Enrichment_book %>% as.data.frame() %>%
  rownames_to_column(., var = "Geneset") %>%
  gather(contains("TCGA"), key = "Patient.ID", value = "Enrich") %>%
  spread(., key = "Geneset", value = "Enrich")


Enrichment_book1$Patient.ID <- as.factor(Enrichment_book1$Patient.ID)
Enrich <- merge(pat_sub, Enrichment_book1, by = "Patient.ID")

Enrich1 <- Enrich %>% gather(key = "Parameter", value = "Enrichment", -CIRC_Genes, -Patient.ID, -Subtype)
Enrich1$Parameter <- as.factor(Enrich1$Parameter)

for(i in levels(Enrich1$Parameter)){
  print(i)
  work <- droplevels(subset(Enrich1, Parameter == i))
  temp_plot <- ggplot(work, aes(x = Subtype, y = Enrichment)) +
    geom_boxplot(alpha = 0.5, width = 0.2) + 
    geom_violin(aes(Subtype, fill = Subtype),
                scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) +
    labs(x = "Subtype", y = paste(i, "enrichment score")) +
    theme_bw() +
    theme(axis.text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.direction = "horizontal", legend.position = "top") + 
    stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif", method = "wilcox.test")
  filen <- paste0(i, ".pdf")
  ggplot2:: ggsave(filen, plot = temp_plot, device = "pdf",
                   path = "./Figures/Gene_Sets/Enrichment/FAM",
                   height = 6, width = 6)}


# New Th17
# Fatty acid metabolism
FAM  <- read.csv("./Exploratory_Data/Genesets/Castro_collated.csv")

FAM_list <- list()
c <- 1
for(i in levels(FAM$CellType)){
  print(i)
  work <- droplevels(subset(FAM, CellType == i))
  Genes <- toupper(levels(work$Hugo_Symbol))
  FAM_list[[i]] <- Genes
  c <- c + 1
}

library(GSVA)
Enrichment_book <- gsva(Ccqn, FAM_list)
Enrichment_book1 <- Enrichment_book %>% as.data.frame() %>%
  rownames_to_column(., var = "Geneset") %>%
  gather(contains("TCGA"), key = "Patient.ID", value = "Enrich") %>%
  spread(., key = "Geneset", value = "Enrich")


Enrichment_book1$Patient.ID <- as.factor(Enrichment_book1$Patient.ID)
Enrich <- merge(pat_sub, Enrichment_book1, by = "Patient.ID")

Enrich1 <- Enrich %>% gather(key = "Parameter", value = "Enrichment", -CIRC_Genes, -Patient.ID, -Subtype)
Enrich1$Parameter <- as.factor(Enrich1$Parameter)
head(Enrichment_book1)
for(i in levels(Enrich1$Parameter)){
  print(i)
  work <- droplevels(subset(Enrich1, Parameter == i))
  temp_plot <- ggplot(work, aes(x = Subtype, y = Enrichment)) +
    geom_boxplot(alpha = 0.5, width = 0.2) + 
    geom_violin(aes(Subtype, fill = Subtype),
                scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) +
    labs(x = "Subtype", y = paste(i, "enrichment score")) +
    theme_bw() +
    theme(axis.text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.direction = "horizontal", legend.position = "top") + 
    stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif", method = "wilcox.test")
  filen <- paste0(i, ".pdf")
  ggplot2:: ggsave(filen, plot = temp_plot, device = "pdf",
                   path = "./Figures/Gene_Sets/Enrichment/",
                   height = 6, width = 6)}



## Bespoke
book_list <- list()
book_list[["SAAs"]] <- c("TLR4", "LY96"#,
                         #"TIRAP"#, "MYD88"
                         # ,"IRAK4", "IRAK1",
                         # "TAB1", "TAB2", "MAP3K7",
                         # "MAP2K3", "MAP2K6", "MAP2K4", "MAP2K7",
                         # "MAPK14", "MAPK8"
)

a_list <- list()
a_list[["test"]] <- c("AGER", "APP", "ATF1", "ATF2", "BIRC2", "BIRC3", "BPI", "BTK",
                      "BTRC", 	"CASP8", 	"CD14", 	"CD180", 	"CD36", 	"CHUK", 	"CREB1", 	"CUL1",
                      "DHX9", 	"DNM1", 	"DNM2", 	"DNM3", 	"DUSP3", 	"DUSP4", 	"DUSP6", 	"DUSP7",
                      "ECSIT", 	"ELK1", 	"FADD", 	"FBXW11", 	"FOS", 	"HMGB1", 	"IKBKB", 	"IKBKE",
                      "IKBKG", 	"IRAK1", 	"IRAK2", 	"IRAK3", 	"IRAK4", 	"IRF3", 	"IRF7", 	"ITGAM",
                      "ITGB2", 	"JUN", 	"LBP", 	"LY86", 	"LY96", 	"MAP2K1", 	"MAP2K3", 	"MAP2K4",
                      "MAP2K6", 	"MAP2K7", 	"MAP3K1", 	"MAP3K7", 	"MAP3K8", 	"MAPK1", 	"MAPK10", 	"MAPK11",
                      "MAPK14", 	"MAPK3", 	"MAPK7", 	"MAPK8", 	"MAPK9", 	"MAPKAPK2", 	"MAPKAPK3", 	"MEF2A",
                      "MEF2C", 	"MIR6502", 	"MIR718", 	"MYD88", 	"NFKB1", 	"NFKB2", 	"NFKBIA", 	"NFKBIB",
                      "NOD1", 	"NOD2", 	"PELI1", 	"PELI2", 	"PELI3", 	"PLCG2", 	"PPP2CA", 	"PPP2CB",
                      "PPP2R1A", 	"PPP2R1B", 	"PPP2R5D", 	"PTPN11", 	"PTPN4", 	"RELA", 	"RIPK1", 	"RIPK2",
                      "RIPK3", 	"RPS27A", 	"RPS6KA1", 	"RPS6KA2", 	"RPS6KA3", 	"RPS6KA5", 	"S100A12", 	"S100B",
                      "SAA1", 	"SARM1", 	"SIGIRR", 	"SKP1", 	"SOCS1", 	"TAB1", 	"TAB2", 	"TAB3",
                      "TANK", 	"TBK1", 	"TICAM1", 	"TICAM2", 	"TIRAP", 	"TLR1", 	"TLR2", 	"TLR3",
                      "TLR4", 	"TLR6", 	"TNIP2", 	"TRAF3", 	"TRAF6", 	"UBA52", 	"UBB", 	"UBC",
                      "UBE2D1", 	"UBE2D2", 	"UBE2D3", 	"UBE2N", 	"UBE2V1", 	"VRK3")

library(GSVA)
Enrichment_book <- gsva(Ccqn, book_list)
Enrichment_book1 <- Enrichment_book %>% as.data.frame() %>%
  rownames_to_column(., var = "Geneset") %>%
  gather(contains("TCGA"), key = "Patient.ID", value = "Enrich") %>%
  spread(., key = "Geneset", value = "Enrich")

Enrichment_book1$Patient.ID <- as.factor(Enrichment_book1$Patient.ID)
Enrich <- merge(pat_sub, Enrichment_book1, by = "Patient.ID")

Enrich1 <- Enrich %>% gather(key = "Parameter", value = "Enrichment", -CIRC_Genes, -Patient.ID, -Subtype)
Enrich1$Parameter <- as.factor(Enrich1$Parameter)

pdf("./Figures/Gene_Sets/ClassI_response.pdf")
ggplot(Enrich1, aes(x = Subtype, y = Enrichment)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(aes(Subtype, fill = Subtype),
              scale = "width", alpha = 0.8) +
  scale_fill_manual(values = cbcols) +
  labs(x = "MSI Status", y = "Enrichment of Class I Genes") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top") + 
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test")
dev.off()






# INDIVIDUAL GENES INTERROGATION ----
TCcqn <- Ccqn %>% as.data.frame() %>%
  rownames_to_column(., var = "Gene") %>%
  gather(contains("TCGA"), key = "Patient.ID", value = "Enrich") %>%
  spread(., key = "Gene", value = "Enrich") %>%
  merge(., pat_sub, by = "Patient.ID")

head(TCcqn)[, 1:10]
LongCQN <- TCcqn %>% gather(-Patient.ID, -Subtype, key = "Gene", value = "CQN")

genes_of_interest <- c("IL6", "IL1B", "IL23A", "TGFB1",
                       "CCL2", "CCL5", "CXCL10", "CCL20",
                       "CCR6", "TLR4", "TLR2", "CIITA",
                       "RORC", "IL17A", "IL23R",
                       "CDC42", "FCAR", "IGHA1", "EHMT2", "OX40",
                       "OX40L", "CD1D", "CCL17", "TLR3")



GOI <- droplevels(LongCQN[LongCQN$Gene %in% genes_of_interest, ]) %>%
  .[, c("Patient.ID", "Gene", "Subtype", "CQN")]
GOI$Gene <- as.factor(GOI$Gene)

## Remove for Pearson correlation
# GOI1 <- droplevels(subset(GOI, Subtype != "MSI-H"))
# for(i in levels(GOI1$SYMBOL)){
#   print(i)
#   work <- droplevels(subset(GOI1, SYMBOL == i))
#   work$Logged <- log2(work$FPKM + 1)
#   temp_plot <- ggplot(work, aes(y = Logged, x = CIRC_Genes))+
#     geom_point(alpha = 0.8, size = 4, colour = "slategray") +
#     labs(x = "CIRC_Genes Enrichment", y = paste0("Log2(", i, " + 1)")) +
#     theme_bw() +
#     # geom_text(aes(x = -0.3, y = .75, label = lm_eqn(lm(CIRC_Genes ~ Enrichment, work))), parse = T) +
#     # scale_color_manual(values = cbcols) +
#     geom_smooth(method = "lm", se = F) +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#     theme(legend.position = "none") +
#     stat_cor()
#   filen <- paste0(i, ".pdf")
#   ggplot2:: ggsave(filen, plot = temp_plot, device = "pdf",
#                    path = "./Figures/Genes_of_interest/Correlations",
#                    height = 6, width = 6)
# }

## Compare across Subtypes
for(i in 1:length(genes_of_interest)){
  gene <- genes_of_interest[i]
  print(gene)
  GOI <- droplevels(subset(LongCQN, Gene == gene))
  temp_plot <- ggplot(GOI, aes(x = Subtype, y = CQN)) +
    geom_boxplot(alpha = 0.5, width = 0.2) + 
    geom_violin(aes(Subtype, fill = Subtype),
                scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) +
    labs(x = "MSI Status", y = paste(gene, "CQN Count", sep = " ")) +
    theme_bw() +
    theme(axis.text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.direction = "horizontal", legend.position = "top") + 
    stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif", method = "wilcox.test")
  filen <- paste0(gene, ".pdf")
  ggplot2:: ggsave(filen, plot = temp_plot, device = "pdf",
                   path = "./Figures/Genes_of_interest",
                   height = 6, width = 6)}

## Bespoke genes
head(LongCQN)

LongCQN$Gene[grepl("CYBB", LongCQN$Gene)] # Check whether your gene exists in the dataset
GOI <- droplevels(subset(LongCQN, Gene == "CD14")) 


ggplot(GOI, aes(x = Subtype, y = CQN)) +
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



# THIS PART OF THE SCRIPT IS TO INVESTIGATE WHATEVER GARY WANTS ----
work <- droplevels(subset(LongCQN, Gene == "CYBB" | Gene == "CD14"))
work1 <- spread(work[, c("Patient.ID", "Subtype", "Gene", "CQN")], key = "Gene", value = "CQN")

ggplot(work1, aes(y = CD14, x = CYBB))+
  geom_point(alpha = 0.8, size = 4, colour = "slategray") +
  labs(x = "cqn(CYBB)", y = "cqn(CD14)") +
  theme_bw() +
  # scale_color_manual(values = cbcols) +
  geom_smooth(method = "lm", se = F) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  stat_cor()


filen <- paste0(i, ".pdf")
ggplot2:: ggsave("CIITA_EHMT2_correl.pdf", plot = temp_plot, device = "pdf",
                 path = "./Figures/Genes_of_interest/Correlations",
                 height = 6, width = 6)


