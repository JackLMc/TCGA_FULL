# A script to find MSS-hiCIRC patients
## Packages
library(UsefulFunctions)
library(tidyverse)
library(ggpubr)
library(reshape2)
library(GSVA)

# BiocManager::install("devtools")
# devtools:: install_github("https://github.com/JinmiaoChenLab/Rphenograph")
library(Rphenograph)

library(ggbiplot)

# install.packages("Rtsne")
library(Rtsne)

# install.packages("umap")
library(umap)

my_comparisons <- list(c("MSS-hiCIRC", "MSI-H"),
                       c("MSS-hiCIRC", "MSS"),
                       c("MSI-H", "MSS"))
cbcols <- c("MSS-hiCIRC" = "#999999",
            "MSI-H" = "#56B4E9",
            "MSS" = "#009E73")

## Load required RData objects
load("./R_Data/Counts_clean.RData")

# Set the seed.
set.seed(123)

## CIRC Score calculation (enrichment of the CIRC gene CQN) ----
SigGen <- read.csv("./Exploratory_Data/Genesets/Signature_Genesets.csv")
CIRC_IG <- droplevels(subset(SigGen, Signature == "CIRC"))
CIRC_IG$HUGO.symbols <- as.factor(CIRC_IG$HUGO.symbols)

# Label
CIRC_genes <- droplevels(subset(CIRC_IG, Signature == "CIRC"))$HUGO.symbols %>%
  levels() %>% list()
names(CIRC_genes) <- "CIRC_Genes"

# Calculate Enrichment of CIRC
Enrichment_CIRC <- gsva(Counts_cqn, CIRC_genes) 
Enrichment_CIRC1 <- Enrichment_CIRC %>% as.data.frame() %>%
  rownames_to_column(., var = "Geneset") %>%
  gather(contains("TCGA"), key = "Patient.ID", value = "Enrich") %>%
  spread(., key = "Geneset", value = "Enrich")

# START MSS-hi-CIRC Caclulation --------
# Read Clinical Stuff in ----
Clin_542 <- read.csv("./Output/Clinical_Data_542.csv")[, c("Patient.ID", "MSI_STATUS")]
CIRC_clin <- merge(Enrichment_CIRC1, Clin_542, by = "Patient.ID")

shapiro.test(Enrichment_CIRC1$CIRC_Genes) # Sig different from normal distribution

# pdf("./Figures/1_Redefinition/CIRC_MSS_MSI.pdf", height = 6, width = 6)
ggplot(CIRC_clin, aes(x = MSI_STATUS, y = CIRC_Genes)) +
  geom_boxplot(alpha = 0.5, width = 0.2) +
  geom_violin(aes(MSI_STATUS, fill = MSI_STATUS),
              scale = "width", alpha = 0.8) +
  scale_fill_manual(values = cbcols) +
  labs(x = "MSI Status", y = "CIRC Enrichment Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top") +
  stat_compare_means(comparisons = list((c("MSI-H", "MSS"))),
                     label = "p.signif", method = "wilcox.test")
dev.off()

# Number of patients per group
count(CIRC_clin, vars = c("MSI_STATUS"))

# Analysing variance differences ----
# Levene-test (analysis of variance)
CIRC_clin$MSI_STATUS <- as.factor(CIRC_clin$MSI_STATUS)
fligner.test(CIRC_clin$CIRC_Genes, g = CIRC_clin$MSI_STATUS)$p.value # Statistical test of variance, variance is the measure of how far numbers differ from the mean

MSI_H <- droplevels(subset(CIRC_clin, MSI_STATUS == "MSI-H"))
MSS <- droplevels(subset(CIRC_clin, MSI_STATUS == "MSS"))
ks.test(MSI_H$CIRC_Genes, MSS$CIRC_Genes, na.rm = T)

# devtools:: install_github("tpepler/nonpar")
library(nonpar)
cucconi.test(MSI_H$CIRC_Genes, MSS$CIRC_Genes, method = c("permutation", "bootstrap"))

# Clustering ----
# Partitioning clustering - use MSI_H patients to find highest clusters for cutoff of CIRC score
## Remove uneeded stuff
CIRC_for_Cluster <- Counts_cqn[rownames(Counts_cqn) %in% CIRC_genes$CIRC_Genes, ]
CIRC_for_Cluster <- rownames_to_column(as.data.frame(Counts_cqn), var = "SYMBOL")

pca1 <- CIRC_for_Cluster %>%
  tidyr:: gather(contains("TCGA"), key = "Patient.ID", value = "CQN")  %>%
  spread(key = "SYMBOL", value = "CQN") %>% 
  merge(Clin_542, by = "Patient.ID") # Merge with cleaned clinical

my_data <- pca1 %>%
  droplevels() %>%
  #select(matches("HLA|MSI|Patient")) %>%
  na.omit()
rownames(my_data) <- NULL

my_data1 <- column_to_rownames(my_data, var = "Patient.ID")
my_data2 <- my_data1

## Subset the data 
IG_genes <- droplevels(subset(CIRC_IG, IG == T)) %>% 
  takegenelevels(.) %>% list()
names(IG_genes) <- "IG_Genes"
these_genes <- append(IG_genes$IG_Genes, "MSI_STATUS")
my_data3 <- my_data2[, colnames(my_data2) %in% these_genes]

MSI_H <- droplevels(subset(my_data3, MSI_STATUS == "MSI-H"))

## Determine optimal number of clusters for kmeans
library(NbClust)
Nb <- NbClust(data = MSI_H[, !('%in%'(colnames(MSI_H), c("MSI_STATUS")))], diss = NULL, distance = "euclidean",
              min.nc = 2, max.nc = 15, method  = "kmeans")


## Perform Phenograph and kmeans
a <- cbind(MSI_H, Phenograph_Clusters = factor(Rphenograph(MSI_H[, !('%in%'(colnames(MSI_H), c("MSI_STATUS")))])[[2]]$membership), 
           kmeans_Clusters = factor(kmeans(MSI_H[, !('%in%'(colnames(MSI_H), c("MSI_STATUS")))], centers = 2, iter.max = 1000)$cluster))

# Clusters with highest expression
## Count patients
library(reshape2)
dcast(a, MSI_STATUS ~ Phenograph_Clusters, length)

df1a <- rownames_to_column(a, var = "Patient.ID") %>% 
  dplyr:: select(., c("Patient.ID", "kmeans_Clusters", "MSI_STATUS", "Phenograph_Clusters")) %>%
  merge(., Enrichment_CIRC1, by = "Patient.ID")

# Phenograph
pdf("./Figures/1_Redefinition/CIRC_Pheno_Clusters_MSI.pdf", height = 6, width = 6)
ggplot(df1a, aes(x = Phenograph_Clusters, y = CIRC_Genes)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(aes(Phenograph_Clusters, fill = Phenograph_Clusters),
              scale = "width", alpha = 0.8) +
  scale_fill_manual(values = c("1" = "#E69F00",
                               "2" = "#CC79A7", "3" = "#0072B2", "6" = "#999999",
                               "7" = "#F0E442", "8" = "#D55E00")) +
  labs(x = "Phenograph Clusters", y = "CIRC Enrichment Score")+
  theme_bw()+
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top") + 
  stat_compare_means(comparisons = list(c("1", "2"), c("1", "3"), c("2", "3")),
                     label = "p.signif", method = "wilcox.test")
dev.off()


## Choosing visualisation method
# Calculate PCs 
pca1a <- data.frame(a[, names(a) != "MSI_STATUS" & 
                        names(a) != "Phenograph_Clusters" & 
                        names(a) != "kmeans_Clusters"])

# Make the Patient ID the row names
df <- MSI_H[, !'%in%'(colnames(MSI_H), "MSI_STATUS")]

# df[, colnames(df) == "CIRC"]
head(df)[, 1:10]

### PCA
prin_comp <- prcomp(df, scale. = F)
for_CIRC <- Enrichment_CIRC1[Enrichment_CIRC1$Patient.ID %in% rownames(df), ]
CIRC_score <- for_CIRC[, "CIRC_Genes"]

### Phenograph
Pcluster <- a[, "Phenograph_Clusters"]
# pdf("./Figures/1_Redefinition/PCA_CIRC_Score.pdf", height = 6, width = 6)
# ggbiplot(prin_comp, obs.scale = 1, var.scale = 1,
#          groups = CIRC_score, circle = T, var.axes = F) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme(legend.direction = "horizontal", legend.position = "top") +
#   scale_colour_gradient(low = "#E3E3E3", high = "#413CFF",
#                         space = "Lab", na.value = "grey50", guide = "colourbar",
#                         aesthetics = "colour")
# dev.off()

# pdf("./Figures/1_Redefinition/PCA_CIRC_Pheno_MSI.pdf", height = 6, width = 6)
# ggbiplot(prin_comp, obs.scale = 1, var.scale = 1,
#          groups = Pcluster, circle = T, var.axes = F) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme(legend.direction = "horizontal", legend.position = "top") +
#   scale_colour_manual(values = c( "1" = "#E69F00",
#                                  "2" = "#CC79A7", "3" = "#0072B2", "6" = "#999999",
#                                  "7" = "#F0E442", "8" = "#D55E00"))
# dev.off()

## UMAP
umap_out <- umap(as.matrix(MSI_H[, !('%in%'(colnames(MSI_H), c("MSI_STATUS", "CIRC")))]))
umap_dimensions <- as.data.frame(umap_out$layout)
colnames(umap_dimensions) <- c("Dim1", "Dim2")

for_CIRC <- Enrichment_CIRC1[Enrichment_CIRC1$Patient.ID %in% rownames(df), ]
CIRC_score <- for_CIRC[, "CIRC_Genes"]

# pdf("./Figures/1_Redefinition/UMAP_CIRC_Score.pdf")
# ggplot(umap_dimensions, aes(x = Dim1, y = Dim2, colour = CIRC_score)) +
#   geom_point(size = 4, alpha = 0.8, pch = 20) +
#   # scale_colour_manual(values = c("1" = "#009E73", "2" = "#56B4E9", "3" = "#E69F00",
#   #                                "4" = "#CC79A7", "5" = "#0072B2", "6" = "#999999",
#   #                                "7" = "#F0E442", "8" = "#D55E00")) +
#   theme_bw() +
#   theme(axis.text = element_text(size = 16)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme(legend.direction = "horizontal", legend.position = "top") +
#   scale_colour_gradient(low = "#E3E3E3", high = "#413CFF",
#                         space = "Lab", na.value = "grey50", guide = "colourbar",
#                         aesthetics = "colour")
# dev.off()


## UMAP
pdf("./Figures/1_Redefinition/UMAP_CIRC_Pheno.pdf")
ggplot(umap_dimensions, aes(x = Dim1, y = Dim2, colour = Pcluster)) +
  geom_point(size = 8, alpha = 0.8, pch = 20) +
  scale_colour_manual(values = c("1" = "#E69F00",
                                 "2" = "#CC79A7", "3" = "#0072B2", "6" = "#999999",
                                 "7" = "#F0E442", "8" = "#D55E00")) +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top")
dev.off()



# Determine CIRC patients based on CIRC score - Phenograph clustering and UMAP visualisation
## The number of MSI-H patients in cluster 1 or 2 is 
count(df1a, vars = c("Phenograph_Clusters"))

## 27 + 32 = 59
## Total patients = 59/86 = 0.6860465 (69%)
# Take the cutoff as 31% quantile
cutoff <- unname(quantile(df1a$CIRC_Genes, 1 - (59/86)))

Seg_hiCIRC <- merge(Enrichment_CIRC1, Clin_542, by = "Patient.ID")
hiCIRC <- droplevels(subset(Seg_hiCIRC, MSI_STATUS == "MSS" & CIRC_Genes >= cutoff))
Patient_pool <- Seg_hiCIRC
Patient_pool$Subtype <- ifelse((Patient_pool$Patient.ID %in% hiCIRC$Patient.ID), "MSS-hiCIRC", as.character(Patient_pool$MSI_STATUS))

# Remove the MSS-hiCIRC patients who are POLE mutants
POLE <- read.csv("Output/POLE_mutants.csv")
Patient_pool$Subtype <- ifelse((Patient_pool$Patient.ID %in% POLE$.), as.character(Patient_pool$MSI_STATUS), as.character(Patient_pool$Subtype))

droplevels(subset(Patient_pool, Subtype == "MSS-hiCIRC")) %>% dim
count(Patient_pool, vars = c("Subtype"))


pdf("./Figures/1_Redefinition/CIRC_Score_hiCIRC.pdf")
ggplot(Patient_pool, aes(x = Subtype, y = CIRC_Genes)) +
  geom_boxplot(alpha = 0.5, width = 0.2) +
  geom_violin(aes(Subtype, fill = Subtype),
              scale = "width", alpha = 0.8) +
  scale_fill_manual(values = cbcols) +
  labs(x = "MSI Status", y = "CIRC Enrichment Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top") +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test")
dev.off()

# WRITE OUT THE hiCIRC PATIENT SUBTYPES ----
write.csv("./Output/Patient_Subtypes_09_03.csv", x = Patient_pool[, c("Patient.ID", "CIRC_Genes", "Subtype")], row.names = F)

Patient_pool <- read.csv("./Output/Patient_Subtypes_09_03.csv")

#### END ####