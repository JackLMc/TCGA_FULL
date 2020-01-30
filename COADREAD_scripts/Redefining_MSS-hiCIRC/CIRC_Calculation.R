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
load("./R_Data/FPKM_clean.RData")

# Set the seed.
set.seed(123)

## CIRC Score calculation (enrichment of the CIRC gene FPKM) ----
CIRC_IG <- read.csv("./Exploratory_Data/Genesets/CIRC.csv")
CIRC_IG$SYMBOL <- as.factor(CIRC_IG$SYMBOL)

# Label
CIRC_genes <- droplevels(subset(CIRC_IG, CIRC == T)) %>% 
  takegenelevels(.) %>% list()
names(CIRC_genes) <- "CIRC_Genes"

# Calculate Enrichment of CIRC
Enrichment_CIRC <- gsva(FPKM3, CIRC_genes) 

Enrichment_CIRC1 <- Enrichment_CIRC %>% as.data.frame() %>%
  rownames_to_column(., var = "Geneset") %>%
  gather(contains("TCGA"), key = "Patient.ID", value = "Enrich") %>%
  spread(., key = "Geneset", value = "Enrich")

# START CIRC Enrichment --------
# Read Clinical Stuff in ----
Clin_614 <- read.csv("./Output/Clinical_Data_614.csv")
CIRC_clin <- merge(Enrichment_CIRC1, Clin_614, by = "Patient.ID")

# pdf("./Figures/1_Redefinition/Violin Plot of Mean CIRC.pdf", height = 6, width = 6)
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
dcast(CIRC_clin, MSI_STATUS ~ ., length)

# Analysing variance differences ----
# Levene-test (analysis of variance)
CIRC_clin$MSI_STATUS <- as.factor(CIRC_clin$MSI_STATUS)
# car:: leveneTest(CIRC_clin$CIRC_Genes, group = CIRC_clin$MSI_STATUS, center = "median") # Assumes normality

# var.test(CIRC_Genes ~ MSI_STATUS, data = CIRC_clin)
flig_test <- fligner.test(CIRC_clin$CIRC_Genes, g = CIRC_clin$MSI_STATUS)$p.value

Test <- cbind("Fligner Test",
              round(flig_test, 6)) %>% as.data.frame()
colnames(Test) <- c("Method", "P Value")
rownames(Test) <- NULL

# bartlett.test(CIRC_clin$CIRC_Genes, g = CIRC_clin$MSI_STATUS) # Assumes normality

# Coefficient of variance
# install_github("benmarwick/cvequality")
library(cvequality)
cv_test <- with(CIRC_clin, asymptotic_test(CIRC_Genes, MSI_STATUS)) # unequal sample sizes
# cv_test_MSLRT <- with(CIRC_clin, mslr_test(nr = 1e4, CIRC_Genes, MSI_STATUS)) # Equal sample sizes
asym <- cbind("Asymptotic Test", round(cv_test$p_value, 4))
colnames(asym) <- c("Method", "P Value")

rbind(Test, asym)




# Clustering ----
# Partitioning clustering
## Remove uneeded stuff
CIRC_for_Cluster <- droplevels(subset(FPKM2, SYMBOL %in% CIRC_genes$CIRC_Genes))
pca1 <- CIRC_for_Cluster %>% 
  gather(contains("TCGA"), key = "Patient.ID", value = "FPKM") %>%
  spread(key = "SYMBOL", value = "FPKM") %>% 
  merge(Clin_614, by = "Patient.ID") # Merge with cleaned clinical

my_data <- pca1 %>%
  droplevels() %>%
  #select(matches("HLA|MSI|Patient")) %>%
  na.omit()
rownames(my_data) <- NULL

# Standardise the data
my_data1 <- column_to_rownames(my_data, var = "Patient.ID")
my_data2 <- droplevels(subset(my_data1, MSI_STATUS == "MSS"))
my_data3 <- my_data2[, !('%in%'(colnames(my_data2), c("MSI_STATUS")))]

my_data3 <- log(my_data3 + 1)

## Determine optimal number of clusters for kmeans
library(factoextra)
pdf("./Figures/1_Redefinition/Number_kmean_cluster.pdf", width = 6, height = 6)
fviz_nbclust(my_data3, 
             kmeans, method = c("gap_stat"), nboot = 500, nstart = 25, iter.max = 1000)
dev.off()
library(NbClust)
Nb <- NbClust(data = my_data3, diss = NULL, distance = "euclidean",
        min.nc = 2, max.nc = 15, method  = "kmeans")

## Perform Phenograph and kmeans
a <- cbind(my_data3, Phenograph_Clusters = factor(Rphenograph(my_data3[, !('%in%'(names(my_data3), c("Patient.ID", "MSI_STATUS")))])[[2]]$membership), 
           kmeans_Clusters = factor(kmeans(my_data3[, !('%in%'(names(my_data3), c("Patient.ID", "MSI_STATUS")))], centers = 3, iter.max = 1000)$cluster))




## Choosing visualisation method
# Calculate PCs
pca1a <- data.frame(a[, names(a) != "MSI_STATUS" & 
                        names(a) != "Phenograph_Clusters" & 
                        names(a) != "kmeans_Clusters"])

# Make the Patient ID the row names
df <- data.frame(pca1a[, names(pca1a) != "Patient.ID"], row.names = pca1a[, names(pca1a) == "Patient.ID"])

newdata <- log(df + 1)

head(newdata)

### PCA
prin_comp <- prcomp(newdata, scale. = F)

## Plot them 
# ### Phenograph
Pcluster <- a[, "Phenograph_Clusters"]
pdf("./Figures/1_Redefinition/PCA/PhenoG_CIRC.pdf", height = 6, width = 6)
ggbiplot(prin_comp, obs.scale = 1, var.scale = 1,
         groups = Pcluster, circle = T, var.axes = F) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top")
dev.off()

### kmeans
Kcluster <- a[, "kmeans_Clusters"]
pdf("./Figures/1_Redefinition/PCA/kmeans_CIRC.pdf", height = 6, width = 6)
ggbiplot(prin_comp, obs.scale = 1, var.scale = 1,
         groups = as.factor(Nb$Best.partition), circle = T, var.axes = F) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top")
dev.off()


View(newdata)

# # RtSNE
tsne_out <- Rtsne(as.matrix(my_data3))
tsne_dimensions <- as.data.frame(tsne_out$Y)
colnames(tsne_dimensions) <- c("Dim1", "Dim2")

## tSNE plot
for_CIRC <- Enrichment_CIRC1[Enrichment_CIRC1$Patient.ID %in% rownames(df), ]
CIRC_score <- for_CIRC[, "CIRC_Genes"]
# 
# pdf("./Figures/1_Redefinition/tSNE_CIRC_Score.pdf")
# ggplot(tsne_dimensions, aes(x = Dim1, y = Dim2, colour = CIRC_score)) +
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

# ### kmeans
# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# pdf("./Figures/1_Redefinition/tSNE_CIRC_kmeans.pdf")

Nb$Best.partition[Nb$Best.partition == 1] %>% length




ggplot(tsne_dimensions, aes(x = Dim1, y = Dim2, colour = as.factor(Nb$Best.partition))) +
  geom_point(size = 4, alpha = 0.8, pch = 20) +
  # scale_colour_manual(values = c("1" = "#009E73", "2" = "#56B4E9", "3" = "#E69F00",
  #                                "4" = "#CC79A7", "5" = "#0072B2", "6" = "#999999",
  #                                "7" = "#F0E442", "8" = "#D55E00")) +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top")
# dev.off()

# pdf("./Figures/1_Redefinition/tSNE_CIRC_PhenoG.pdf")
# ggplot(tsne_dimensions, aes(x = Dim1, y = Dim2, colour = Pcluster)) +
#   geom_point(size = 4, alpha = 0.8, pch = 20) +
#   scale_colour_manual(values = c("1" = "#009E73", "2" = "#56B4E9", "3" = "#E69F00",
#                                  "4" = "#CC79A7", "5" = "#0072B2", "6" = "#999999",
#                                  "7" = "#F0E442", "8" = "#D55E00")) +
#   theme_bw() +
#   theme(axis.text = element_text(size = 16)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme(legend.direction = "horizontal", legend.position = "top")
# dev.off()

## UMAP
umap_out <- umap(as.matrix(my_data3))
umap_dimensions <- as.data.frame(umap_out$layout)
colnames(umap_dimensions) <- c("Dim1", "Dim2")

for_CIRC <- Enrichment_CIRC1[Enrichment_CIRC1$Patient.ID %in% rownames(df), ]
CIRC_score <- for_CIRC[, "CIRC_Genes"]

pdf("./Figures/1_Redefinition/UMAP_CIRC_Score.pdf")
ggplot(umap_dimensions, aes(x = Dim1, y = Dim2, colour = CIRC_score)) +
  geom_point(size = 4, alpha = 0.8, pch = 20) +
  # scale_colour_manual(values = c("1" = "#009E73", "2" = "#56B4E9", "3" = "#E69F00",
  #                                "4" = "#CC79A7", "5" = "#0072B2", "6" = "#999999",
  #                                "7" = "#F0E442", "8" = "#D55E00")) +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top") +
  scale_colour_gradient(low = "#E3E3E3", high = "#413CFF",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour")
dev.off()

# # kmeans
# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# pdf("./Figures/1_Redefinition/UMAP_CIRC_kmeans.pdf")
ggplot(umap_dimensions, aes(x = Dim1, y = Dim2, colour = Kcluster)) +
  geom_point(size = 4, alpha = 0.8, pch = 20) +
  # scale_colour_manual(values = c("1" = "#009E73", "2" = "#56B4E9", "3" = "#E69F00",
  #                                "4" = "#CC79A7", "5" = "#0072B2", "6" = "#999999",
  #                                "7" = "#F0E442", "8" = "#D55E00")) +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top")
dev.off()

## Phenograph
pdf("./Figures/1_Redefinition/UMAP_CIRC_Pheno.pdf")
ggplot(umap_dimensions, aes(x = Dim1, y = Dim2, colour = as.factor(Nb$Best.partition))) +
  geom_point(size = 4, alpha = 0.8, pch = 20) +
  scale_colour_manual(values = c("1" = "#009E73", "2" = "#000000", "3" = "#E69F00",
                                 "4" = "#CC79A7", "5" = "#0072B2", "6" = "#999999",
                                 "7" = "#F0E442", "8" = "#D55E00")) +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top")
dev.off()





# Count patients
library(knitr)
library(reshape2)
# install.packages("rmarkdown")
# library("rmarkdown")

dcast(a, MSI_STATUS ~ Phenograph_Clusters, length)
# dcast(a, MSI_STATUS ~ kmeans_Clusters, length)
# Calculate CIRC Expression for clusters ----
j <- a[, c("Patient.ID", "kmeans_Clusters", "Phenograph_Clusters", "MSI_STATUS")]
df1a <- merge(Enrichment_CIRC1, j, by = "Patient.ID")
df1 <- droplevels(subset(df1a, MSI_STATUS == "MSS"))

clusters <- Nb$Best.partition %>% as.data.frame() %>% rownames_to_column(., var = "Patient.ID")

colnames(clusters) <- c("Patient.ID", "cluster")

df1 <- merge(Enrichment_CIRC1, clusters, by = "Patient.ID")
df1$cluster <- as.factor(df1$cluster)


# Phenograph
pdf("./Figures/1_Redefinition/CIRC_Pheno_Clusters.pdf", height = 6, width = 6)
ggplot(df1, aes(x = cluster, y = CIRC_Genes)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(aes(cluster, fill = cluster),
              scale = "width", alpha = 0.8) +
  scale_fill_manual(values = c("1" = "#009E73", "2" = "#56B4E9", "3" = "#E69F00",
                               "4" = "#CC79A7", "5" = "#0072B2", "6" = "#999999",
                               "7" = "#F0E442", "8" = "#D55E00")) +
  labs(x = "Phenograph Clusters", y = "CIRC Enrichment Score")+
  theme_bw()+
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top") + 
  stat_compare_means(comparisons = list(c("1", "2")#, #c("1", "3"), c("1", "4"), c("1", "5"), c("1", "6"), c("1", "7"), c("1", "8"),
                                       # c("2", "3"), c("2", "4"), c("2", "5"), c("2", "6"), c("2", "7"), c("2", "8"),
                                       # c("3", "4"), c("3", "5"), c("3", "6"), c("3", "7"), c("3", "8")
                                        # c("4", "5"), c("4", "6"), c("4", "7"), c("4", "8"),
                                        # c("5", "6"), c("5", "7"), c("5", "8"),
                                        # c("6", "7"), c("6", "8"),
                                        # c("7", "8")
  ),
  label = "p.signif", method = "wilcox.test")
dev.off()




# Determine CIRC patients based on CIRC score - Phenograph clustering and tSNE visualisation
## Take patients in cluster 1, 2, or 3 with a CIRC score greater than 0
Patient_pool <- merge(a[, c("Patient.ID", "MSI_STATUS", "Phenograph_Clusters")], Enrichment_CIRC1, by = "Patient.ID")

MSI_H_cutoff <- droplevels(subset(Patient_pool, MSI_STATUS == "MSI-H"))$CIRC_Genes %>% mean()

MSS <- droplevels(subset(Patient_pool, MSI_STATUS == "MSS"))
hiCIRC <- droplevels(subset(MSS, Phenograph_Clusters == "1" & CIRC_Genes >= MSI_H_cutoff | 
                    Phenograph_Clusters == "2" & CIRC_Genes > MSI_H_cutoff | 
                    Phenograph_Clusters == "3" & CIRC_Genes > MSI_H_cutoff))
head(hiCIRC)
poten <- droplevels(subset(MSS, Phenograph_Clusters == "1" | Phenograph_Clusters == "2" | Phenograph_Clusters == "3"))


df1a$Subtype_cut <- ifelse((df1a$Patient.ID %in% hiCIRC$Patient.ID), "MSS-hiCIRC", as.character(df1a$MSI_STATUS))
df1a$Subtype_clust <- ifelse((df1a$Patient.ID %in% poten$Patient.ID), "MSS-hiCIRC", as.character(df1a$MSI_STATUS))


# ## tSNE
# remove_low_CIRC <- df[tsne_dimensions$Dim1 >= (-25) & tsne_dimensions$Dim1 <= -5, ]
# lessthan <- rownames_to_column(remove_low_CIRC, var = "Patient.ID")
# lt <- merge(lessthan, Clin_614, by = "Patient.ID")
# hiCIRC_pats <- droplevels(subset(lt, MSI_STATUS == "MSS"))$Patient.ID
# df1a$Subtype_tSNE <- ifelse((df1a$Patient.ID %in% hiCIRC_pats), "MSS-hiCIRC", as.character(df1a$MSI_STATUS))

# Remove the MSS-hiCIRC patients who are POLE mutants
POLE <- read.csv("Output/POLE_mutants.csv")
df1a$Subtype_cut <- ifelse((df1a$Patient.ID %in% POLE), as.character(df1a$MSI_STATUS), as.character(df1a$Subtype_cut))
df1a$Subtype_clust <- ifelse((df1a$Patient.ID %in% POLE), as.character(df1a$MSI_STATUS), as.character(df1a$Subtype_clust))

ggplot(df1a, aes(x = Subtype_cut, y = CIRC_Genes)) +
  geom_boxplot(alpha = 0.5, width = 0.2) +
  geom_violin(aes(Subtype_cut, fill = Subtype_cut),
              scale = "width", alpha = 0.8) +
  scale_fill_manual(values = cbcols) +
  labs(x = "MSI Status", y = "CIRC Enrichment Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top") +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test")


colnames(df1a)[colnames(df1a) == "Subtype_cut"] <- "Subtype"

# write.csv("./Output/Patient_Subtypes.csv", x = df1a[, c("Patient.ID", "CIRC_Genes", "Subtype")], row.names = F)



hiCIRC_old <- read.csv("./Output/Old_pat_sub.csv")
write.csv("./Output/Old_pat_sub.csv", x = hiCIRC_old, row.names = F)

hiCIRC_oldss <- droplevels(subset(hiCIRC_old, Subtype == "MSS-hiCIRC"))
df1a[hiCIRC_oldss$Patient.ID %in% hiCIRC$Patient.ID, ]

df2 <- merge(df1a, hiCIRC_old[, c("Patient.ID", "Subtype")], by = "Patient.ID")

mismatches <- droplevels(subset(df2, Subtype_cut != Subtype))


# Hierarchical clustering
# library(gplots)
# library(RColorBrewer)
# mypalette <- brewer.pal(11,"RdYlBu")
# morecols <- colorRampPalette(mypalette)
# mycol <- colorpanel(1000,"blue","white","red")
# distCor <- function(x) as.dist(1-cor(t(x)))
# hclustAvg <- function(x) hclust(x, method = "average")
# heatmap.2(as.matrix(df),
#           col = mycol,
#           trace = "none",
#           density.info = "none",
#           cex.main = 1.5,
#           #ColSideColors = col.cell,
#           scale = "column",
#           margin = c(10,5), lhei = c(2,10),
#           hclustfun = hclustAvg,
#           RowSideColors = col.cell1)


# WRITE OUT THE hiCIRC PATIENT SUBTYPES ----
# write.csv("./Output/Patient_Subtypes.csv", x = df1a[, c("Patient.ID", "CIRC_Genes", "Subtype")], row.names = F)
 