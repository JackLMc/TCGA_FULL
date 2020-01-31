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
CIRC_IG <- read.csv("./Exploratory_Data/Genesets/CIRC.csv")
CIRC_IG$SYMBOL <- as.factor(CIRC_IG$SYMBOL)

# Label
CIRC_genes <- droplevels(subset(CIRC_IG, CIRC == T)) %>% 
  takegenelevels(.) %>% list()
names(CIRC_genes) <- "CIRC_Genes"

# Calculate Enrichment of CIRC
Enrichment_CIRC <- gsva(Counts_cqn, CIRC_genes) 

Enrichment_CIRC1 <- Enrichment_CIRC %>% as.data.frame() %>%
  rownames_to_column(., var = "Geneset") %>%
  gather(contains("TCGA"), key = "Patient.ID", value = "Enrich") %>%
  spread(., key = "Geneset", value = "Enrich")

# START CIRC Enrichment --------
# Read Clinical Stuff in ----
Clin_614 <- read.csv("./Old_Output/Clinical_Data_614.csv") # Now it's in old_output
CIRC_clin <- merge(Enrichment_CIRC1, Clin_614, by = "Patient.ID")

shapiro.test(Enrichment_CIRC1$CIRC_Genes) # Sig different from normal distributiomn

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

CIRC_for_Cluster <- Counts_cqn[rownames(Counts_cqn) %in% CIRC_genes$CIRC_Genes, ]
CIRC_for_Cluster <- rownames_to_column(as.data.frame(CIRC_for_Cluster), var = "SYMBOL")
head(CIRC_for_Cluster)

pca1 <- CIRC_for_Cluster %>%
  tidyr:: gather(contains("TCGA"), key = "Patient.ID", value = "CQN")  %>%
  spread(key = "SYMBOL", value = "CQN") %>% 
  merge(Clin_614, by = "Patient.ID") # Merge with cleaned clinical

head(pca1)[, 1:10]

my_data <- pca1 %>%
  droplevels() %>%
  #select(matches("HLA|MSI|Patient")) %>%
  na.omit()
rownames(my_data) <- NULL

# Standardise the data
my_data1 <- column_to_rownames(my_data, var = "Patient.ID")
my_data2 <- my_data1
# my_data2 <- droplevels(subset(my_data1, MSI_STATUS == "MSS"))

## Determine optimal number of clusters for kmeans
library(NbClust)
Nb <- NbClust(data = my_data2[, !('%in%'(colnames(my_data2), c("MSI_STATUS")))], diss = NULL, distance = "euclidean",
              min.nc = 2, max.nc = 15, method  = "kmeans")

## Perform Phenograph and kmeans
a <- cbind(my_data2, Phenograph_Clusters = factor(Rphenograph(my_data2[, !('%in%'(colnames(my_data2), c("MSI_STATUS")))])[[2]]$membership), 
           kmeans_Clusters = factor(kmeans(my_data2[, !('%in%'(colnames(my_data2), c("MSI_STATUS")))], centers = 3, iter.max = 1000)$cluster))

# Clusters with highest expression
## Count patients
library(reshape2)
dcast(a, MSI_STATUS ~ Phenograph_Clusters, length)

df1a <- rownames_to_column(a, var = "Patient.ID") %>% 
  dplyr:: select(., c("Patient.ID", "kmeans_Clusters", "MSI_STATUS", "Phenograph_Clusters")) %>%
  merge(., Enrichment_CIRC1, by = "Patient.ID")

MSS <- droplevels(subset(df1a, MSI_STATUS == "MSS"))

  
# Phenograph
pdf("./Figures/1_Redefinition/CIRC_Pheno_Clusters.pdf", height = 6, width = 6)
ggplot(MSS, aes(x = Phenograph_Clusters, y = CIRC_Genes)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(aes(Phenograph_Clusters, fill = Phenograph_Clusters),
              scale = "width", alpha = 0.8) +
  scale_fill_manual(values = c("1" = "#009E73", "2" = "#56B4E9", "3" = "#E69F00",
                               "4" = "#CC79A7", "5" = "#0072B2", "6" = "#999999",
                               "7" = "#F0E442", "8" = "#D55E00")) +
  labs(x = "Phenograph Clusters", y = "CIRC Enrichment Score")+
  theme_bw()+
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top") + 
  stat_compare_means(comparisons = list(c("1", "2"), c("1", "3"), c("1", "4"), c("1", "5"), c("1", "6"),# c("1", "7"), c("1", "8"),
                                        c("2", "3"), c("2", "4"), c("2", "5"), c("2", "6"),# c("2", "7"), c("2", "8"),
                                        c("3", "4"), c("3", "5"), c("3", "6"),# c("3", "7"), c("3", "8")
                                        c("4", "5"), c("4", "6"),# c("4", "7"), c("4", "8"),
                                        c("5", "6")#, c("5", "7"), c("5", "8"),
                                        #c("6", "7"), c("6", "8"),
                                        #c("7", "8")
  ),
  label = "p.signif", method = "wilcox.test")
dev.off()


## Choosing visualisation method
# Calculate PCs 
pca1a <- data.frame(a[, names(a) != "MSI_STATUS" & 
                        names(a) != "Phenograph_Clusters" & 
                        names(a) != "kmeans_Clusters"])

# Make the Patient ID the row names
df <- pca1a[, !'%in%'(colnames(pca1a), "CIRC")]


head(df)

### PCA
prin_comp <- prcomp(df, scale. = F)
for_CIRC <- Enrichment_CIRC1[Enrichment_CIRC1$Patient.ID %in% rownames(df), ]
CIRC_score <- for_CIRC[, "CIRC_Genes"]
### Phenograph
Pcluster <- a[, "Phenograph_Clusters"]
pdf("./Figures/1_Redefinition/PCA_CIRC_Score.pdf", height = 6, width = 6)
ggbiplot(prin_comp, obs.scale = 1, var.scale = 1,
         groups = CIRC_score, circle = T, var.axes = F) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top") +
  scale_colour_gradient(low = "#E3E3E3", high = "#413CFF",
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        aesthetics = "colour")
dev.off()

pdf("./Figures/1_Redefinition/PCA_CIRC_Pheno.pdf", height = 6, width = 6)
ggbiplot(prin_comp, obs.scale = 1, var.scale = 1,
         groups = Pcluster, circle = T, var.axes = F) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top") +
  scale_colour_manual(values = c("1" = "#009E73", "2" = "#56B4E9", "3" = "#E69F00",
                                 "4" = "#CC79A7", "5" = "#0072B2", "6" = "#999999",
                                 "7" = "#F0E442", "8" = "#D55E00"))
dev.off()


# RtSNE
head(my_data2)
tsne_out <- Rtsne(as.matrix(my_data2[, !('%in%'(colnames(my_data2), c("MSI_STATUS", "CIRC")))]))
tsne_dimensions <- as.data.frame(tsne_out$Y)
colnames(tsne_dimensions) <- c("Dim1", "Dim2")

## tSNE plot
pdf("./Figures/1_Redefinition/tSNE_CIRC_Score.pdf")
ggplot(tsne_dimensions, aes(x = Dim1, y = Dim2, colour = CIRC_score)) +
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


pdf("./Figures/1_Redefinition/tSNE_CIRC_PhenoG.pdf")
ggplot(tsne_dimensions, aes(x = Dim1, y = Dim2, colour = Pcluster)) +
  geom_point(size = 4, alpha = 0.8, pch = 20) +
  scale_colour_manual(values = c("1" = "#009E73", "2" = "#56B4E9", "3" = "#E69F00",
                                 "4" = "#CC79A7", "5" = "#0072B2", "6" = "#999999",
                                 "7" = "#F0E442", "8" = "#D55E00")) +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top")
dev.off()

## UMAP
umap_out <- umap(as.matrix(my_data2[, !('%in%'(colnames(my_data2), c("MSI_STATUS", "CIRC")))]))
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


## Phenograph
pdf("./Figures/1_Redefinition/UMAP_CIRC_Pheno.pdf")
ggplot(umap_dimensions, aes(x = Dim1, y = Dim2, colour = Pcluster)) +
  geom_point(size = 4, alpha = 0.8, pch = 20) +
  scale_colour_manual(values = c("1" = "#009E73", "2" = "#000000", "3" = "#E69F00",
                                 "4" = "#CC79A7", "5" = "#0072B2", "6" = "#999999",
                                 "7" = "#F0E442", "8" = "#D55E00")) +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top")
dev.off()





# Determine CIRC patients based on CIRC score - Phenograph clustering and tSNE visualisation
## Take patients in cluster 3 or 4 with a CIRC score greater than 0
hiCIRC <- droplevels(subset(MSS, Phenograph_Clusters == "4" & CIRC_Genes >= 0 | 
                              Phenograph_Clusters == "3" & CIRC_Genes >= 0))

Patient_pool <- df1a

Patient_pool$Subtype <- ifelse((Patient_pool$Patient.ID %in% hiCIRC$Patient.ID), "MSS-hiCIRC", as.character(Patient_pool$MSI_STATUS))

# ## tSNE
# remove_low_CIRC <- df[tsne_dimensions$Dim1 >= (-25) & tsne_dimensions$Dim1 <= -5, ]
# lessthan <- rownames_to_column(remove_low_CIRC, var = "Patient.ID")
# lt <- merge(lessthan, Clin_614, by = "Patient.ID")
# hiCIRC_pats <- droplevels(subset(lt, MSI_STATUS == "MSS"))$Patient.ID
# Patient_pool$Subtype_tSNE <- ifelse((Patient_pool$Patient.ID %in% hiCIRC_pats), "MSS-hiCIRC", as.character(Patient_pool$MSI_STATUS))

# Remove the MSS-hiCIRC patients who are POLE mutants
POLE <- read.csv("Output/POLE_mutants.csv")
Patient_pool$Subtype <- ifelse((Patient_pool$Patient.ID %in% POLE$.), as.character(Patient_pool$MSI_STATUS), as.character(Patient_pool$Subtype))

droplevels(subset(Patient_pool, Subtype == "MSS-hiCIRC")) %>% dim



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


# write.csv("./Output/Patient_Subtypes.csv", x = df1a[, c("Patient.ID", "CIRC_Genes", "Subtype")], row.names = F)



hiCIRC_old <- read.csv("./Output/Old_pat_sub.csv")
write.csv("./Output/Old_pat_sub.csv", x = hiCIRC_old, row.names = F)

hiCIRC_oldss <- droplevels(subset(hiCIRC_old, Subtype == "MSS-hiCIRC"))
Patient_pool[hiCIRC_oldss$Patient.ID %in% hiCIRC$Patient.ID, ]

df2 <- merge(Patient_pool, hiCIRC_old[, c("Patient.ID", "Subtype")], by = "Patient.ID")
head(df2)
mismatches <- droplevels(subset(df2, Subtype.y != Subtype.x))


head(Patient_pool)

# WRITE OUT THE hiCIRC PATIENT SUBTYPES ----
write.csv("./Output/Patient_Subtypes_30_01.csv", x = Patient_pool[, c("Patient.ID", "CIRC_Genes", "Subtype")], row.names = F)


## Investigate the PC
# PCA interrogation (Remember which set you're investigating with, both inclusion of NN and exclusion are prin_comp)
## Find the Eigenvalues
library(factoextra)
fviz_screeplot(prin_comp, ncp = 10, choice = "eigenvalue")


eig <- (prin_comp$sdev)^2

## Variances in percentage
variance <- eig*100/sum(eig)

## Cumulative variances
cumvar <- cumsum(variance)

Eigenvalues.All <- data.frame(eig = eig, variance = variance,
                              cumvariance = cumvar)
# writeCsvO(Eigenvalues.All)

# Store the variances
var <- get_pca_var(prin_comp)

## Find the correlation between variables and principal components
loadings <- prin_comp$rotation
loads <- as.data.frame(loadings)
load1 <- tibble:: rownames_to_column(loads, "Parameter")
load1$Parameter <- as.factor(load1$Parameter)

# Biggest contributions to a given PC
load1[order(load1$PC1, decreasing = T)[1:4],]
sdev <- prin_comp$sdev

### Find the correlataions
var.coord <- t(apply(loadings, 1, var_cor_func, sdev))

## Calculate the Cos2 (square of the coordinates)
var.cos2 <- var.coord^2

## Contribution of variables to the Components ((cos2) * 100) / (total cos2 of the component)
comp.cos2 <- apply(var.cos2, 2, sum)
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
contrib.var <- as.data.frame(var.contrib)

contrib.var[order(contrib.var$PC1, decreasing = T),]

## Find the most contributing variable
head(contrib.var)
load2 <- column_to_rownames(load1, var = "Parameter") %>% t() %>% as.data.frame()

dplyr:: select(load2, matches("HLA|PC"))
t_con <- t(contrib.var) %>% as.data.frame()



dplyr:: select(t_con, matches("HLA|PC")) %>% rowSums()






