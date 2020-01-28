# A script to find MSS-hiCIRC patients
## Packages
library(UsefulFunctions)
library(tidyverse)

## Load required RData objects
load("./R_Data/FPKMs.RData")



## CIRC Score calculation (enrichment of the CIRC gene FPKM) ----
CIRC_IG <- read.csv("./Exploratory_Data/Genesets/CIRC.csv")
CIRC_IG$SYMBOL <- as.factor(CIRC_IG$SYMBOL)

# Label
CIRC_genes <- droplevels(subset(CIRC_IG, CIRC == T)) %>% 
  takegenelevels(.) %>% list()
names(CIRC_genes) <- "CIRC_Genes"

# Calculate Enrichment of CIRC
library(GSVA)
Enrichment_CIRC <- gsva(FPKM3, CIRC_genes) 

Enrichment_CIRC1 <- Enrichment_CIRC %>% as.data.frame() %>%
  rownames_to_column(., var = "Geneset") %>%
  gather(contains("TCGA"), key = "Patient.ID", value = "Enrich") %>%
  spread(., key = "Geneset", value = "Enrich")


# START CIRC Enrichment --------
# Read Clinical Stuff in ----
Clin_614 <- read.csv("./Output/Clinical_Data_614.csv")
CIRC_clin <- merge(Enrichment_CIRC1, Clin_614, by = "Patient.ID")

# pdf("./Figures/Clustering/Violin Plot of Mean CIRC.pdf", height = 6, width = 6)
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
# dev.off()

kable(dcast(CIRC_clin, MSI_STATUS ~ ., length))


# Analysing variance differences
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

library(knitr)
library(kableExtra)
library(magrittr)
rbind(Test, asym)

# kable_out <- knitr::kable(rbind(Test, asym), "html") %>% kableExtra::kable_styling(bootstrap_options = c("striped", "hover"))%>%
#   kable_styling()
# readr::write_file(kable_out, "kable_out.html")

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

## Determine optimal number of clusters for kmeans
library(factoextra)
pdf("./Figures/Clustering/Number_kmean_cluster.pdf", width = 6, height = 6)
fviz_nbclust(my_data[, !('%in%'(names(my_data), c("Patient.ID", "MSI_STATUS")))], kmeans, method = "gap_stat")
dev.off()

## Perform Phenograph and kmeans
# BiocManager::install("devtools")
devtools:: install_github("https://github.com/JinmiaoChenLab/Rphenograph")
library(Rphenograph)
set.seed(1)
a <- cbind(my_data, Phenograph_Clusters = factor(Rphenograph(my_data[, !('%in%'(names(my_data), c("Patient.ID", "MSI_STATUS")))])[[2]]$membership), 
           kmeans_Clusters = factor(kmeans(my_data[, !('%in%'(names(my_data), c("Patient.ID", "MSI_STATUS")))], 10)$cluster))

# Calculate PCs
pca1a <- data.frame(a[, names(a) != "MSI_STATUS" & 
                        names(a) != "Phenograph_Clusters" & 
                        names(a) != "kmeans_Clusters"])

# Make the Patient ID the row names
df <- data.frame(pca1a[, names(pca1a) != "Patient.ID"], row.names = pca1a[, names(pca1a) == "Patient.ID"])

## Do PCA
prin_comp <- prcomp(df, scale. = T)

## Plot them 
### kmeans
# Kcluster <- a[, "kmeans_Clusters"]
library(ggbiplot)
### Phenograph
Pcluster <- a[, "Phenograph_Clusters"]
pdf("./Figures/Clustering/PhenoG_CIRC.pdf", height = 6, width = 6)
ggbiplot(prin_comp, obs.scale = 1, var.scale = 1, 
         groups = Pcluster, circle = T, var.axes = F) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top") 
dev.off()


# RTsne
# set.seed(1)
# library(Rtsne)
# tsne_out <- Rtsne(as.matrix(df))
# tsne_dimensions <- as.data.frame(tsne_out$Y)
# colnames(tsne_dimensions) <- c("Dim1", "Dim2")
# 
# head(tsne_dimensions)
# 
# ## tSNE plot - looks the same as PCA just on a different axis
# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# # pdf("./Figures/Clustering/CIRC_PhenoGraph_tsNE.pdf")
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
# 
# for_CIRC <- Enrichment_CIRC1[Enrichment_CIRC1$Patient.ID %in% rownames(df), ]
# CIRC_score <- for_CIRC[, "CIRC_Genes"]
# 
# pdf("./Figures/Clustering/CIRC_Score_tsNE.pdf")
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
#                      space = "Lab", na.value = "grey50", guide = "colourbar",
#                      aesthetics = "colour")
# dev.off()


## UMAP
install.packages("umap")
library(umap)
set.seed(1)

umap_out <- umap(as.matrix(df))
str(umap_out)

umap_dimensions <- as.data.frame(umap_out$layout)
colnames(umap_dimensions) <- c("Dim1", "Dim2")

## umap plot - looks the same as PCA just on a different axis
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# pdf("./Figures/Clustering/CIRC_PhenoGraph_umap.pdf")
ggplot(umap_dimensions, aes(x = Dim1, y = Dim2, colour = Pcluster)) +
  geom_point(size = 4, alpha = 0.8, pch = 20) +
  scale_colour_manual(values = c("1" = "#009E73", "2" = "#56B4E9", "3" = "#E69F00",
                                 "4" = "#CC79A7", "5" = "#0072B2", "6" = "#999999",
                                 "7" = "#F0E442", "8" = "#D55E00")) +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top")
dev.off()

for_CIRC <- Enrichment_CIRC1[Enrichment_CIRC1$Patient.ID %in% rownames(df), ]
CIRC_score <- for_CIRC[, "CIRC_Genes"]

pdf("./Figures/Clustering/CIRC_Score_umap.pdf")
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





# Count patients
library(knitr)
library(reshape2)
# install.packages("rmarkdown")
# library("rmarkdown")

dcast(a, MSI_STATUS ~ Phenograph_Clusters, length)
head(tsne_out)

# Calculate CIRC Expression for clusters ----
j <- a[, c("Patient.ID", "kmeans_Clusters", "Phenograph_Clusters", "MSI_STATUS")]
df1a <- merge(Enrichment_CIRC1, j, by = "Patient.ID")
df1 <- droplevels(subset(df1a, MSI_STATUS == "MSS"))

# Phenograph
pdf("./Figures/Clustering/CIRC_Pheno.pdf", height = 6, width = 6)
ggplot(df1, aes(x = Phenograph_Clusters, y = CIRC_Genes)) +
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
  stat_compare_means(comparisons = list(c("1", "2"), c("1", "3"), c("1", "4"), c("1", "5"), c("1", "6"), c("1", "7"), c("1", "8"),
                                        c("2", "3"), c("2", "4"), c("2", "5"), c("2", "6"), c("2", "7"), c("2", "8"),
                                        c("3", "4"), c("3", "5"), c("3", "6"), c("3", "7"), c("3", "8")
                                        # c("4", "5"), c("4", "6"), c("4", "7"), c("4", "8"),
                                        # c("5", "6"), c("5", "7"), c("5", "8"),
                                        # c("6", "7"), c("6", "8"),
                                        # c("7", "8")
  ),
  label = "p.signif", method = "wilcox.test")
dev.off()

# Determine CIRC patients based on CIRC score
## UMAP
hiCIRC_all <-  df[umap_dimensions$Dim2 < 0, ]
hiCIRC_all1 <-  rownames_to_column(remove_low_CIRC, var = "Patient.ID") %>%
  merge(., Clin_614, by = "Patient.ID") %>% subset(., MSI_STATUS == "MSS") %>% droplevels()
df1a$Subtype_UMAP <- ifelse((df1a$Patient.ID %in% hiCIRC_all1$Patient.ID), "MSS-hiCIRC", as.character(df1a$MSI_STATUS))

# ## tSNE
# remove_low_CIRC <- df[tsne_dimensions$Dim1 >= (-25) & tsne_dimensions$Dim1 <= -5, ]
# lessthan <- rownames_to_column(remove_low_CIRC, var = "Patient.ID")
# lt <- merge(lessthan, Clin_614, by = "Patient.ID")
# hiCIRC_pats <- droplevels(subset(lt, MSI_STATUS == "MSS"))$Patient.ID
# df1a$Subtype_tSNE <- ifelse((df1a$Patient.ID %in% hiCIRC_pats), "MSS-hiCIRC", as.character(df1a$MSI_STATUS))

# Remove the MSS-hiCIRC patients who are POLE mutants
POLE <- read.csv("Output/POLE_mutants.csv")$.
df1a$Subtype_UMAP <- ifelse((df1a$Patient.ID %in% POLE), as.character(df1a$MSI_STATUS), as.character(df1a$Subtype_UMAP))
df1a$Subtype_tSNE <- ifelse((df1a$Patient.ID %in% POLE), as.character(df1a$MSI_STATUS), as.character(df1a$Subtype_tSNE))

ggplot(df1a, aes(x = Subtype_UMAP, y = CIRC_Genes)) +
  geom_boxplot(alpha = 0.5, width = 0.2) +
  geom_violin(aes(Subtype_UMAP, fill = Subtype_UMAP),
              scale = "width", alpha = 0.8) +
  scale_fill_manual(values = cbcols) +
  labs(x = "MSI Status", y = "CIRC Enrichment Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top") +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test")


head(df1a)
View(df1a)

hiCIRC_old <- read.csv("./Output/Patient_Subtypes.csv")

dim(hiCIRC_old)
dim(df1a)

all.equal(df1a$Subtype_UMAP, df1a$Subtype_tSNE)

all.equal(df1a$Subtype_UMAP, hiCIRC_old$Subtype)

df1a$Subtype_UMAP <- as.factor(df1a$Subtype_UMAP)




head(df1a)
View(df1a)
tSNE_hiCIRC <- droplevels(subset(df1a, Subtype_UMAP == "MSS-hiCIRC"))
dim(tSNE_hiCIRC)

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
