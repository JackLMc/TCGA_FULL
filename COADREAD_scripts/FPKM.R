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

# source("Clinical.R") # Run to gain the clinical dataframe that's in Output (Clin_614)
load("./R_Data/FPKMs.RData")

# Looking for the files
thousand.folders <- list.dirs(path = "./Data/FPKM", full.names = T)
filelist1 <- sapply(thousand.folders[-1], function(x){
  list.files(x, pattern = "FPKM.txt.gz$", full.names = T)})
filelist = unlist(filelist1)

# Read in files and combine
lists <- lapply(filelist, read.delim, header = F)
listsDF <- lists
lists <- listsDF

for (i in names(lists)){
  p <- gsub("./Data/FPKM/", "", i)
  colnames(lists[[i]]) <- c("Gene", p)}

names(lists) <- gsub("./Data/FPKM/", "", names(lists))


# GDC sample sheet
converter <- read.delim("./Data/Important/gdc_sample_sheet_FPKM.tsv")
converter$Patient.ID <- gsub("-", ".", converter$Case.ID)
# converter1 <- converter[converter$Patient.ID %in% pat_sub$Patient.ID, ]
converter$File.Name <- gsub(".FPKM.txt.gz", "", converter$File.Name)
# clin <- read.delim("./Data/Clinical/clinical.tsv")
# clin$Patient.ID <- gsub("-", ".", clin$submitter_id)
# clin1 <- clin[clin$Patient.ID %in% pat_sub$Patient.ID, ]

lists1 <- lists[names(lists) %in% converter$File.ID]
converter <- converter[converter$File.ID %in% names(lists), ]

multi_join <- function(list_of_loaded_data, join_func, ...){
  require("dplyr")
  output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
  return(output)
}

combined_df <- multi_join(lists1, full_join, by = "Gene")

temp_df <- combined_df %>% gather(key = "File.ID", value = "FPKM", -Gene)
converter2 <- droplevels(subset(converter,  Sample.Type == "Primary Tumor"))

# this <- data.frame(stringsAsFactors = F)
# converter2$Patient.ID <- as.factor(converter2$Patient.ID)
# c <- 1
# for(i in levels(converter2$Patient.ID)){
#   work <- droplevels(subset(converter2, Patient.ID == i))
#   this[c, "Patient.ID"] <- i
#   this[c, "Sample_lev"] <- nlevels(work$Sample.ID)
#   this[c, "File_ID_lev"] <- nlevels(work$File.ID)
#   this[c, "File_Name_lev"] <- nlevels(work$File.Name)
#   c <- c + 1
#   }


temp_df1 <- merge(converter2[, c("Patient.ID", "File.ID")], temp_df, by = "File.ID")

library(reshape2)
head(temp_df1)
FPKMs <- dcast(temp_df1, Gene ~ Patient.ID, sum, value.var = "FPKM")

## Should these be used to talk about library sizes?
# Patient_list <- colnames(FPKMs)[!'%in%'(colnames(FPKMs), "Gene")]
# write.table(Patient_list, "./Output/Patient_list.txt", row.names = F)

# Gaining second dataframe (Symbols)
# biocLite("Homo.sapiens", dependencies = T)
library(Homo.sapiens)
FPKMs$Gene <- gsub("\\..*", "", FPKMs$Gene) #Removes version from Ensembl gene ID
geneid <- FPKMs$Gene

genes <- select(Homo.sapiens, keys = geneid, columns = c("SYMBOL", "TXCHROM", "ENTREZID"), 
                keytype = "ENSEMBL")
genes <- genes[!duplicated(genes$SYMBOL),]

colnames(FPKMs)[colnames(FPKMs) %in% "Gene"] <- "ENSEMBL"
FPKM <- merge(genes[, c("ENSEMBL", "SYMBOL")], FPKMs, by = "ENSEMBL") %>% 
  gather(contains("TCGA"), key = "Patient.ID", value = "FPKM")
library(reshape2)
FPKM1 <- dcast(FPKM, SYMBOL ~ Patient.ID, sum, value.var = "FPKM")
FPKM2 <- FPKM1[!is.na(FPKM1$SYMBOL), ]
rownames(FPKM2) <- NULL
FPKM3 <- FPKM2 %>%
  column_to_rownames(., var = "SYMBOL") %>%
  as.matrix()

# CIRC Score calculation ----
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

Enrichment_CIRC1$Patient.ID <- as.factor(Enrichment_CIRC1$Patient.ID)
# write.table("./Output/Patient_list.txt", x = Enrichment_CIRC1$Patient.ID, row.names = F)
# save.image(file = "./R_Data/FPKMs.RData")


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
set.seed(1)
library(Rtsne)
tsne_out <- Rtsne(as.matrix(df))
tsne_dimensions <- as.data.frame(tsne_out$Y)
colnames(tsne_dimensions) <- c("Dim1", "Dim2")

head(tsne_dimensions)

## tSNE plot - looks the same as PCA just on a different axis
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# pdf("./Figures/Clustering/CIRC_PhenoGraph_tsNE.pdf")
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

for_CIRC <- Enrichment_CIRC1[Enrichment_CIRC1$Patient.ID %in% rownames(df), ]
CIRC_score <- for_CIRC[, "CIRC_Genes"]

pdf("./Figures/Clustering/CIRC_Score_tsNE.pdf")
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


# Count patients
library(knitr)
library(reshape2)
# install.packages("rmarkdown")
# library("rmarkdown")

kable(dcast(a, MSI_STATUS ~ Phenograph_Clusters, length), caption = "Phenograph clusters")
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
remove_low_CIRC <- df[tsne_dimensions$Dim1 >= (-25) & tsne_dimensions$Dim1 <= -5, ]
lessthan <- rownames_to_column(remove_low_CIRC, var = "Patient.ID")
lt <- merge(lessthan, Clin_614, by = "Patient.ID")
hiCIRC_pats <- droplevels(subset(lt, MSI_STATUS == "MSS"))$Patient.ID
df1a$Subtype <- ifelse((df1a$Patient.ID %in% hiCIRC_pats), "MSS-hiCIRC", as.character(df1a$MSI_STATUS))

# Remove the MSS-hiCIRC patients who are POLE mutants
POLE <- read.csv("Output/POLE_mutants.csv")$.
df1a$Subtype <- ifelse((df1a$Patient.ID %in% POLE), as.character(df1a$MSI_STATUS), as.character(df1a$Subtype))

ggplot(df1a, aes(x = Subtype, y = CIRC_Genes)) +
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

# START ----
load("./R_Data/FPKMs.RData")
pat_sub <- read.csv("./Output/Patient_Subtypes.csv")
library(reshape2)
dcast(pat_sub, Subtype ~., length)

# GENESET Interrogation ----
# BiocManager::install("GSVA")
library(GSVA)

# Ping
Ping <- read.csv("./Exploratory_Data/Genesets/Ping_Chih_Ho.csv")
Ping <- factorthese(Ping, c("Name", "Gene"))

Ping_List <- list()
c <- 1
for(i in levels(Ping$Name)){
  print(i)
  work <- droplevels(subset(Ping, Name == i))
  Genes <- levels(work$Gene)
  Ping_List[[i]] <- Genes
  c <- c + 1
  }

Enrichment_Ping <- gsva(FPKM3, Ping_List)
Enrichment_Ping1 <- Enrichment_Ping %>% as.data.frame() %>%
  rownames_to_column(., var = "Geneset") %>%
  gather(contains("TCGA"), key = "Patient.ID", value = "Enrich") %>%
  spread(., key = "Geneset", value = "Enrich")

Enrichment_Ping1$Patient.ID <- as.factor(Enrichment_Ping1$Patient.ID)
Enrich <- merge(pat_sub, Enrichment_Ping1, by = "Patient.ID")

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
## Pearson correlation across genesets.
CTGenesets <- read.csv("./Exploratory_Data/Genesets/Cell_Type_Geneset.csv")
SigGenesets <- read.csv("./Exploratory_Data/Genesets/Signature_Geneset.csv")
Genesets <- rbind(CTGenesets, SigGenesets)

dd <- deduplicate(Genesets)

geneset_list <- list()
for(i in levels(Genesets$Parameter)){
  print(i)
  work <- droplevels(subset(Genesets, Parameter == i))
  genes <- levels(work$Hugo_Symbol)
  geneset_list[[i]] <- genes
  }

Enrichments <- gsva(FPKM3, geneset_list) %>% as.data.frame() %>%
  rownames_to_column(., "Geneset") %>% gather(contains("TCGA"), key = "Patient.ID", value = "Enrichment") %>%
  merge(., pat_sub, by = "Patient.ID")

Enrichments$Geneset <- as.factor(Enrichments$Geneset)

## Pearson
for(i in levels(Enrichments$Geneset)){
  print(i)
  work <- droplevels(subset(Enrichments, Geneset == i))
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

Enrichment_book <- gsva(FPKM3, ROS_list)
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
Enrichment_book <- gsva(FPKM3, FAM_list)
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


# Phagocytosis
Phago  <- read.csv("./Exploratory_Data/Genesets/GO_term_summary_20190603_065405.csv")

Phago_list <- list()
c <- 1
for(i in levels(Phago$Annotated.Term)){
  print(i)
  work <- droplevels(subset(Phago, Annotated.Term == i))
  Genes <- toupper(levels(work$Symbol))
  Phago_list[[i]] <- Genes
  c <- c + 1
}

library(GSVA)
Enrichment_book <- gsva(FPKM3, Phago_list)
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
                   path = "./Figures/Gene_Sets/Enrichment/Phago",
                   height = 6, width = 6)}



## Bespoke
FPKM2$SYMBOL[grepl("S100A", FPKM2$SYMBOL)]
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
Enrichment_book <- gsva(FPKM3, a_list)
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
MA <- merge(FPKM, pat_sub[, c("Patient.ID", "Subtype", "CIRC_Genes")], by = "Patient.ID")
genes_of_interest <- c("IL6", "IL1B", "IL23A", "TGFB1",
                       "CCL2", "CCL5", "CXCL10", "CCL20",
                       "CCR6", "TLR4", "TLR2", "CIITA",
                       "RORC", "IL17A", "IL23R",
                       "CDC42", "FCAR", "IGHA1", "EHMT2", "OX40",
                       "OX40L", "CD1D")

GOI <- droplevels(MA[MA$SYMBOL %in% genes_of_interest, ]) %>%
  .[, c("Patient.ID", "SYMBOL", "FPKM", "Subtype", "CIRC_Genes")]
GOI$SYMBOL <- as.factor(GOI$SYMBOL)

## Remove for Pearson correlation
GOI1 <- droplevels(subset(GOI, Subtype != "MSI-H"))
for(i in levels(GOI1$SYMBOL)){
  print(i)
  work <- droplevels(subset(GOI1, SYMBOL == i))
  work$Logged <- log2(work$FPKM + 1)
  temp_plot <- ggplot(work, aes(y = Logged, x = CIRC_Genes))+
    geom_point(alpha = 0.8, size = 4, colour = "slategray") +
    labs(x = "CIRC_Genes Enrichment", y = paste0("Log2(", i, " + 1)")) +
    theme_bw() +
    # geom_text(aes(x = -0.3, y = .75, label = lm_eqn(lm(CIRC_Genes ~ Enrichment, work))), parse = T) +
    # scale_color_manual(values = cbcols) +
    geom_smooth(method = "lm", se = F) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.position = "none") +
    stat_cor()
  filen <- paste0(i, ".pdf")
  ggplot2:: ggsave(filen, plot = temp_plot, device = "pdf",
                   path = "./Figures/Genes_of_interest/Correlations",
                   height = 6, width = 6)
}

## Compare across Subtypes
for(i in 1:length(genes_of_interest)){
  gene <- genes_of_interest[i]
  print(gene)
  GOI <- droplevels(subset(MA, SYMBOL == gene))
  GOI$Rank <- rank(GOI$FPKM)
  temp_plot <- ggplot(GOI, aes(x = Subtype, y = Rank)) +
    geom_boxplot(alpha = 0.5, width = 0.2) + 
    geom_violin(aes(Subtype, fill = Subtype),
                scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) +
    labs(x = "MSI Status", y = paste("Rank transformed", gene, "FPKM", sep = " ")) +
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
FPKM2$SYMBOL[grepl("B2M", FPKM2$SYMBOL)] # Check whether your gene exists in the dataset
GOI <- droplevels(subset(MA, SYMBOL == "ADGRB1")) 

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



# THIS PART OF THE SCRIPT IS TO INVESTIGATE WHATEVER GARY WANTS ----
work <- droplevels(subset(GOI1, SYMBOL == "EHMT2" | SYMBOL == "CIITA"))
work$Logged <- log2(work$FPKM + 1)

work1 <- spread(work[, c("Patient.ID", "Subtype", "SYMBOL", "Logged")], key = "SYMBOL", value = "Logged")

temp_plot <- ggplot(work1, aes(y = EHMT2, x = CIITA))+
  geom_point(alpha = 0.8, size = 4, colour = "slategray") +
  labs(x = "Log2(CIITA FPKM + 1)", y = "Log2(EHMT2 FPKM + 1)") +
  theme_bw() +
  # geom_text(aes(x = -0.3, y = .75, label = lm_eqn(lm(CIRC_Genes ~ Enrichment, work))), parse = T) +
  # scale_color_manual(values = cbcols) +
  geom_smooth(method = "lm", se = F) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  stat_cor()


filen <- paste0(i, ".pdf")
ggplot2:: ggsave("CIITA_EHMT2_correl.pdf", plot = temp_plot, device = "pdf",
                 path = "./Figures/Genes_of_interest/Correlations",
                 height = 6, width = 6)


