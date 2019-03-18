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

load("FPKMs.RData")

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
  colnames(lists[[i]]) <- c("Gene", p)
}

names(lists) <- gsub("./Data/FPKM/", "", names(lists))


# Patients I have CIRC scores and microsatellite status for.
# pat_sub <- read.csv("./Data/Important/patient_subtypes.csv")
converter <- read.delim("./Data/Important/gdc_sample_sheet_FPKM.tsv")
converter$Patient.ID <- gsub("-", ".", converter$Case.ID)
# converter1 <- converter[converter$Patient.ID %in% pat_sub$Patient.ID, ]
converter$File.Name <- gsub(".FPKM.txt.gz", "", converter$File.Name)
# clin <- read.delim("./Data/Clinical/clinical.tsv")
# clin$Patient.ID <- gsub("-", ".", clin$submitter_id)
# clin1 <- clin[clin$Patient.ID %in% pat_sub$Patient.ID, ]

lists1 <- lists[names(lists) %in% converter$File.ID]

multi_join <- function(list_of_loaded_data, join_func, ...){
  require("dplyr")
  output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
  return(output)
}

combined_df <- multi_join(lists1, full_join)

temp_df <- combined_df %>% gather(key = "File.ID", value = "FPKM", -Gene)
converter2 <- droplevels(subset(converter, Sample.Type != "Solid Tissue Normal" &
                                  Sample.Type != "Blood Derived Normal"))



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
FPKMs <- dcast(temp_df1, Gene ~ Patient.ID, sum, value.var = "FPKM")

## Should these be used to talk about library sizes?
# Patient_list <- colnames(FPKMs_cleaned)[!'%in%'(colnames(FPKMs_cleaned), "Gene")]
# write.table(Patient_list, "./Output/Patient_list.txt", row.names = F)

# ##### START EDGER
# library(edgeR)
# # Replenish "x"
# row.names(FPKMs_cleaned) <- NULL
# FPKMs_cleaned <- FPKMs_cleaned %>% column_to_rownames(., var = "Gene")
# x <- DGEList(FPKMs_cleaned)
# 
# # Get groups
# temp <- x$samples %>% 
#   rownames_to_column(., "Patient.ID")
# colnames(pat_sub)[colnames(pat_sub) == "Subtype"] <- "group"
# x$samples <- merge(temp[, c("Patient.ID", "lib.size", "norm.factors")], pat_sub[, c("Patient.ID", "group")], by = "Patient.ID") %>% column_to_rownames(., var = "Patient.ID")
# group <- x$samples$group
# 
# samplenames <- colnames(x)

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

#### Read in data and process ####
FPKM2 <- FPKM1[!is.na(FPKM1$SYMBOL), ]
rownames(FPKM2) <- NULL
FPKM3 <- FPKM2 %>%
  column_to_rownames(., var = "SYMBOL") %>%
  as.matrix()

# CIRC Geneset
CIRC_IG <- read.csv("./Exploratory_Data/Genesets/CIRC.csv")
CIRC_IG$Hugo_Symbol <- as.factor(CIRC_IG$Hugo_Symbol)

#### CIRC ####
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
save.image(file = "FPKMs.RData")
# CIRC_clin <- merge(Enrichment_CIRC1, tcga_pub_clinical, by = "Patient.ID")

# Plot
# pdf("./Figures/Clustering/Violin Plot of Mean CIRC.pdf", height = 6, width = 6)
# ggplot(CIRC_clin, aes(x = MSI_STATUS, y = CIRC_Genes)) +
#   geom_boxplot(alpha = 0.5, width = 0.2) + 
#   geom_violin(aes(MSI_STATUS, fill = MSI_STATUS),
#               scale = "width", alpha = 0.8) +
#   scale_fill_manual(values = cbcols) +
#   labs(x = "MSI Status", y = "CIRC Enrichment Score") +
#   theme_bw() +
#   theme(axis.text = element_text(size = 16)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme(legend.direction = "horizontal", legend.position = "top") + 
#   stat_compare_means(comparisons = list((c("MSI-H", "MSS"))),
#                      label = "p.signif", method = "wilcox.test")
# dev.off()
# 
# 
# ##### Analysing variance differences ####
# # Levene-test (analysis of variance)
# car:: leveneTest(CIRC_clin$CIRC_Genes, group = CIRC_clin$MSI_STATUS, center = "median")
# # var.test(CIRC_Genes ~ MSI_STATUS, data = CIRC_clin)
# fligner.test(CIRC_clin$CIRC_Genes, g = CIRC_clin$MSI_STATUS)
# bartlett.test(CIRC_clin$CIRC_Genes, g = CIRC_clin$MSI_STATUS)
# 
# # Coefficient of variance
# # install_github("benmarwick/cvequality")
# library(cvequality)
# # cv_test <- with(CIRC_clin, asymptotic_test(CIRC_Genes, MSI_STATUS)) # Not significant
# cv_test_MSLRT <- with(CIRC_clin, mslr_test(nr = 1e4, CIRC_Genes, MSI_STATUS))

#### Clustering ####
# Partitioning clustering
## Remove uneeded stuff
CIRC_for_Cluster <- droplevels(subset(FPKM2, SYMBOL %in% CIRC_genes$CIRC_Genes))
pca1 <- CIRC_for_Cluster %>% 
  gather(contains("TCGA"), key = "Patient.ID", value = "FPKM") %>%
  spread(key = "SYMBOL", value = "FPKM") %>% 
  merge(FD1, by = "Patient.ID") # Merge with cleaned clinical
head(pca1)

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
# pdf("./Figures/Clustering/Kmeans_CIRC.pdf", height = 6, width = 6)
# ggbiplot(prin_comp, obs.scale = 1, var.scale = 1, 
#          groups = Kcluster, circle = T, var.axes = F) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme(legend.direction = "horizontal", legend.position = "top") +
#   ylim (-10,10)
# dev.off()

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
set.seed(98)
library(Rtsne)
tsne_out <- Rtsne(as.matrix(df))
tsne_dimensions <- as.data.frame(tsne_out$Y)
colnames(tsne_dimensions) <- c("Dim1", "Dim2")

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

#### Count patients + produce tables (nice) ####
library(knitr)
# install.packages("rmarkdown")
# library("rmarkdown")

kable(dcast(a, MSI_STATUS ~ Phenograph_Clusters, length), caption = "Phenograph clusters")
head(tsne_out)

#### Calculate CIRC Expression for cluster ####
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

head(df)



#### Doesn't seem to work nicely as before ####
less_than_neg10_dim1 <- df[tsne_dimensions$Dim1 <= (-15), ]
lessthan <- rownames_to_column(less_than_neg10_dim1, var = "Patient.ID")
lt <- merge(lessthan, FD1, by = "Patient.ID")
hiCIRC_pats <- droplevels(subset(lt, MSI_STATUS == "MSS"))$Patient.ID


df1a$Subtype <- ifelse((df1a$Patient.ID %in% hiCIRC_pats), "MSS-hiCIRC", as.character(df1a$MSI_STATUS))
head(df1a)



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






thee <- droplevels(subset(df1a, Subtype == "MSS-hiCIRC" & CIRC_Genes <= 0.3))
thee
MSIs <- droplevels(subset(df1a, Subtype == "MSI-H" & CIRC_Genes <= 0.3))
head(MSIs)

top75 <- unname(quantile(MSIs$CIRC_Genes, .25))
top75


these <- merge(Enrichment_CIRC1, pat_sub[, c("Patient.ID", "Subtype")], by = "Patient.ID")
these

range(these$CIRC_Genes)





# find_pheno_comb <- function(Phenotypes){
#   Phenotype1 <- Phenotypes
#   Phenotype2 <- Phenotypes
#   combin <- as.data.frame(cbind(Phenotype1, Phenotype2))
#   combin$Phenotype1 <- as.factor(Phenotype1)
#   combin$Phenotype2 <- as.factor(Phenotype2)
#   combinations <- expand.grid(levels(combin$Phenotype1), levels(combin$Phenotype2))
#   uniqcomb <- combinations[!combinations$Var1==combinations$Var2,]
#   pairs = list()
#   for(i in seq_len(nrow(uniqcomb))){
#     here <- as.vector(unlist(uniqcomb[i,]))
#     pairs[[length(pairs)+1]] <- here
#   }
#   return(pairs)
# }


# Hierarchical clustering
library(gplots)
library(RColorBrewer)
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
mycol <- colorpanel(1000,"blue","white","red")
distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method = "average")
heatmap.2(as.matrix(df),
          col = mycol,
          trace = "none",
          density.info = "none",
          cex.main = 1.5,
          #ColSideColors = col.cell,
          scale = "column",
          margin = c(10,5), lhei = c(2,10),
          hclustfun = hclustAvg,
          RowSideColors = col.cell1)


## WORKS quite well
FD1$MSI_STATUS <- as.factor(FD1$MSI_STATUS)
col.cell1 <- c("#999999","#56B4E9")[FD1$MSI_STATUS]


head(df1)
CIRCs <- droplevels(subset(df1, Subtype != "MSS-hiCIRC"))


range(CIRCs$CIRC_Genes)
MSIs <- droplevels(subset(df1, Subtype == "MSI-H"))

df1$Subtype <- as.factor(df1$Subtype)
levels(df1$Subtype)

#### FPKM
b <- a[, c("Patient.ID", "MSI_STATUS", "Phenograph_Clusters")]
droplevels(subset(b, MSI_STATUS == "MSI-H"))

pat_sub <- read.csv("./Data/Important/patient_subtypes.csv")

these <- merge(b, pat_sub, by = "Patient.ID")
head(these)


View(these)

# Calculate Enrichment of CIRC
library(GSVA)
Enrichment_MSI <- gsva(FPKM3, MSI_genes) 

Enrichment_MSI1 <- Enrichment_MSI %>% as.data.frame() %>%
  rownames_to_column(., var = "Geneset") %>%
  gather(contains("TCGA"), key = "Patient.ID", value = "Enrich") %>%
  spread(., key = "Geneset", value = "Enrich")

Enrichment_MSI1$Patient.ID <- as.factor(Enrichment_MSI1$Patient.ID)