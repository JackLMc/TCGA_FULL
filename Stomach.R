# Stomach
pap_clin <- read.csv("./Data/Clinical_data/Paper_Clinical.csv")
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

# source("Clinical.R") # Run to gain the clinical dataframe that's in Output (Clin_STAD)
# load("FPKMs.RData")

# Looking for the files
thousand.folders <- list.dirs(path = "./STAD/Data/FPKM", full.names = T)
filelist1 <- sapply(thousand.folders[-1], function(x){
  list.files(x, pattern = "FPKM.txt.gz$", full.names = T)})
filelist = unlist(filelist1)

# Read in files and combine
lists <- lapply(filelist, read.delim, header = F)
listsDF <- lists
lists <- listsDF

for (i in names(lists)){
  p <- gsub("./STAD/Data/FPKM/", "", i)
  colnames(lists[[i]]) <- c("Gene", p)
}

names(lists) <- gsub("./STAD/Data/FPKM/", "", names(lists))


# Patients I have CIRC scores and microsatellite status for.
# pat_sub <- read.csv("./Data/Important/patient_subtypes.csv")
converter <- read.delim("./STAD/gdc_sample_sheet.2019-03-21.tsv")
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
# write.table("./Output/Patient_list.txt", x = Enrichment_CIRC1$Patient.ID, row.names = F)
# save.image(file = "FPKMs.RData")


########## START ##########
##### Get the Clinical Stuff in ####
Clin_STAD <- pap_clin[, c("TCGA.ID", "Microsatellite.instability") ]
Clin_STAD$Patient.ID <- samptopat(pap_clin$TCGA.ID)
Clin_STAD$Patient.ID <- gsub("-", ".", Clin_STAD$Patient.ID)
Clin_STAD <- Clin_STAD[!is.na(Clin_STAD$Microsatellite.instability),]

Clin_STAD$MSI_STATUS <- ifelse((Clin_STAD$Microsatellite.instability == "MSI-H"), "MSI-H", "MSS")


CIRC_clin <- merge(Enrichment_CIRC1, Clin_STAD, by = "Patient.ID")



# pdf("./STAD/Figures/Clustering/Violin Plot of Mean CIRC.pdf", height = 6, width = 6)
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

##### Analysing variance differences ####
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

#### Clustering ####
# Partitioning clustering
## Remove uneeded stuff
CIRC_for_Cluster <- droplevels(subset(FPKM2, SYMBOL %in% CIRC_genes$CIRC_Genes))
pca1 <- CIRC_for_Cluster %>% 
  gather(contains("TCGA"), key = "Patient.ID", value = "FPKM") %>%
  spread(key = "SYMBOL", value = "FPKM") %>% 
  merge(Clin_STAD, by = "Patient.ID") # Merge with cleaned clinical

my_data <- pca1 %>%
  droplevels() %>%
  #select(matches("HLA|MSI|Patient")) %>%
  na.omit()
head(my_data)
my_data <- my_data[, !'%in%'(colnames(my_data), c("TCGA.ID", "Microsatellite.instability")) ]

## Determine optimal number of clusters for kmeans
library(factoextra)
# pdf("./STAD/Figures/Clustering/Number_kmean_cluster.pdf", width = 6, height = 6)
fviz_nbclust(my_data[, !('%in%'(names(my_data), c("Patient.ID", "MSI_STATUS")))], kmeans, method = "gap_stat")
dev.off()

## Perform Phenograph and kmeans
library(Rphenograph)
set.seed(1)
a <- cbind(my_data, Phenograph_Clusters = factor(Rphenograph(my_data[, !('%in%'(names(my_data), c("Patient.ID", "MSI_STATUS")))])[[2]]$membership), 
           kmeans_Clusters = factor(kmeans(my_data[, !('%in%'(names(my_data), c("Patient.ID", "MSI_STATUS")))], 6)$cluster))

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
# pdf("./STAD/Figures/Clustering/PhenoG_CIRC.pdf", height = 6, width = 6)
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
# pdf("./STAD/Figures/Clustering/CIRC_PhenoGraph_tsNE.pdf")
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

# pdf("./STAD/Figures/Clustering/CIRC_Score_tsNE.pdf")
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
# pdf("./STAD/Figures/Clustering/CIRC_Pheno.pdf", height = 6, width = 6)
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

remove_low_CIRC <- df[tsne_dimensions$Dim1 <= (-5), ]
lessthan <- rownames_to_column(remove_low_CIRC, var = "Patient.ID")
lt <- merge(lessthan, Clin_STAD, by = "Patient.ID")
hiCIRC_pats <- droplevels(subset(lt, MSI_STATUS == "MSS"))$Patient.ID

df1a$Subtype <- ifelse((df1a$Patient.ID %in% hiCIRC_pats), "MSS-hiCIRC", as.character(df1a$MSI_STATUS))

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
dev.off()
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


#### WRITE OUT THE hiCIRC PATIENT SUBTYPES ####
write.csv("./Output/Patient_Subtypes.csv", x = df1a[, c("Patient.ID", "CIRC_Genes", "Subtype")], row.names = F)

load("FPKMs.RData")
pat_sub <- read.csv("./Output/Patient_Subtypes.csv")
library(reshape2)
dcast(pat_sub, Subtype ~.)

##### GENESETS #####
# Bact
rm(book1)
book_list <- list()

FPKM2$SYMBOL[grepl("EP4", FPKM2$SYMBOL)]


book_list[["SAAs"]] <- c("TLR4", "LY96"#,
                         #"TIRAP"#, "MYD88"
                         # ,"IRAK4", "IRAK1",
                         # "TAB1", "TAB2", "MAP3K7",
                         # "MAP2K3", "MAP2K6", "MAP2K4", "MAP2K7",
                         # "MAPK14", "MAPK8"
)
a_list <- list()
a_list[["test"]] <- c("COX2", "PGE2", "EP4")


## ROS
ROS <- read.csv("~/Downloads/GO_term_summary_20190320_151206.csv")
head(ROS)

pat_sub <- df1a[, c("Patient.ID", "CIRC_Genes", "Subtype")]

ROS1 <- droplevels(subset(ROS, Annotated.Term == "cellular response to reactive oxygen species"))

ROS_list <- list()
ROS_list[["ROS"]] <- toupper(levels(ROS1$Symbol))
head(ROS_list)
# try <- FPKM2[FPKM2$SYMBOL %in% book1$SYMBOL, ]

library(GSVA)
Enrichment_book <- gsva(FPKM3, ROS_list)
Enrichment_book1 <- Enrichment_book %>% as.data.frame() %>%
  rownames_to_column(., var = "Geneset") %>%
  gather(contains("TCGA"), key = "Patient.ID", value = "Enrich") %>%
  spread(., key = "Geneset", value = "Enrich")

Enrichment_book1$Patient.ID <- as.factor(Enrichment_book1$Patient.ID)
Enrich <- merge(df1a[, c("Patient.ID", "CIRC_Genes", "Subtype")], Enrichment_book1, by = "Patient.ID")

Enrich1 <- Enrich %>% gather(key = "Parameter", value = "Enrichment", -CIRC_Genes, -Patient.ID, -Subtype)
Enrich1$Parameter <- as.factor(Enrich1$Parameter)

# pdf("./STAD/Figures/Gene_Sets/ROS_response.pdf")
ggplot(Enrich1, aes(x = Subtype, y = Enrichment)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(aes(Subtype, fill = Subtype),
              scale = "width", alpha = 0.8) +
  scale_fill_manual(values = cbcols) +
  labs(x = "MSI Status", y = "Enrichment of GO Pathway 0034614 (cellular response to ROS))") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top") + 
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test")
dev.off()

# Ping
Ping <- read.csv("./Exploratory_Data/Genesets/Ping_Chih_Ho.csv")

Ping_List <- list()
Ping$Name <- as.factor(Ping$Name)
Ping$Gene <- as.factor(Ping$Gene)

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
  ggplot2:: ggsave(filen, plot = temp_plot, device = "pdf", path = "./STAD/Figures/Gene_Sets/Pearson",
                   height = 6, width = 6)
}

## Collated
library(GSVA)
# Pearson correlation across genesets.
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

# Pearson
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
  ggplot2:: ggsave(filen, plot = temp_plot, device = "pdf", path = "./STAD/Figures/",
                   height = 6, width = 6)
}

# Enrichment for Subtype
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
                   path = "./STAD/Figures/",
                   height = 6, width = 6)}


## Individual genes
# how many pats
kable(dcast(pat_sub, Subtype ~., length))

MA <- merge(FPKM, pat_sub[, c("Patient.ID", "Subtype", "CIRC_Genes")], by = "Patient.ID")
genes_of_interest <- c("IL6", "IL1B", "IL23A", "TGFB1",
                       "CCL2", "CCL5", "CXCL10", "CCL20",
                       "CCR6", "TLR4", "TLR2", "CIITA",
                       "RORC", "IL17A", "IL23R",
                       "CDC42", "FCAR", "IGHA1")

GOI <- droplevels(MA[MA$SYMBOL %in% genes_of_interest, ]) %>%
  .[, c("Patient.ID", "SYMBOL", "FPKM", "Subtype", "CIRC_Genes")]
GOI$SYMBOL <- as.factor(GOI$SYMBOL)

### IMPORTANT BIT
GOI1 <- droplevels(subset(GOI, Subtype != "MSI-H"))
### IMPORTANT BIT

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
                   path = "./STAD/Figures/Genes_of_interest/Correlations",
                   height = 6, width = 6)
}




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
                   path = "./STAD/Figures/Genes_of_interest",
                   height = 6, width = 6)}

FPKM2$SYMBOL[grepl("CIAP", FPKM2$SYMBOL)]
GOI <- droplevels(subset(MA, SYMBOL == "TLR4"))
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

Clin_STAD

STADs <- droplevels(subset(pap_clin, Cancer.type == "STAD"))

STADs$Patient.ID <- samptopat(STADs$TCGA.ID)
STADs$Patient.ID <- gsub("-", ".", STADs$Patient.ID)


invest <- merge(STADs, pat_sub, by = "Patient.ID")


head(invest)
invest$Viral <- as.factor(invest$Viral)

levels(invest$Viral)
head(invest)


droplevels(subset(invest, Virus.detection. == "EBV"))


head(invest)
invest$Viral <- ifelse((is.na(invest$Virus.detection.)), "NC", as.character(invest$Virus.detection.))




