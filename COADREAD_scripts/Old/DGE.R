# A script to look at differential gene expression of COADREAD TCGA project
library(UsefulFunctions)
library(tidyverse)
library(ggpubr)

my_comparisons <- list(c("MSS-hiCIRC", "MSI-H"),
                       c("MSS-hiCIRC", "MSS"),
                       c("MSI-H", "MSS"))
cbcols <- c("MSS-hiCIRC" = "#999999",
            "MSI-H" = "#56B4E9",
            "MSS" = "#009E73")

load("./R_Data/Counts_clean.RData")
# rm(list = setdiff(ls(), c("cqn_Counts"))) # Clean environment

##### START EDGER
library(edgeR)

# Replenish "x"
Counts_clean <- cqn_Counts$counts
x <- DGEList(Counts_clean)

# Get groups
temp <- x$samples %>% 
  rownames_to_column(., "Patient.ID")
pat_sub <- read.csv("./Output/Patient_Subtypes_09_03.csv")

colnames(pat_sub)[colnames(pat_sub) == "Subtype"] <- "group"
x$samples <- merge(temp[, c("Patient.ID", "lib.size", "norm.factors")], pat_sub[, c("Patient.ID", "group")], by = "Patient.ID") %>%
  column_to_rownames(., var = "Patient.ID") # Pat_sub has less patients in than the raw count dataframe
group <- x$samples$group

x$counts <- x$counts[, colnames(x$counts) %in% pat_sub$Patient.ID] # Only include those with the subtypes
samplenames <- colnames(x)

# Calculate counts per million, log counts per million
cpms <- cpm(x)
lcpm <- cpm(x, log = T)

# Remove lowly expressed transcripts
dim(x)
table(rowSums(x$counts==0)==10) # Show me the amount of transcripts that are zero for 10 samples
CPM_scaling <- min(x$samples$lib.size)/1000000
cpm_scale <- 6.5/CPM_scaling

keep.exprs <- rowSums(cpms>cpm_scale)>=1 # Genes must count of 6.5 in lowest library
x <- x[keep.exprs,, keep.lib.sizes = F]

# Normalising gene expression
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors
x <- estimateCommonDisp(x)
x <- estimateTagwiseDisp(x)


# ### Double check...
# lcpm <- cpm(x, log = T)
# col.cell <- c("#999999","#56B4E9","#E69F00")[x$samples$group]
# plotMDS(cpms, pch = 16, cex = 2, col = cbcols)
# legend("top",
#        fill = c("#999999", "#56B4E9",
#                 "#E69F00"),
#        legend = levels(x$samples$group))
# ## Doesn't really separate them based on top 500 varying genes...
# 
# # Look at variable genes
# var_genes <- apply(lcpm, 1, var)
# select_var <- names(sort(var_genes, decreasing = T))[1:150]
# 
# # Subset logcounts matrix
# highly_variable_lcpm <- lcpm[select_var,]
# # dim(highly_variable_lcpm)
# # head(highly_variable_lcpm)
# 
# ## Get some nicer colours
# library(gplots)
# library(RColorBrewer)
# mypalette <- brewer.pal(11,"RdYlBu")
# morecols <- colorRampPalette(mypalette)
# mycol <- colorpanel(1000,"blue","white","red")
# 
# ## Get the SYMBOLS instead of ENTREZ ID
# HvL <- rownames_to_column(as.data.frame(highly_variable_lcpm), var = "ENSEMBL")
# distCor <- function(x) as.dist(1-cor(t(x)))
# hclustAvg <- function(x) hclust(x, method = "average")
# col.cell1 <- c("#999999","#56B4E9","#E69F00")[x$samples$group]
# 
# HvL2 <- highly_variable_lcpm_sym[, names(highly_variable_lcpm_sym) != "SYMBOL"] %>%
#   column_to_rownames(., "ENTREZID")
# 
# heatmap.2(as.matrix(highly_variable_lcpm),
#           col = mycol,
#           trace = "none",
#           density.info = "none",
#           main = "B. Heatmap of the Top 150 most Variable Genes",
#           cex.main = 1.5,
#           #ColSideColors = col.cell,
#           scale = "row",
#           margin = c(10,5), lhei = c(2,10),
#           labCol = colnames(highly_variable_lcpm),
#           labRow = rownames(highly_variable_lcpm),
#           hclustfun = hclustAvg,
#           ColSideColors = col.cell1)




# Differential gene expression
group <- gsub("-", "_", group)
x$samples$group <- gsub("-", "_", x$samples$group)
design <- model.matrix(~0 + group)
colnames(design) <- gsub("group", "", colnames(design))

contr.matrix <- makeContrasts(
  MSI_HvsMSS_hiCIRC = MSI_H - MSS_hiCIRC, 
  MSI_HvsMSS = MSI_H - MSS, 
  MSS_hiCIRCvsMSS = MSS_hiCIRC - MSS, 
  levels = colnames(design))

## Removing heteroscedascity
# biocLite("limma", dependencies = T)
library(limma)
# pdf("../Figures/Proof/Mean-variance-tred.pdf")
# par(mfrow = c(1,2))
v <- voom(x, design, plot = F)

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit, robust = T)
# plotSA(efit, main = "Final model: Mean−variance trend")

## Number of differentially expressed genes
summary(decideTests(efit))

### More stringently selected DE genes
tfit <- treat(vfit, lfc = 1)
summary(decideTests(tfit))

dt <- decideTests(tfit)
# summary(dt)
# de.common <- which(dt[,2]!=0 & dt[,3]!=0)
# length(de.common)


# Finding the genes
## MSS-hiCIRC versus MSS
# MSS_hiCIRC_MSS <- which(dt[, "MSS_hiCIRCvsMSS"] != 0)
# ENSEMBL_MSSs <- efit$genes$SYMBOL[MSS_hiCIRC_MSS]

MSS_hiCIRC_MSI <- which(dt[, "MSI_HvsMSS_hiCIRC"] != 0)

str(efit)
str(tfit)

efit$t

Genes_MSI_hiCIRC <- names(MSS_hiCIRC_MSI)

Genes_MSI_hiCIRC[grepl("ROR", Genes_MSI_hiCIRC)]
Genes_MSI_hiCIRC[grepl("IL17", Genes_MSI_hiCIRC)]
Genes_MSI_hiCIRC[grepl("CCL", Genes_MSI_hiCIRC)]
Genes_MSI_hiCIRC[grepl("IL", Genes_MSI_hiCIRC)]
Genes_MSI_hiCIRC[grepl("CCR", Genes_MSI_hiCIRC)]
Genes_MSI_hiCIRC[grepl("TLR", Genes_MSI_hiCIRC)]
Genes_MSI_hiCIRC[grepl("MUC", Genes_MSI_hiCIRC)]

# i <- which(v$genes$ENSEMBL %in% Symbol_MSI)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale = "row",
#           labRow = v$genes$SYMBOL[i], labCol = colnames(v), #can also do labCol = groups
#           col = mycol, trace = "none", density.info = "none",
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column",
#           main = "Shared DE genes (Naive v Effectors)",
#           hclustfun = hclustAvg)
# dev.off()

# Looking at the data
MSI_H.vs.MSS_hiCIRC <- topTreat(tfit, coef = 1, n = Inf)
MSI_H.vs.MSS <- topTreat(efit, coef = 2, n = Inf)
MSS_hiCIRC.vs.MSS <- topTreat(efit, coef = 3, n = Inf)

# Volcano plots...
myData <- as.data.frame(MSI_H.vs.MSS_hiCIRC)
myData$padjThresh <- as.factor(myData$adj.P.Val < 0.05)

take_a_look <- myData[myData$adj.P.Val <= 0.01, ]
View(take_a_look)

head(myData)
myData$Labels <- rownames(myData)
labelled_genes <- c("RORC", "TLR4", "CCR6")
myData$Labels[!rownames(myData) %in% labelled_genes] <- ""
library(ggrepel)

# pdf(file = "../Figures/Figure_2/Figure2B.pdf", width = 6, height = 6)
ggplot(data = myData, aes(x = logFC, y = -log10(P.Value))) + 
  geom_point(alpha = 0.2, size = 1) +
  theme(legend.position = "none", text = element_text(size = 10)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ## Colours for significance
  geom_point(data = subset(myData, padjThresh==T & logFC<(-0.5)), aes(logFC, -log10(P.Value)), alpha = 0.6, size = 0.6, colour = "blue4") +
  geom_point(data = subset(myData, padjThresh==T & logFC>0.5), aes(logFC, -log10(P.Value)), alpha = 0.6, size = 0.6, colour = "orangered3") +
  ## Label these
  geom_label_repel(mapping = aes(label = Labels), box.padding = 0.35, point.padding = 0.2) +
  # geom_text_repel(data=subset(myData,padjThresh==TRUE & logFC>0.5), aes(logFC,-log10(P.Value),label=SYMBOL), nudge_x = 0.05, colour="black",force=0.1,size=1.25,segment.size=0.1,segment.alpha = 0.5) +
  labs(x = expression(Log[2]*" fold change"), y = expression(-Log[10]*" p-value"))
dev.off()





MSI_H.vs.MSS[grepl("RORC", rownames(MSI_H.vs.MSS)), ]
MSI_H.vs.MSS_hiCIRC[grepl("RORC", rownames(MSI_H.vs.MSS_hiCIRC)), ]
MSS_hiCIRC.vs.MSS[grepl("RORC", rownames(MSS_hiCIRC.vs.MSS)), ]

MSI_H.vs.MSS[grepl("IL17A", rownames(MSI_H.vs.MSS)), ]
MSI_H.vs.MSS_hiCIRC[grepl("IL17A", rownames(MSI_H.vs.MSS_hiCIRC)), ]
MSS_hiCIRC.vs.MSS[grepl("IL17A", rownames(MSS_hiCIRC.vs.MSS)), ]

MSI_H.vs.MSS[grepl("IL23R", rownames(MSI_H.vs.MSS)), ]
MSI_H.vs.MSS_hiCIRC[grepl("IL23R", rownames(MSI_H.vs.MSS_hiCIRC)), ]
MSS_hiCIRC.vs.MSS[grepl("IL23R", rownames(MSS_hiCIRC.vs.MSS)), ]

MSI_H.vs.MSS[grepl("CCL20", rownames(MSI_H.vs.MSS)), ]
MSI_H.vs.MSS_hiCIRC[grepl("CCL20", rownames(MSI_H.vs.MSS_hiCIRC)), ]
MSS_hiCIRC.vs.MSS[grepl("CCL20", rownames(MSS_hiCIRC.vs.MSS)), ]

MSI_H.vs.MSS[grepl("CCR6", rownames(MSI_H.vs.MSS)), ]
MSI_H.vs.MSS_hiCIRC[grepl("CCR6", rownames(MSI_H.vs.MSS_hiCIRC)), ]
MSS_hiCIRC.vs.MSS[grepl("CCR6", rownames(MSS_hiCIRC.vs.MSS)), ]


MSI_H.vs.MSS[grepl("CIITA", rownames(MSI_H.vs.MSS)), ]
MSI_H.vs.MSS_hiCIRC[grepl("CIITA", rownames(MSI_H.vs.MSS_hiCIRC)), ]
MSS_hiCIRC.vs.MSS[grepl("CIITA", rownames(MSS_hiCIRC.vs.MSS)), ]

MSI_H.vs.MSS[grepl("HLA-DR", rownames(MSI_H.vs.MSS)), ]
MSI_H.vs.MSS_hiCIRC[grepl("HLA-DR", rownames(MSI_H.vs.MSS_hiCIRC)), ]
MSS_hiCIRC.vs.MSS[grepl("HLA-DR", rownames(MSS_hiCIRC.vs.MSS)), ]

# save.image(file = "./R_Data/Counts.RData")

see_pats <- merge(temp, pat_sub, by = "Patient.ID")
head(see_pats)

nrow(droplevels(subset(see_pats, group.y == "MSS-hiCIRC")))
nrow(droplevels(subset(see_pats, group.y == "MSS")))
nrow(droplevels(subset(see_pats, group.y == "MSI-H")))


##################################################################

# Geneset enrichment analysis
## Open connections
tryCatch(expr = { library("limma")}, 
         error = function(e) { source("https://bioconductor.org/biocLite.R")
           biocLite("limma")}, 
         finally = library("limma"))

tryCatch(expr = { library("GSA")}, 
         error = function(e) { source("https://bioconductor.org/biocLite.R")
           biocLite("GSA")}, 
         finally = library("GSA"))

tryCatch(expr = { library("RCurl")}, 
         error = function(e) { 
           install.packages("RCurl")}, 
         finally = library("RCurl"))


library(qusage)

## Read in the Genesets
All_gmt <- read.gmt("./Data/Genesets/GSEA/Symbol/msigdb.v7.0.symbols.gmt")
KEGG_gmt <- read.gmt("./Data/Genesets/GSEA/Symbol/c2.cp.kegg.v7.0.symbols.gmt")
GO_gmt <- read.gmt("./Data/Genesets/GSEA/Symbol/c5.all.v7.0.symbols.gmt")
immuno_gmt <- read.gmt("./Data/Genesets/GSEA/Symbol/c7.all.v7.0.symbols.gmt")
hallmark_gmt <- read.gmt("./Data/Genesets/GSEA/Symbol/h.all.v7.0.symbols.gmt")

## Filter genesets that appear in only KEGG and GO databases (6103 genesets)
# filter_gmt <- All_gmt[names(All_gmt) %in% names(KEGG_gmt) | names(All_gmt) %in% names(GO_gmt)]
# filter_gmt <- All_gmt[names(All_gmt) %in% names(GO_gmt)]

## Filter genesets for very small/very big sizes (reduces multiple comparison deficit) (2918 genesets)
geneset_sizes <- unlist(lapply(GO_gmt, length))
geneset_indices <- which(geneset_sizes>=30 & geneset_sizes<200)
filtered_set <- immuno_gmt[geneset_indices]

library(tidyverse)
?contains()


vars_select(names(filtered_set), contains("VS"))
names(filtered_set) %>% contains("VS")

names(filtered_set)[grep("th17", ignore.case = T, names(filtered_set))] ### THIS MIGHT BE AN IDEA?!?!



length(filtered_set)

str(All_gmt)


names(All_gmt)[base:: grepl("REACTIVE_OXYGEN", names(All_gmt))]


## Perform camera analysis on filtered geneset
idx <- ids2indices(filtered_set, id = rownames(v))

camera_results <- camera(v, idx, design, contrast = contr.matrix[, "MSS_hiCIRCvsMSS"])
head(camera_results)
View(camera_results)
droplevels(subset(camera_results, PValue <= 0.001))


camera_results[grepl("REACTIVE_OXYGEN", rownames(camera_results)), ]

# BiocManager::install("qusage")
library(qusage)
camera_results_a <- camera_results[rownames(camera_results) %in% names(filtered_set),]

genesets_filtered <- idx
data_for_gs_analysis <- v

camera_descr <- unlist(lapply(rownames(camera_results_a), 
                              function(x){unlist(strsplit(x,"\\%"))[1]}))
camera_Phenotype <- unlist(lapply(camera_results_a[, "Direction"], 
                                  function(x){if(x=="Up"){1}else{(-1)}}))

camera_genes <- c()
for(i in 1:length(rownames(camera_results_a))){
  current_geneset <- unlist( 
    genesets_filtered[which(names(genesets_filtered) %in% 
                              rownames(camera_results_a)[i])])
  current_genes <- c()
  for(j in 1:length(current_geneset)){
    if(j==length(current_geneset)){
      current_genes <- paste(current_genes, 
                             rownames(data_for_gs_analysis)[current_geneset[j]],
                             sep = "")
    } else {
      current_genes <- paste(current_genes, 
                             rownames(data_for_gs_analysis)[current_geneset[j]], ",", 
                             sep = "")
    }
  }
  camera_genes <- rbind(camera_genes, current_genes)
}
rownames(camera_genes) <- rownames(camera_results_a)

camera_results_generic_em <- data.frame(rownames(camera_results_a), camera_descr, 
                                        PValue = camera_results_a[, "PValue"],
                                        FDR = camera_results_a[, "FDR"],
                                        camera_Phenotype,
                                        camera_genes)

camera_results_file <- "./Data/Genesets/camera_results_generic_hiCIRC_MSI.txt"
write.table(camera_results_generic_em, file.path(camera_results_file), 
            col.name = T, sep = "\t", row.names = F, quote = F)
# 
# camera_results_CD8 <- droplevels(subset(camera_results_generic_em, FDR <= 0.001))
# camera_results_CD8$camera_Phenotype <- ifelse((camera_results_CD8$camera_Phenotype == 1), "CD8 EMRA", "CD8 Naive")
# write.csv("../../Genesets/CD8_results.csv", x = camera_results_CD8, row.names = F)
# 
# camera_results_VD1$camera_Phenotype <- ifelse((camera_results_VD1$camera_Phenotype == 1), "VD1 CD27LO", "VD1 CD27HI")
# write.csv("../../Genesets/VD1_results.csv", x = camera_results_VD1, row.names = F)

expression_file <- "./Data/Genesets/expression_file.txt"
exp_fil <- as.data.frame(v$E)
write.table(exp_fil, file.path(expression_file),
            col.name = T, sep = "\t", row.names = F, quote = F)


#use easy cyRest library to communicate with cytoscape.
tryCatch(expr = { library(RCy3)}, 
         error = function(e) { install_github("cytoscape/RCy3")}, finally = library(RCy3))

#defined threshold for GSEA enrichments (need to be strings for cyrest call)
pvalue_threshold <- "0.5"
qvalue_threshold <- "0.1"

similarity_threshold <- "0.25"
similarity_metric <- "JACCARD"

# generic_gmt_file <- file.path(getwd(), gmt_file)
analysis_name <- "MSI_HvsMSS_hiCIRC"
cur_model_name <- paste("camera", analysis_name, sep="_")
results_filename <- file.path(getwd(),  camera_results_file)


current_network_name <- paste(cur_model_name, pvalue_threshold, qvalue_threshold, sep = "_")


results_filename <- "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/6_TCGA_Full/TCGA_FULL/Data/Genesets/camera_results_generic_hiCIRC_MSI.txt"
expression_filename <- "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/6_TCGA_Full/TCGA_FULL/Data/Genesets/expression_file.txt"
generic_gmt_file <- "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/6_TCGA_Full/TCGA_FULL/Data/Genesets/msigdb.v6.2.entrez.gmt"
em_command = paste('enrichmentmap build analysisType=generic',
                   "gmtFile=", generic_gmt_file,
                   "pvalue=", pvalue_threshold,
                   "qvalue=", qvalue_threshold,
                   "similaritycutoff=", similarity_threshold,
                   "coefficients=", similarity_metric,
                   "enrichmentsDataset1=", results_filename,
                   "expressionDataset1=", expression_filename)

# Above had sep = " "

#enrichment map command will return the suid of newly created network.
# install.packages("BiocManager")
# BiocManager::install("RCy3")
library(RCy3)
response <- commandsGET(em_command)

current_network_suid <- 0
#enrichment map command will return the suid of newly created network unless it Failed.  
#If it failed it will contain the word failed
if(grepl(pattern="Failed", response)){
  paste(response)
} else {
  current_network_suid <- response
}
response <- renameNetwork(current_network_name, as.numeric(current_network_suid))


# Bespoke
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

## Filter genesets for very small/very big sizes (reduces multiple comparison deficit) (4326 genesets)
geneset_sizes <- unlist(lapply(geneset_list, length))
geneset_indices <- which(geneset_sizes>=15 & geneset_sizes<200)
filtered_set <- geneset_list[geneset_indices]

th17 <- geneset_list[["Th17"]]
## Perform camera analysis on filtered geneset
idx <- ids2indices(th17, id = rownames(v))

camera_results <- camera(v, idx, design, contrast = contr.matrix[, "MSI_HvsMSS_hiCIRC"])
head(camera_results)
View(camera_results)
droplevels(subset(camera_results, FDR <= 0.01))

rm(listDF)
save.image(file = "./R_Data/Counts.RData")
