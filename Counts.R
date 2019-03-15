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


# Looking for the files
thousand.folders <- list.dirs(path = "./Data/Counts", full.names = T)
filelist1 <- sapply(thousand.folders[-1], function(x){
  list.files(x, pattern = "htseq.counts.gz$", full.names = T)})
filelist = unlist(filelist1)

# Read in files and combine
lists <- lapply(filelist, read.delim, header = F)
listsDF <- lists
lists <- listsDF

for (i in names(lists)){
  p <- gsub("./Data/Counts/", "", i)
  colnames(lists[[i]]) <- c("Gene", p)
}

names(lists) <- gsub("./Data/Counts/", "", names(lists))


# Patients I have CIRC scores for.
pat_sub <- read.csv("./Data/patient_subtypes.csv")
converter <- read.delim("./Data/Sample_Map.tsv")
converter$Patient.ID <- gsub("-", ".", converter$Case.ID)
converter1 <- converter[converter$Patient.ID %in% pat_sub$Patient.ID, ]
converter1$File.Name <- gsub(".htseq.counts.gz", "", converter1$File.Name)
# clin <- read.delim("./Data/Clinical/clinical.tsv")
# clin$Patient.ID <- gsub("-", ".", clin$submitter_id)
# clin1 <- clin[clin$Patient.ID %in% pat_sub$Patient.ID, ]

lists1 <- lists[names(lists) %in% converter1$File.ID]

multi_join <- function(list_of_loaded_data, join_func, ...){
  require("dplyr")
  output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
  return(output)
}

combined_df <- multi_join(lists1, full_join)

temp_df <- combined_df %>% gather(key = "File.ID", value = "Count", -Gene)
converter2 <- droplevels(subset(converter1, Sample.Type != "Solid Tissue Normal" &
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
Counts <- dcast(temp_df1, Gene ~ Patient.ID, sum, value.var = "Count")
# droplevels(subset(temp_df1, Patient.ID == "TCGA.A6.2672" & Gene == "ENSG00000000003.13")) # Matches the figure in "Try"


Counts_cleaned <- droplevels(Counts[!'%in%'(Counts$Gene, 
                                            c("__alignment_not_unique", "__no_feature",
                                              "__not_aligned", "__too_low_aQual", "__ambiguous")), ])


##### START EDGER
library(edgeR)
# Replenish "x"
row.names(Counts_cleaned) <- NULL
Counts_cleaned <- Counts_cleaned %>% column_to_rownames(., var = "Gene")
x <- DGEList(Counts_cleaned)

# Get groups
temp <- x$samples %>% rownames_to_column(., "Patient.ID")
colnames(pat_sub)[colnames(pat_sub) == "Subtype"] <- "group"
x$samples <- merge(temp[, c("Patient.ID", "lib.size", "norm.factors")], pat_sub[, c("Patient.ID", "group")], by = "Patient.ID") %>% column_to_rownames(., var = "Patient.ID")
group <- x$samples$group

samplenames <- colnames(x)

# Gaining second dataframe (Symbols)
# biocLite("Homo.sapiens", dependencies = T)
library(Homo.sapiens)
rownames(x) <- gsub("\\..*", "", rownames(x)) #Removes version from Ensembl gene ID
geneid <- rownames(x)

genes <- select(Homo.sapiens, keys = geneid, columns = c("SYMBOL", "TXCHROM", "ENTREZID"), 
                keytype = "ENSEMBL")
genes <- genes[!duplicated(genes$SYMBOL),]
x$genes <- genes

# Calculate counts per million, log counts per million (can also do RPKM [function rpkm()])
cpm <- cpm(x)
lcpm <- cpm(x, log = T)

# Remove lowly expressed transcripts
dim(x)
table(rowSums(x$counts==0)==10) # Show me the amount of transcripts that are zero for all 40 samples
CPM_scaling <- min(x$samples$lib.size)/1000000
cpm_scale <- 6.5/CPM_scaling

keep.exprs <- rowSums(cpm>cpm_scale)>=2 # Genes must have a cpm above 0.44 (count of 6.5 in lowest library) and be expressed in at least 2 groups (1 population)
x <- x[keep.exprs,, keep.lib.sizes = F]

# Normalising gene expression
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors
x <- estimateCommonDisp(x)
x <- estimateTagwiseDisp(x)


### Double check...
lcpm <- cpm(x, log = T)

## Online launch of this.
# biocLite("Glimma", dependencies = T)
# library(Glimma)
# glMDSPlot(lcpm, labels = paste(colnames(x), sep = "_"),
#           groups = x$samples[,c(2)], launch = T, top = 13564)
# col.cell <- c("#999999","#56B4E9","#E69F00","#009E73","#CC79A7")[x$samples$group]
# data.frame(x$samples$group, col.cell)
# pdf("../Figures/Proof/PCA_of_all_genes.pdf")
cbols
col.cell <- c("#999999","#56B4E9","#E69F00")[x$samples$group]
dim(lcpm)
plotMDS(lcpm, pch = 16, cex = 2, col = col.cell)
legend("top",
       fill = c("#999999", "#56B4E9",
                "#E69F00"),
       legend = levels(x$samples$group))


var_genes <- apply(lcpm, 1, var)


select_var <- names(sort(var_genes, decreasing = T))[1:150]

# Subset logcounts matrix
highly_variable_lcpm <- lcpm[select_var,]
# dim(highly_variable_lcpm)
# head(highly_variable_lcpm)

## Get some nicer colours
library(gplots)
library(RColorBrewer)
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
mycol <- colorpanel(1000,"blue","white","red")

## Get the SYMBOLS instead of ENTREZ ID
this <- rownames_to_column(as.data.frame(highly_variable_lcpm), var = "ENSEMBL")

this1 <- merge(this,x$genes, by = "ENSEMBL")
highly_variable_lcpm_sym <- within(this1, rm(TXCHROM, ENSEMBL))
distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method = "average")
col.cell1 <- c("#999999","#56B4E9","#E69F00")[x$samples$group]

# heatmap.2(as.matrix(highly_variable_lcpm_sym[, names(highly_variable_lcpm_sym) != "SYMBOL"]),
#           col = mycol,
#           trace = "none",
#           density.info = "none",
#           main = "B. Heatmap of the Top 150 most Variable Genes",
#           cex.main = 1.5,
#           #ColSideColors = col.cell,
#           scale = "row",
#           margin = c(10,5), lhei = c(2,10),
#           labCol = colnames(highly_variable_lcpm_sym),
#           labRow = highly_variable_lcpm_sym$SYMBOL,
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
# save.image(file = "half.RData")
# load("Bulk/Counts/half.RData")

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit, robust = T)
# plotSA(efit, main = "Final model: Mean−variance trend")
# dev.off()

## Number of differentially expressed genes
summary(decideTests(efit))

### More stringently selected DE genes
tfit <- treat(vfit, lfc = 1)
summary(decideTests(tfit))

dt <- decideTests(efit)
# summary(dt)

# de.common <- which(dt[,2]!=0 & dt[,3]!=0)
# length(de.common)


MSS_hiCIRC <- which(dt[, "MSS_hiCIRCvsMSS"] != 0)
Ensem_MSS <- efit$genes$SYMBOL[MSS_hiCIRC]


MSS_hiCIRC_ <- which(dt[, "MSI_HvsMSS_hiCIRC"] != 0)
Ensem_MSI <- efit$genes$ENSEMBL[MSS_hiCIRC_]
Symbol_MSI <- efit$genes$SYMBOL[MSS_hiCIRC_]
Symbol_MSI[grepl("ROR", Symbol_MSI)]


# i <- which(v$genes$ENSEMBL %in% Symbol_MSI)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale = "row",
#           labRow = v$genes$SYMBOL[i], labCol = colnames(v), #can also do labCol = groups
#           col = mycol, trace = "none", density.info = "none",
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column",
#           main = "Shared DE genes (Naive v Effectors)",
#           hclustfun = hclustAvg)
# dev.off()

MSI_H.vs.MSS_hiCIRC <- topTreat(efit, coef = 1, n = Inf)
MSI_H.vs.MSS <- topTreat(efit, coef = 2, n = Inf)
MSS_hiCIRC.vs.MSS <- topTreat(efit, coef = 3, n = Inf)


myData <- as.data.frame(MSI_H.vs.MSS_hiCIRC)
myData$padjThresh <- as.factor(myData$adj.P.Val < 0.05)
myData$Labels <- myData$SYMBOL
labelled_genes <- c("RORC", "TLR4")
myData$Labels[!myData$SYMBOL %in% labelled_genes] <- ""
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
All_gmt <- read.gmt("./Data/Genesets/msigdb.v6.2.entrez.gmt")
KEGG_gmt <- read.gmt("./Data/Genesets/c2.cp.kegg.v6.2.entrez.gmt")
GO_terms <- read.gmt("./Data/Genesets/c5.all.v6.2.entrez.gmt")


## Filter genesets that appear in only KEGG and GO databases (6103 genesets)
# filter_gmt <- All_gmt[names(All_gmt) %in% names(KEGG_gmt) | names(All_gmt) %in% names(GO_terms)]
filter_gmt <- All_gmt[names(All_gmt) %in% names(GO_terms)]

## Filter genesets for very small/very big sizes (reduces multiple comparison deficit) (4326 genesets)
geneset_sizes <- unlist(lapply(All_gmt, length))
geneset_indices <- which(geneset_sizes>=15 & geneset_sizes<200)
filtered_set <- All_gmt[geneset_indices]

head(filtered_set)

## Perform camera analysis on filtered geneset
idx <- ids2indices(filtered_set, id = rownames(v))

camera_results <- camera(v, idx, design, contrast = contr.matrix[, "MSI_HvsMSS_hiCIRC"])
nrow(camera_results)
droplevels(subset(camera_results, FDR <= 0.001))


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

camera_results_file <- "../../Genesets/camera_results_generic_CD8.txt"
write.table(camera_results_generic_em, file.path(camera_results_file), 
            col.name = T, sep = "\t", row.names = F, quote = F)
# 
# camera_results_CD8 <- droplevels(subset(camera_results_generic_em, FDR <= 0.001))
# camera_results_CD8$camera_Phenotype <- ifelse((camera_results_CD8$camera_Phenotype == 1), "CD8 EMRA", "CD8 Naive")
# write.csv("../../Genesets/CD8_results.csv", x = camera_results_CD8, row.names = F)
# 
# camera_results_VD1$camera_Phenotype <- ifelse((camera_results_VD1$camera_Phenotype == 1), "VD1 CD27LO", "VD1 CD27HI")
# write.csv("../../Genesets/VD1_results.csv", x = camera_results_VD1, row.names = F)

expression_file <- "../../Genesets/expression_file.txt"
exp_fil <- as.data.frame(v$E)
write.table(exp_fil, file.path(expression_file),
            col.name = T, sep = "\t", row.names = F, quote = F)


#use easy cyRest library to communicate with cytoscape.
tryCatch(expr = { library(RCy3)}, 
         error = function(e) { install_github("cytoscape/RCy3")}, finally = library(RCy3))

#defined threshold for GSEA enrichments (need to be strings for cyrest call)
pvalue_threshold <- "0.05"
qvalue_threshold <- "0.001"

similarity_threshold <- "0.25"
similarity_metric <- "JACCARD"

# generic_gmt_file <- file.path(getwd(), gmt_file)
analysis_name <- "EMRA_vs_Naive"
cur_model_name <- paste("camera", analysis_name, sep="_")
results_filename <- file.path(getwd(),  camera_results_file)


current_network_name <- paste(cur_model_name, pvalue_threshold, qvalue_threshold, sep = "_")


results_filename <- "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/Cyto/camera_results_generic_CD8.txt"
expression_filename <- "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/Cyto/expression_file.txt"
generic_gmt_file <- "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/4_Gamma_Delta/Cyto/msigdb.v6.2.entrez.gmt"
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

