# A script to investigate the differences in immunome between MSI-H and MSS-hiCIRC

## Packages
required <- c("tidyverse", "ggpubr", "ggbiplot", "devtools", "gplots", "UsefulFunctions")
for (lib in required)
{
  if (!require(lib, character.only = T))
  {
    install.packages(lib)
    suppressMessages(library(lib, character.only = T, quietly = T))
  }
}

# Comparisons and Colours
my_comparisons <- list(c("MSS-hiCIRC", "MSI-H"),
                       c("MSS-hiCIRC", "MSS"),
                       c("MSI-H", "MSS"))
cbcols <- c("MSS-hiCIRC" = "#999999",
            "MSI-H" = "#56B4E9",
            "MSS" = "#009E73")


## Performing differential gene expression ##
load("./R_Data/Counts_clean.RData")
# rm(list = setdiff(ls(), c("cqn_Counts"))) # Clean environment

##### START EDGER
BiocManager::install("edgeR")
library(edgeR)

# Replenish "x"
Counts_clean <- cqn_Counts$counts


# Make rownames ENTREZID again
ensembl_DB <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
Gene_Map <- getBM(attributes = c("entrezgene_id", "hgnc_symbol"),
                  filters = "hgnc_symbol", values = rownames(Counts_clean), ensembl_DB)
head(Gene_Map) # Might throw warnings if a mirror is down

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

# # Remove lowly expressed transcripts #### DON'T DO THIS, REMOVES GENES TAHT ARE INCLUDED ####
# dim(x)
# table(rowSums(x$counts==0)==10) # Show me the amount of transcripts that are zero for 10 samples
# CPM_scaling <- min(x$samples$lib.size)/1000000
# cpm_scale <- 6.5/CPM_scaling
# 
# keep.exprs <- rowSums(cpms>cpm_scale)>=1 # Genes must count of 6.5 in lowest library
# x <- x[keep.exprs,, keep.lib.sizes = F]

# Normalising gene expression
x <- calcNormFactors(x, method = "TMM")
# x$samples$norm.factors
x <- estimateCommonDisp(x) ###########ERROR
x <- estimateTagwiseDisp(x)

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

### Removing heteroscedascity
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
Genes_MSI_hiCIRC <- names(MSS_hiCIRC_MSI)

Genes_MSI_hiCIRC[grepl("ROR", Genes_MSI_hiCIRC)]
Genes_MSI_hiCIRC[grepl("IL17", Genes_MSI_hiCIRC)]
Genes_MSI_hiCIRC[grepl("CCL", Genes_MSI_hiCIRC)]
Genes_MSI_hiCIRC[grepl("IL", Genes_MSI_hiCIRC)]
Genes_MSI_hiCIRC[grepl("CCR", Genes_MSI_hiCIRC)]
Genes_MSI_hiCIRC[grepl("TLR", Genes_MSI_hiCIRC)]
Genes_MSI_hiCIRC[grepl("MUC", Genes_MSI_hiCIRC)]


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
# All_gmt <- read.gmt("./Data/Genesets/GSEA/Symbol/msigdb.v7.0.symbols.gmt")
KEGG_gmt <- read.gmt("./Data/Genesets/GSEA/Symbol/c2.cp.kegg.v7.0.symbols.gmt")
# GO_gmt <- read.gmt("./Data/Genesets/GSEA/Symbol/c5.all.v7.0.symbols.gmt")
# immuno_gmt <- read.gmt("./Data/Genesets/GSEA/Symbol/c7.all.v7.0.symbols.gmt")
hallmark_gmt <- read.gmt("./Data/Genesets/GSEA/Symbol/h.all.v7.0.symbols.gmt")
bio_gmt <- read.gmt("./Data/Genesets/GSEA/Symbol/c5.bp.v7.1.symbols.gmt")

filtered_set <- read.delim("./Data/Genesets/Filtered_set_PCGSE.txt")$x %>% c()

bio_filter <- bio_gmt[names(bio_gmt) %in% filtered_set]



#### ENSURE 03/06/09 is MARCHF6 ###
GoI <- read.csv("./Output/Genesets/Genesets_of_interest.csv", stringsAsFactors = F)

GoI1 <- droplevels(subset(GoI, Type_of_data == "CellType"))[, c("CellType", "Symbol")]

imm_list <- split(GoI1$Symbol, GoI1$CellType)


GoI1$Symbol[!('%in%'(GoI1$Symbol, rownames(v)))]

## Filter genesets that appear in only KEGG and GO databases (6103 genesets)
# filter_gmt <- All_gmt[names(All_gmt) %in% names(KEGG_gmt) | names(All_gmt) %in% names(GO_gmt)]
# filter_gmt <- All_gmt[names(All_gmt) %in% names(GO_gmt)]

## Filter genesets for very small/very big sizes (reduces multiple comparison deficit) (2918 genesets)
geneset_sizes <- unlist(lapply(bio_filter, length))
geneset_sizes

geneset_indices <- which(geneset_sizes>=50 & geneset_sizes<200)
filtered_set <- bio_gmt[geneset_indices]
filtered_set1 <- filtered_set[!'%in%'(names(filtered_set), c("Th17_origin", "SW480 cancer cells", "T cells", "T helper cells", "Cytotoxic cells"))]

hallmark_gmt1 <- hallmark_gmt[!'%in%'(names(hallmark_gmt), c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_MITOTIC_SPINDLE",
                                                             "HALLMARK_MYOGENESIS", "HALLMARK_ESTROGEN_RESPONSE_EARLY",
                                                             "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_ESTROGEN_RESPONSE_LATE",
                                                             "HALLMARK_HEME_METABOLISM", "HALLMARK_BILE_ACID_METABOLISM",
                                                             "HALLMARK_SPERMATOGENESIS", "HALLMARK_PANCREAS_BETA_CELLS",
                                                             "HALLMARK_G2M_CHECKPOINT", "HALLMARK_ANDROGEN_RESPONSE", 
                                                             "HALLMARK_APICAL_SURFACE", "HALLMARK_XENOBIOTIC_METABOLISM"
                                                             ))]

## Perform camera analysis on filtered geneset



idx <- ids2indices(bio_filter, id = rownames(v))
camera_results <- camera(v, idx, design, contrast = contr.matrix[, "MSI_HvsMSS_hiCIRC"])

## Use the camera_result table?

droplevels(subset(camera_results, FDR < 0.05)) %>% View()



# BiocManager::install("qusage")
library(qusage)
camera_results_a <- camera_results[rownames(camera_results) %in% names(filtered_set1),]

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

camera_results_file <- "./Output/Genesets/DGE_Cyto/camera_results_generic_hiCIRC_MSI.txt"
write.table(camera_results_generic_em, file.path(camera_results_file), 
            col.name = T, sep = "\t", row.names = F, quote = F)




# Getting it ready for Cytoscape running...


expression_file <- "./Output/Genesets/DGE_Cyto/expression_file.txt"
exp_fil <- as.data.frame(v$E)
write.table(exp_fil, file.path(expression_file),
            col.name = T, sep = "\t", row.names = F, quote = F)


#use easy cyRest library to communicate with cytoscape.
tryCatch(expr = { library(RCy3)}, 
         error = function(e) { install_github("cytoscape/RCy3")}, finally = library(RCy3))

#defined threshold for GSEA enrichments (need to be strings for cyrest call)
pvalue_threshold <- "0.05"
qvalue_threshold <- "0.05"

similarity_threshold <- "0.25"
similarity_metric <- "JACCARD"

# generic_gmt_file <- file.path(getwd(), gmt_file)
analysis_name <- "MSI_H_vs_MSS-hiCIRC"
cur_model_name <- paste("camera", analysis_name, sep="_")
results_filename <- file.path(getwd(),  camera_results_file)


current_network_name <- paste(cur_model_name, pvalue_threshold, qvalue_threshold, sep = "_")



results_filename <- "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/6_TCGA_Full/TCGA_FULL/Output/Genesets/DGE_Cyto/camera_results_generic_hiCIRC_MSI.txt"
expression_filename <- "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/6_TCGA_Full/TCGA_FULL/Output/Genesets/DGE_Cyto/expression_file.txt"
generic_gmt_file <- "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/6_TCGA_Full/TCGA_FULL/Output/Genesets/Genesets_of_interest.csv"
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



### Barcode plot
barcodeplot(tfit$t[, "MSI_HvsMSS_hiCIRC"], index = idx$Th17_extend, main = ".")
barcodeplot(tfit$t[, "MSI_HvsMSS_hiCIRC"], index = idx$`B cells`, main = ".")
barcodeplot(tfit$t[, "MSI_HvsMSS_hiCIRC"], index = idx$`Mast cells`, main = ".")

tfit$t

heatmap.2(tfit$t[, "MSI_HvsMSS_hiCIRC"][idx$Th17_extend])

heatmap.2()












