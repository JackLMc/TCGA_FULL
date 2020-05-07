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
x <- estimateCommonDisp(x)
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
# KEGG_gmt <- read.gmt("./Data/Genesets/GSEA/Symbol/c2.cp.kegg.v7.0.symbols.gmt")
# GO_gmt <- read.gmt("./Data/Genesets/GSEA/Symbol/c5.all.v7.0.symbols.gmt")
# immuno_gmt <- read.gmt("./Data/Genesets/GSEA/Symbol/c7.all.v7.0.symbols.gmt")
# hallmark_gmt <- read.gmt("./Data/Genesets/GSEA/Symbol/h.all.v7.0.symbols.gmt")
#### ENSURE 03/06/09 is MARCHF6 ###
GoI <- read.csv("./Output/Genesets/Genesets_of_interest.csv", stringsAsFactors = F)

GoI1 <- droplevels(subset(GoI, Type_of_data == "CellType"))[, c("CellType", "Symbol")]

imm_list <- split(GoI1$Symbol, GoI1$CellType)


GoI1$Symbol[!('%in%'(GoI1$Symbol, rownames(v)))]

## Filter genesets that appear in only KEGG and GO databases (6103 genesets)
# filter_gmt <- All_gmt[names(All_gmt) %in% names(KEGG_gmt) | names(All_gmt) %in% names(GO_gmt)]
# filter_gmt <- All_gmt[names(All_gmt) %in% names(GO_gmt)]

## Filter genesets for very small/very big sizes (reduces multiple comparison deficit) (2918 genesets)
geneset_sizes <- unlist(lapply(imm_list, length))
geneset_indices <- which(geneset_sizes>=15 & geneset_sizes<200)
filtered_set <- imm_list[geneset_indices]
head(filtered_set)
filtered_set1 <- filtered_set[!'%in%'(names(filtered_set), c("Th17_origin", "SW480 cancer cells", "T cells", "T helper cells", "Cytotoxic cells"))]


## Perform camera analysis on filtered geneset
idx <- ids2indices(filtered_set1, id = rownames(v))
camera_results <- camera(v, idx, design, contrast = contr.matrix[, "MSI_HvsMSS_hiCIRC"])
head(camera_results)
View(camera_results)
droplevels(subset(camera_results, PValue <= 0.01))








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




















############################


CD27LO.vs.Naive <- topTreat(tfit, coef = 6, n = Inf)
CD27LO.vs.VD2 <- topTreat(tfit, coef = 7, n = Inf)

EMRA.vs.Naive <- topTreat(tfit, coef = 8, n = Inf)
# write.csv(EMRA.vs.Naive, file = "./Bulk/Output/EMRA.vs.Naive.csv")

EMRA.vs.VD2 <- topTreat(tfit, coef = 9, n = Inf)

Naive.vs.VD2 <- topTreat(tfit, coef = 10, n = Inf)


######### AME
sig <- droplevels(subset(CD27LO.vs.CD27HI, adj.P.Val <= 0.01))
sig$Direction <- as.factor(ifelse((sig$logFC < 0), "CD27HI", "CD27LO"))

CD27HI <- droplevels(subset(sig, Direction == "CD27HI"))
genes_27HI <- CD27HI$ENTREZID

CD27LO <- droplevels(subset(sig, Direction == "CD27LO"))
genes_27LO <- CD27LO$ENTREZID

library(biomaRt)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

# CD27HI
CD27HI_genes <- getBM(attributes = c("transcript_start", "transcript_end", "chromosome_name", "ensembl_gene_id"),
                      filters = "entrezgene", values = genes_27HI, mart = mart)


CD27HI_genes$chromosome_name_right <- paste("chr", CD27HI_genes$chromosome_name, sep = "") 
CD27HI_genes$chromosome_name_right <- as.factor(CD27HI_genes$chromosome_name_right)
CD27HI_genes1 <- CD27HI_genes[!grepl("chrCHR", CD27HI_genes$chromosome_name_right),]

CD27HI_ensembl1 <- CD27HI_genes1[!duplicated(CD27HI_genes1$ensembl_gene_id), ]

CD27HI_ensembl1$ensembl_gene_id <- as.factor(CD27HI_ensembl1$ensembl_gene_id)

CD27HI_ensembl2 <- data.frame(ensembl_gene = character(),
                              transcript_start = double(),
                              transcript_end= double(),
                              chromosome_name_right = character(),
                              stringsAsFactors = F)
c <- 1
for(i in levels(CD27HI_ensembl1$ensembl_gene_id)){
  print(i)
  working <- droplevels(subset(CD27HI_ensembl1, ensembl_gene_id == i))
  avg_start <- mean(working$transcript_start)
  avg_end <- mean(working$transcript_end)
  working$chromosome_name_right <- as.factor(working$chromosome_name_right)
  CD27HI_ensembl2[c, "ensembl_gene"] <- i
  CD27HI_ensembl2[c, "transcript_start"] <- avg_start
  CD27HI_ensembl2[c, "transcript_end"] <- avg_end
  CD27HI_ensembl2[c, "chromosome_name_right"] <- as.character(levels(working$chromosome_name_right))
  c <- c + 1
}

CD27HI_ensembl2$transcript_start_ext <- CD27HI_ensembl2$transcript_start - 2000
CD27HI_ensembl2$transcript_end_ext <- CD27HI_ensembl2$transcript_end + 1000
library(GenomicRanges)

gr_27HI <- GRanges(seqnames = Rle(CD27HI_ensembl2$chromosome_name_right),
                   ranges = IRanges(CD27HI_ensembl2$transcript_start_ext,
                                    end = CD27HI_ensembl2$transcript_end_ext),
                   strand = Rle(strand(c(rep("*", length(CD27HI_ensembl2$chromosome_name_right))))),
                   names = CD27HI_ensembl2$ensembl_gene)

df_27HI <- data.frame(seqnames = seqnames(gr_27HI),
                      starts = start(gr_27HI)-1,
                      ends = end(gr_27HI),
                      names = gr_27HI$names
)

df_27HI1 <- droplevels(subset(df_27HI, names != "ENSG00000189283"))

subset(df_27HI, starts == 63637707)

write.table(df_27HI1, file = "Bulk/Output/CD27HI_list1.bed", quote = F, sep = "\t", row.names = F, col.names = F)


library(Homo.sapiens)
promoters(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), 2000, -2000)


# CD27LO
CD27LO_genes <- getBM(attributes = c("transcript_start", "transcript_end", "chromosome_name", "ensembl_gene_id"),
                      filters = "entrezgene", values = genes_27LO, mart = mart)


CD27LO_genes$chromosome_name_right <- paste("chr", CD27LO_genes$chromosome_name, sep = "") 
CD27LO_genes$chromosome_name_right <- as.factor(CD27LO_genes$chromosome_name_right)
CD27LO_genes1 <- CD27LO_genes[!grepl("chrCHR", CD27LO_genes$chromosome_name_right),]

CD27LO_ensembl1 <- CD27LO_genes1[!duplicated(CD27LO_genes1$ensembl_gene_id), ]

CD27LO_ensembl1$ensembl_gene_id <- as.factor(CD27LO_ensembl1$ensembl_gene_id)

CD27LO_ensembl2 <- data.frame(ensembl_gene = character(),
                              transcript_start = double(),
                              transcript_end= double(),
                              chromosome_name_right = character(),
                              stringsAsFactors = F)
c <- 1
for(i in levels(CD27LO_ensembl1$ensembl_gene_id)){
  print(i)
  working <- droplevels(subset(CD27LO_ensembl1, ensembl_gene_id == i))
  avg_start <- mean(working$transcript_start)
  avg_end <- mean(working$transcript_end)
  working$chromosome_name_right <- as.factor(working$chromosome_name_right)
  CD27LO_ensembl2[c, "ensembl_gene"] <- i
  CD27LO_ensembl2[c, "transcript_start"] <- avg_start
  CD27LO_ensembl2[c, "transcript_end"] <- avg_end
  CD27LO_ensembl2[c, "chromosome_name_right"] <- as.character(levels(working$chromosome_name_right))
  c <- c + 1
}

CD27LO_ensembl2$transcript_start_ext <- CD27LO_ensembl2$transcript_start - 2000
CD27LO_ensembl2$transcript_end_ext <- CD27LO_ensembl2$transcript_end + 1000
library(GenomicRanges)

gr_27LO <- GRanges(seqnames = Rle(CD27LO_ensembl2$chromosome_name_right),
                   ranges = IRanges(CD27LO_ensembl2$transcript_start_ext,
                                    end = CD27LO_ensembl2$transcript_end_ext),
                   strand = Rle(strand(c(rep("*", length(CD27LO_ensembl2$chromosome_name_right))))),
                   names = CD27LO_ensembl2$ensembl_gene)

df_27LO <- data.frame(seqnames = seqnames(gr_27LO),
                      starts = start(gr_27LO)-1,
                      ends = end(gr_27LO),
                      names = gr_27LO$names
)

write.table(df_27LO, file = "Bulk/Output/CD27LO_list.bed", quote = F, sep = "\t", row.names = F, col.names = F)

### end ####

length(df_27LO$seqnames)

write.table(EMRA.vs.Naive, "Bulk/Output/EMRA_vs_Naive_DGE.txt", sep = "\t", quote = F, row.names = F)

# These are the outputs of EdgeR!!
# head(CD27LO.vs.EMRA)
# head(CD27HI.vs.VD2)
# head(CD27HI.vs.EMRA)
# head(CD27HI.vs.Naive)

## Highlight this for notable genes... (CD28 one example)
pdf("Bulk/Figures/Volcano_Bulk.pdf")
plotMD(efit, column = 1, status = dt[,1], main = colnames(tfit)[1], 
       xlim = c(-2,13), col = c("#009E73", "#999999"))
legend("topright", fill = c("#999999", "black", "#009E73"),
       legend = c("VD1 CD27LO", "No Change", "VD1 CD27HI"))

dev.off()
glMDPlot(tfit, coef = 1, status = dt, main = colnames(efit)[1],
         side.main = "ENTREZID", counts = x$counts, groups = group, launch = T)

# Heatmaps
## Change the top genes, for various heatmaps looking at different sets
library(gplots)

# CD27HI versus CD27LO
CD27LO.vs.CD27HI1 <- CD27LO.vs.CD27HI[!grepl("^TRAV", CD27LO.vs.CD27HI$SYMBOL), ]
CD27LO.vs.CD27HI2 <- CD27LO.vs.CD27HI1[!grepl("^TRBV", CD27LO.vs.CD27HI1$SYMBOL), ]
CD27LO.vs.CD27HI2 <- Take_Sigs(CD27LO.vs.CD27HI2)

CD27LO.vs.CD27HI.topgenes <- CD27LO.vs.CD27HI2$ENTREZID[1:100]
# i <- which(v$genes$ENTREZID %in% CD27LO.vs.CD27HI.topgenes)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale = "row",
#           labRow = v$genes$SYMBOL[i], labCol = group, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column",
#           main = "CD27LO vs CD27HI")

# gaining a dataframe with the differentially expressed
# this <- CD27LO.vs.CD27HI2[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "CD27HI", "CD27LO") # THIS WILL CHANGE ONCE YOU SET THE LFC = 2
# that <- droplevels(subset(this, Group == "CD27HI"))
# 
# CD27LO_CD27HI_DEgenes <- this
# writeCsvO(CD27LO_CD27HI_DEgenes)

# CD27HI versus EMRA
CD27HI.vs.EMRA1 <- CD27HI.vs.EMRA[!grepl("^TRAV", CD27HI.vs.EMRA$SYMBOL),]
CD27HI.vs.EMRA2 <- CD27HI.vs.EMRA1[!grepl("^TRBV", CD27HI.vs.EMRA1$SYMBOL),]
CD27HI.vs.EMRA2 <- Take_Sigs(CD27HI.vs.EMRA2)

CD27HI.vs.EMRA.topgenes <- CD27HI.vs.EMRA2$ENTREZID[1:100]
# i <- which(v$genes$ENTREZID %in% CD27HI.vs.EMRA.topgenes)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale = "row",
#           labRow = v$genes$SYMBOL[i], labCol = group, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column", 
#           main = "CD27HI vs EMRA")
## gaining a dataframe with the differentially expressed
# this <- CD27HI.vs.EMRA2[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "EMRA", "CD27HI")
# 
# CD27HI_EMRA_DEgenes <- this
# writeCsvO(CD27HI_EMRA_DEgenes)

# CD27HI versus Naive
CD27HI.vs.Naive1 <- CD27HI.vs.Naive[!grepl("^TRAV", CD27HI.vs.Naive$SYMBOL),]
CD27HI.vs.Naive2 <- CD27HI.vs.Naive1[!grepl("^TRBV", CD27HI.vs.Naive1$SYMBOL),]
CD27HI.vs.Naive2 <- Take_Sigs(CD27HI.vs.Naive2)

CD27HI.vs.Naive.topgenes <- CD27HI.vs.Naive2$ENTREZID[1:100]
# i <- which(v$genes$ENTREZID %in% CD27HI.vs.Naive.topgenes)
# 
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale = "row",
#           labRow = v$genes$SYMBOL[i], labCol = group, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column",
#           main = "CD27HI vs Naive")

## gaining a dataframe with the differentially expressed
# this <- CD27HI.vs.Naive2[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "Naive", "CD27HI")
# CD27HI_Naive_DEgenes <- this
# writeCsvO(CD27HI_Naive_DEgenes)


# CD27HI versus VD2
CD27HI.vs.VD21 <- CD27HI.vs.VD2[!grepl("^TRAV", CD27HI.vs.VD2$SYMBOL),]
CD27HI.vs.VD22 <- CD27HI.vs.VD21[!grepl("^TRBV", CD27HI.vs.VD21$SYMBOL),]
CD27HI.vs.VD22 <- Take_Sigs(CD27HI.vs.VD22)

CD27HI.vs.VD2.topgenes <- CD27HI.vs.VD22$ENTREZID[1:100]
# i <- which(v$genes$ENTREZID %in% CD27HI.vs.VD2.topgenes)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale = "row",
#           labRow = v$genes$SYMBOL[i], labCol = group, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column", 
#           main = "CD27HI vs VD2")

## gaining a dataframe with the differentially expressed
# this <- CD27HI.vs.VD22[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "VD2", "CD27HI")
# CD27HI_VD2_DEgenes <- this
# writeCsvO(CD27HI_VD2_DEgenes)


# CD27LO versus VD2
CD27LO.vs.VD21 <- CD27LO.vs.VD2[!grepl("^TRAV", CD27LO.vs.VD2$SYMBOL),]
CD27LO.vs.VD22 <- CD27LO.vs.VD21[!grepl("^TRBV", CD27LO.vs.VD21$SYMBOL),]
CD27LO.vs.VD22 <- Take_Sigs(CD27LO.vs.VD22)

CD27LO.vs.VD2.topgenes <- CD27LO.vs.VD22$ENTREZID[1:100]
# i <- which(v$genes$ENTREZID %in% CD27LO.vs.VD2.topgenes)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale = "row",
#           labRow = v$genes$SYMBOL[i], labCol = group, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column",
#           main = "CD27LO vs VD2")

## gaining a dataframe with the differentially expressed
# this <- CD27LO.vs.VD22[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "VD2", "CD27LO")
# CD27LO_VD2_DEgenes <- this
# writeCsvO(CD27LO_VD2_DEgenes)


# CD27LO versus EMRA
CD27LO.vs.EMRA1 <- CD27LO.vs.EMRA[!grepl("^TRAV", CD27LO.vs.EMRA$SYMBOL),]
CD27LO.vs.EMRA2 <- CD27LO.vs.EMRA1[!grepl("^TRBV", CD27LO.vs.EMRA1$SYMBOL),]
CD27LO.vs.EMRA2 <- Take_Sigs(CD27LO.vs.EMRA2)

CD27LO.vs.EMRA.topgenes <- CD27LO.vs.EMRA2$ENTREZID[1:100]
# i <- which(v$genes$ENTREZID %in% CD27LO.vs.EMRA.topgenes)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale = "row",
#           labRow = v$genes$SYMBOL[i], labCol = group, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column",
#           main = "CD27LO vs EMRA")

## gaining a dataframe with the differentially expressed
# this <- CD27LO.vs.EMRA2[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "EMRA", "CD27LO")
# CD27LO_EMRA_DEgenes <- this
# writeCsvO(CD27LO_EMRA_DEgenes)


# CD27LO versus Naive
CD27LO.vs.Naive1 <- CD27LO.vs.Naive[!grepl("^TRAV", CD27LO.vs.Naive$SYMBOL),]
CD27LO.vs.Naive2 <- CD27LO.vs.Naive1[!grepl("^TRBV", CD27LO.vs.Naive1$SYMBOL),]
CD27LO.vs.Naive2 <- Take_Sigs(CD27LO.vs.Naive2)

CD27LO.vs.Naive.topgenes <- CD27LO.vs.Naive2$ENTREZID[1:100]
# i <- which(v$genes$ENTREZID %in% CD27LO.vs.Naive.topgenes)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale = "row",
#           labRow = v$genes$SYMBOL[i], labCol = group, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column", 
#           main = "CD27LO vs Naive")

## gaining a dataframe with the differentially expressed
# this <- CD27LO.vs.Naive2[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "Naive", "CD27LO")
# CD27LO_Naive_DEgenes <- this
# writeCsvO(CD27LO_Naive_DEgenes)


# EMRA versus Naive
EMRA.vs.Naive1 <- EMRA.vs.Naive[!grepl("^TRAV", EMRA.vs.Naive$SYMBOL),]
EMRA.vs.Naive2 <- EMRA.vs.Naive1[!grepl("^TRBV", EMRA.vs.Naive1$SYMBOL),]
EMRA.vs.Naive2 <- Take_Sigs(EMRA.vs.Naive2)

EMRA.vs.Naive.topgenes <- EMRA.vs.Naive2$ENTREZID[1:100]
# i <- which(v$genes$ENTREZID %in% EMRA.vs.Naive.topgenes)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale = "row",
#           labRow = v$genes$SYMBOL[i], labCol = group, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column", 
#           main = "EMRA vs Naive")

## gaining a dataframe with the differentially expressed
# this <- EMRA.vs.Naive2[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "Naive", "EMRA")
# EMRA_Naive_DEgenes <- this
# writeCsvO(EMRA_Naive_DEgenes)


# EMRA versus VD2
EMRA.vs.VD21 <- EMRA.vs.VD2[!grepl("^TRAV", EMRA.vs.VD2$SYMBOL),]
EMRA.vs.VD22 <- EMRA.vs.VD21[!grepl("^TRBV", EMRA.vs.VD21$SYMBOL),]
EMRA.vs.VD22 <- Take_Sigs(EMRA.vs.VD22)

EMRA.vs.VD2.topgenes <- EMRA.vs.VD22$ENTREZID[1:100]
# i <- which(v$genes$ENTREZID %in% EMRA.vs.VD2.topgenes)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale = "row",
#           labRow = v$genes$SYMBOL[i], labCol = group, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column", 
#           main = "EMRA vs VD2")

## gaining a dataframe with the differentially expressed
# this <- EMRA.vs.VD22[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "VD2", "EMRA")
# EMRA_VD2_DEgenes <- this
# writeCsvO(EMRA_VD2_DEgenes)


# Naive versus VD2
Naive.vs.VD21 <- Naive.vs.VD2[!grepl("^TRAV", Naive.vs.VD2$SYMBOL),]
Naive.vs.VD22 <- Naive.vs.VD21[!grepl("^TRBV", Naive.vs.VD21$SYMBOL),]
Naive.vs.VD22 <- Take_Sigs(Naive.vs.VD22)

Naive.vs.VD2.topgenes <- Naive.vs.VD22$ENTREZID[1:100]
# i <- which(v$genes$ENTREZID %in% Naive.vs.VD2.topgenes)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale = "row",
#           labRow = v$genes$SYMBOL[i], labCol = group, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column", 
#           main = "Naive vs VD2")

## gaining a dataframe with the differentially expressed
# this <- Naive.vs.VD22[1:100, c("SYMBOL", "logFC", "adj.P.Val")]
# this$Group <- ifelse((this$logFC<0), "VD2", "Naive")
# Naive_VD2_DEgenes <- this
# writeCsvO(Naive_VD2_DEgenes)


# Shared top 100 DE genes
k <- which(CD27LO.vs.CD27HI.topgenes %in% EMRA.vs.Naive.topgenes)
j <- CD27LO.vs.CD27HI.topgenes[k]
i <- which(v1$genes$ENTREZID %in% j)
mycol <- colorpanel(1000,"blue","white","red")

col.cell1 <- c("#56B4E9","#E69F00","#009E73","#999999")[v1$targets$group]
data.frame(v1$targets$group, col.cell1)

# png(filename = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/Figures/Paper/Shared_top100_DEgenes.png",
#     width = 300, height = 300, units = "mm", res = 300)
heatmap.2(v1$E[i,], scale = "row",
          labRow = v1$genes$SYMBOL[i], labCol = group1, 
          col = mycol, trace = "none", density.info = "none", 
          margin = c(8,6), lhei = c(2,10), ColSideColors = col.cell1)
par(xpd = T)
legend(x = 0.87, y = 1.05, 
       fill = c("#999999", "#56B4E9",
                "#E69F00", "#009E73"),
       legend = levels(v1$targets$group))
# dev.off()


# # Check other way for %in% (same)
# k <- which(EMRA.vs.Naive.topgenes  %in% CD27LO.vs.CD27HI.topgenes)
# j <- EMRA.vs.Naive.topgenes[k]
# i <- which(v1$genes$ENTREZID %in% j)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v1$E[i,], scale = "row",
#           labRow = v1$genes$SYMBOL[i], labCol = group1, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column",
#           main = "Shared top 100 DE genes (Naive v Effectors)")



# Find top 10 genes in all combinations
# colnames(contr.matrix)
# 
# this <- as.vector(rbind(CD27LO.vs.CD27HI2$ENTREZID[1:25], 
#                         CD27HI.vs.EMRA2$ENTREZID[1:25], 
#                         CD27HI.vs.Naive2$ENTREZID[1:25],
#                         CD27HI.vs.Naive2$ENTREZID[1:25],
#                         CD27HI.vs.VD22$ENTREZID[1:25],
#                         CD27LO.vs.EMRA2$ENTREZID[1:25],
#                         CD27LO.vs.Naive2$ENTREZID[1:25],
#                         CD27LO.vs.VD22$ENTREZID[1:25],
#                         EMRA.vs.Naive2$ENTREZID[1:25],
#                         EMRA.vs.VD22$ENTREZID[1:25],
#                         Naive.vs.VD22$ENTREZID[1:25]))
# this <- this[!duplicated(this)]
# 
# i <- which(v$genes$ENTREZID %in% this)
# mycol <- colorpanel(1000,"blue","white","red")
# heatmap.2(v$E[i,], scale = "row",
#           labRow = v$genes$SYMBOL[i], labCol = group, 
#           col = mycol, trace = "none", density.info = "none", 
#           margin = c(8,6), lhei = c(2,10), #dendrogram = "column", 
#           main = "Top 20 DE genes across all populations")
# length(this)

#### Write the common genes



# Geneset Enrichment Analysis
CD27LO.vs.CD27HI.topgenes

# Hallmark_genesets
load("~/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/RData_Objects/Hallmark_genesets.rdata")
idx <- ids2indices(Hs.H, id = rownames(v))
cam.CD27LO.vs.CD27HI <- camera(v, idx, design, contrast = contr.matrix[, "CD27LOvsCD27HI"])
head(cam.CD27LO.vs.CD27HI, 5)

sig.CD27LO.vs.CD27HI <- droplevels(subset(cam.CD27LO.vs.CD27HI, FDR <= 0.01))
sig.CD27LO.vs.CD27HI

cam.EMRA.vs.Naive <- camera(v, idx, design, contrast = contr.matrix[, "EMRAvsNaive"])
head(cam.EMRA.vs.Naive, 5)
sig.EMRA.vs.Naive <- droplevels(subset(cam.EMRA.vs.Naive, FDR <= 0.01))
sig.EMRA.vs.Naive




length(Hs.H[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]])


# KEGG
load("~/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/RData_Objects/kegg_human.rdata")

## Geneset enrichment
idx <- ids2indices(kegg_human, id = rownames(v))
cam.CD27LO.vs.CD27HI <- camera(v, idx, design, contrast = contr.matrix[, "CD27LOvsCD27HI"])
head(cam.CD27LO.vs.CD27HI, 5)
sig.CD27LO.vs.CD27HI <- droplevels(subset(cam.CD27LO.vs.CD27HI, FDR <= 0.05))
KEGG_GD <- rownames_to_column(sig.CD27LO.vs.CD27HI, var = "KEGG_Pathway")
KEGG_GD1 <- KEGG_GD[(KEGG_GD$KEGG_Pathway %in% KEGG_CD8$KEGG_Pathway), ]


head(KEGG_GD)
KEGG_GD$Oneminus <- 1 - KEGG_GD$FDR


p <- ggplot(KEGG_GD, aes(x = KEGG_Pathway, y = Oneminus)) +
  geom_bar(stat = "identity") + coord_flip()
p

kegg_sizes<- as.data.frame(lengths(kegg_human)) %>% rownames_to_column(., var = "kegg_pathway")
colnames(kegg_sizes) <- c("KEGG_Pathway", "Size")
head(kegg_sizes)
GD_K <- merge(KEGG_GD, kegg_sizes, by = "KEGG_Pathway")
GD_K$Num <- GD_K$NGenes/GD_K$Size


Comparison <- "GD"
KEGG_GD2 <- cbind(KEGG_GD, Comparison)
KEGG_GD2

idx <- ids2indices(kegg_human, id = rownames(v))
cam.EMRA.vs.Naive <- camera(v, idx, design, contrast = contr.matrix[, "EMRAvsNaive"])
head(cam.EMRA.vs.Naive, 5)
sig.EMRA.vs.Naive <- droplevels(subset(cam.EMRA.vs.Naive, FDR <= 0.01))
KEGG_CD8 <- rownames_to_column(sig.EMRA.vs.Naive, var = "KEGG_Pathway")
Comparison <- "CD8"
KEGG_CD8a <- cbind(KEGG_CD8, Comparison)


idx <- ids2indices(kegg_human, id = rownames(v))
cam.EMRA.vs.Naive <- camera(v, idx, design, contrast = contr.matrix[, "CD27LOvsNaive"])
head(cam.EMRA.vs.Naive, 5)
sig.EMRA.vs.Naive <- droplevels(subset(cam.EMRA.vs.Naive, FDR <= 0.01))
KEGG_CD8 <- rownames_to_column(sig.EMRA.vs.Naive, var = "KEGG_Pathway")
Comparison <- "CD8"
KEGG_CD8a <- cbind(KEGG_CD8, Comparison)



rbind(KEGG_GD2, KEGG_CD8a)

KEGG <- rbind(KEGG_CD8a, KEGG_GD2)

KEGG$FDR_Dir <- ifelse((KEGG$Direction == "Down"), KEGG$FDR * -1, KEGG$FDR * 1)

KEGG1 <- KEGG[, c("KEGG_Pathway", "Comparison", "FDR_Dir")]

writeCsvO(KEGG1)

KEGG2 <- spread(KEGG1, key = "KEGG_Pathway", value = "FDR_Dir")

KEGG3 <- KEGG2 %>% remove_rownames %>% column_to_rownames(var = "Comparison")

heatmap.2(t(as.matrix(KEGG3)), col = mycol, trace = "none", density.info = "none", scale = "row"
)
mycol <- colorpanel(1000,"blue","white","red")

mycol



writeCsvO(try)
head(try)

thesePath <- try$KEGG_Pathway








this <- as.matrix(kegg_human)
idx <- ids2indices(kegg_human, id = rownames(v))
cam.CD27LO.vs.EMRA <- camera(v, idx, design, contrast = contr.matrix[, "CD27LOvsEMRA"])
sig.CD27LO.vs.EMRA <- droplevels(subset(cam.CD27LO.vs.EMRA, FDR <= 0.01))


idx <- ids2indices(kegg_human, id = rownames(v))
cam.CD27HI.vs.Naive <- camera(v, idx, design, contrast = contr.matrix[, "CD27HIvsNaive"])
sig.CD27HI.vs.Naive <- droplevels(subset(cam.CD27HI.vs.Naive, FDR <= 0.01))

head(cam.CD27HI.vs.Naive, 5)

View(contr.matrix)
par(mfrow = c(1,1))
barcodeplot(efit$t[, 1], index = idx$`hsa04650 Natural killer cell mediated cytotoxicity`, main = "Natural killer cell mediated cytotoxicity", labels = c("VD1.CD27HI", "VD1.CD27LO"))




#### Less Interesting Genesets ####
## GO_genesets
load("~/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/RData_Objects/GO_genesets.rdata")

## Remove genesets over or equal to 300 length and below 10 genes
GO_sizes<- as.data.frame(lengths(Hs.c5)) %>% rownames_to_column(., var = "GO_pathway")
colnames(GO_sizes) <- c("GO_Pathway", "Size")
remove_these <- droplevels(subset(GO_sizes, Size >= 300 & Size <= 10))$GO_Pathway
isNameInIndex <- names(Hs.c5) %in% remove_these
GO_terms <- Hs.c5[!isNameInIndex]


idx <- ids2indices(GO_terms, id = rownames(v))
camera_results <- camera(v, idx, design, contrast = contr.matrix[, "CD27LOvsCD27HI"])
head(camera_results, 5)
sig.CD27LO.vs.CD27HI <- droplevels(subset(camera_results, FDR <= 0.001))
GO_GD <- rownames_to_column(sig.CD27LO.vs.CD27HI, var = "GO_Pathway")


GO_sizes <- as.data.frame(lengths(Hs.c5)) %>% rownames_to_column(., var = "GO_pathway")
colnames(GO_sizes) <- c("GO_Pathway", "Size")
head(GO_sizes)
GD_G <- merge(GO_GD, GO_sizes, by = "GO_Pathway")
head(GD_G)

dim(GD_G)

GD_G$Num <- GD_G$NGenes/GD_G$Size
head(GD_G)

GD_G

# ord_GD_GO <- arrange(GD_G, FDR)
# ord_GD_GO$min_log10_FDR <- -log10(ord_GD_GO$FDR)
# 
# head(ord_GD_GO)
# p <- ggplot(ord_GD_GO, aes(x = GO_Pathway, y = min_log10_FDR)) +
#   geom_bar(stat = "identity") 
# 
# p
# ?sec_axis
# Comparison <- "GD"
# KEGG_GD2 <- cbind(KEGG_GD, Comparison)
# KEGG_GD2

cam.EMRA.vs.Naive <- camera(v, idx, design, contrast = contr.matrix[, "EMRAvsNaive"])
head(cam.EMRA.vs.Naive, 5)
sig.EMRA.vs.Naive <- droplevels(subset(cam.EMRA.vs.Naive, FDR <= 0.001))
GO_CD8 <- rownames_to_column(sig.EMRA.vs.Naive, var = "KEGG_Pathway")
dim(GO_CD8)



Comparison <- "CD8"
KEGG_CD8a <- cbind(KEGG_CD8, Comparison)





head(GD_GENESETS)

common_elements <- combn(GD_GENESETS, 2, 
                         FUN = function(x) intersect(x[[1]], x[[2]]), simplify = F)

names(common_elements) <- vapply(combn(names(GD_GENESETS), 2, simplify = F), 
                                 paste, collapse = "___", FUN.VALUE = character(1))

head(common_elements)

common_num <- lengths(common_elements)

library(tidyverse)
try <- as.data.frame(common_num)
head(try)
try1 <- rownames_to_column(try, var = "geneset")

this <- try1 %>% separate("geneset", into = c("data1", "data2"), sep = "___")
head(this)

length_of_genesets <- data.frame(data1 = character(),
                                 length1 = double(),
                                 stringsAsFactors = F)
c <- 1
for(i in names(GD_GENESETS)){
  work <- GD_GENESETS[[i]]
  that <- length(work)
  length_of_genesets[c, "data1"] <- i
  length_of_genesets[c, "length1"] <- that
  c <- c + 1
}

head(length_of_genesets)
names(GD_GENESETS) %in% levels(as.factor(this$data1))

this1 <- merge(this, length_of_genesets, by = "data1", all.x = T)
subset(this1, data2 == "GO_RECEPTOR_INHIBITOR_ACTIVITY")


length_of_genesets <- data.frame(data2 = character(),
                                 length2 = double(),
                                 stringsAsFactors = F)
c <- 1
for(i in names(GD_GENESETS)){
  work <- GD_GENESETS[[i]]
  that <- length(work)
  length_of_genesets[c, "data2"] <- i
  length_of_genesets[c, "length2"] <- that
  c <- c + 1
}

this2 <- merge(this1, length_of_genesets, by = "data2", all.x = T)

head(this2)
Geneset_overlap <- this2

Geneset_overlap$specific_1 <- Geneset_overlap$length1 - Geneset_overlap$common_num
Geneset_overlap$specific_2 <- Geneset_overlap$length2 - Geneset_overlap$common_num

Geneset_overlap$neg_agreement <- 2350 - (Geneset_overlap$common_num + Geneset_overlap$specific_1 + Geneset_overlap$specific_2)
head(Geneset_overlap)

Geneset_overlap$kappa <- (Geneset_overlap$common_num + 0)/(Geneset_overlap$common_num + Geneset_overlap$specific_1 +Geneset_overlap$specific_2 + 0)
head(Geneset_overlap)

test <- droplevels(subset(Geneset_overlap, data1 == "GO_INTERFERON_GAMMA_MEDIATED_SIGNALING_PATHWAY")) #- NA are where there aren't compared
as.factor(test$data2)

test <- droplevels(subset(Geneset_overlap, data2 == "GO_INTERFERON_GAMMA_MEDIATED_SIGNALING_PATHWAY")) #- NA are where there aren't compared
as.factor(test$data1)

writeCsvO(Geneset_overlap)

kappa_scores_GD <- Geneset_overlap[, c("data1", "data2", "kappa")]

library(GGally)
library(network)
mm.net <- network(test[,1:2], directed = T)

ggnet2(mm.net,
       labelon = TRUE,
       size = 2, vjust = -0.6, mode = "kamadakawai", label.size = 3)


test <- droplevels(subset(Geneset_overlap, kappa != 0))

mm.net

install.packages("geomnet")
data(madmen, package = 'geomnet')

head(madmen)

?ggnetworkmap()


KS_GD <- spread(kappa_scores_GD, "data2", value = "kappa")
head(KS_GD)
str(KS_GD)

KS_GD1 <- column_to_rownames(KS_GD, var = "data1")
KS_GD1
is.na(KS_GD1) <- 0

getwd()
GO_GD <- read.csv("Bulk/Output/DEgenes/GO_GD.csv")

head(GO_GD)

GD <- as.character(GO_GD$GO_Geneset)

GD1 <- Geneset_overlap[Geneset_overlap$data1 %in% GD, ]
GD2 <- GD1[GD1$data2 %in% GD, ]

head(GD2)


## Making of a binary gene-term matrix
### Take only the terms that are significant in GD


GD_GENESETS <- Hs.c5[names(Hs.c5) %in% GD]
that <- as.data.frame(unlist(GD_GENESETS)) %>% rownames_to_column(var = "Geneset")
that$Geneset <- gsub("[0-9]{1,4}$", "", that$Geneset)
colnames(that) <- c("Geneset", "Genes")
length(that$Genes)

gene_term_binary <- as.data.frame.matrix(table(that))

head(gene_term_binary)





library(Matrix)
#https://stackoverflow.com/a/51421029/1412059
fun <- function(x) {
  n <- 0.5 + sqrt(0.25 + 2 * length(x)) #calculate dimension
  i <- sequence(seq_len(n - 1)) #row indices
  j <- rep(seq_len(n - 1), seq_len(n - 1)) + 1 # column indices
  sparseMatrix(i = i, j = j, x = x, triangular = TRUE, dims = c(n, n))
}

output <- fun(common_num)
diag(output) <- lengths(GD_GENESETS)
dimnames(output) <- rep(list(names(GD_GENESETS)), 2)

View(output)

head(output)

output1 <- as.data.frame(output)




idx <- ids2indices(Hs.c5, id = rownames(v))
cam.CD27LO.vs.CD27HI <- camera(v, idx, design, contrast = contr.matrix[, 1])
head(cam.CD27LO.vs.CD27HI, 5)

sig.CD27LO.vs.CD27HI <- droplevels(subset(cam.CD27LO.vs.CD27HI, FDR <= 0.01))
GO_GD <- rownames_to_column(sig.CD27LO.vs.CD27HI, var = "GO_Geneset")
writeCsvO(GO_GD)

barcodeplot(efit$t[, 1], index = idx$GO_CELL_KILLING, main = "CD27LO vs CD27HI")

cam.EMRA.vs.Naive <- camera(v, idx, design, contrast = contr.matrix[, "EMRAvsNaive"])
head(cam.CD27LO.vs.CD27HI, 5)

sig.EMRA.vs.Naive <- droplevels(subset(cam.EMRA.vs.Naive, FDR <= 0.01))
GO_CD8 <- rownames_to_column(sig.EMRA.vs.Naive, var = "GO_Geneset")
writeCsvO(GO_CD8)



nrow(GO_GD)
GO_GD$Directed_FDR <- ifelse((GO_GD$Direction == "Up"), GO_GD$FDR * 1, GO_GD$FDR * -1)

head(GO_GD)


ggplot(GO_GD, aes(x = GO_Geneset, y = Directed_FDR)) + geom_point()

# Immunological Signatures
load("~/OneDrive/UoB/PhD/Projects/GammaDelta/Bulk/RData_Objects/Immunological_Signatures.rdata")
idx <- ids2indices(Hs.c7, id = rownames(v))
cam.CD27LO.vs.CD27HI <- camera(v, idx, design, contrast = contr.matrix[, 1])
head(cam.CD27LO.vs.CD27HI, 10)

barcodeplot(efit$t[, 1], index = idx$GSE26495_NAIVE_VS_PD1LOW_CD8_TCELL_UP , 
            index2 = idx$GSE26495_NAIVE_VS_PD1LOW_CD8_TCELL_DN, main = "NaiveVsCD8")










install.packages("DOSE")
BiocManager::install("PathwaySplice")
library(DOSE)
library(PathwaySplice)

gene.based.table <- makeGeneTable(Hs.c5)
?makeGeneTable
View(featureBasedData)

res <- runPathwaySplice(gene.based.table,genome='hg19',
                        id='ensGene',test.cats=c('GO:BP'),
                        go.size.limit=c(5,30),method='Wallenius')

# labeling each node by gene set name
enmap <- enrichmentMap(res,n=10,similarity.threshold=0.3,
                       label.node.by.index = FALSE)

# labeling each node by gene set index
enmap <- enrichmentMap(res,n=10,similarity.threshold=0.3,
                       label.node.by.index = TRUE)

## Not run: 
# illustrates specification of output file directory
# Enable interactive map and label each node by gene set index
enmap <- enrichmentMap(res,n=10,fixed=FALSE, similarity.threshold=0.3,
                       label.node.by.index = TRUE, output.file.dir=tempdir())

enmap <- enrichmentMap(res,n=10,similarity.threshold=0.3,
                       label.node.by.index = FALSE, output.file.dir=tempdir())
