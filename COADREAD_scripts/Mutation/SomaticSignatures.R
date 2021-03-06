# A script to investigate the existence of mutational signatures within the TCGA COADREAD cohort
# Packages
library(UsefulFunctions)
library(tidyverse)
library(ggpubr)
large <- c("SomaticSignatures",
           "BSgenome.Hsapiens.1000genomes.hs37d5")
for (lib in large){
  if (!require(lib, character.only = T)){
    BiocManager:: install(lib, dependencies = T)
    suppressMessages(library(lib, character.only = T, quietly = T))}}

my_comparisons <- list(c("MSS-hiCIRC", "MSI-H"),
                       c("MSS-hiCIRC", "MSS"),
                       c("MSI-H", "MSS"))
cbcols <- c("MSS-hiCIRC" = "#999999",
            "MSI-H" = "#56B4E9",
            "MSS" = "#009E73")

load("./R_Data/Mutation_clean.RData")
rm("mutation_maf") # Don't need this for this analysis, reduce burden by removing

pat_sub <- read.csv("Output/Patient_Subtypes_09_03.csv")[, c("Patient.ID", "Subtype")]

# Preprocess
##Only take SNPs
somsig <- droplevels(subset(tcga_mut, Variant_Type == "SNP"))

all.equal(somsig$Start_Position, somsig$End_Position)
## Ensure this causes these columns to be equal
all.equal(somsig$Reference_Allele, somsig$Match_Norm_Seq_Allele1) 
all.equal(somsig$Reference_Allele, somsig$Match_Norm_Seq_Allele2) 
all.equal(somsig$Match_Norm_Seq_Allele2, somsig$Match_Norm_Seq_Allele1) 

## Label the base changes
somsig <- somsig %>%
  mutate(BaseChange = ifelse((as.character(Reference_Allele)==as.character(Tumor_Seq_Allele1) & as.character(Reference_Allele) == as.character(Tumor_Seq_Allele2)), "NoCh",
                             ifelse(as.character(Reference_Allele)==as.character(Tumor_Seq_Allele1), paste(Reference_Allele,Tumor_Seq_Allele2, sep = "."), paste(Reference_Allele, Tumor_Seq_Allele1, sep = "."))))


# somsig$BaseChange <- apply(somsig[names(somsig) %in% c("Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2")], 1, 
#                            function(i) paste0(unique(i), collapse = '.')) 
# somsig$BaseChange <- replace(somsig$BaseChange, nchar(somsig$BaseChange) == 1, "NoCh")
somsig$BaseChange <- gsub(".Absent", "", somsig$BaseChange)

## Change the BaseChange to a character
somsig$BaseChange <- as.character(somsig$BaseChange)
somsig$MutatedBase <- as.factor(lastvalue(somsig$BaseChange, 1))
somsig$BaseChange <- as.factor(somsig$BaseChange)

## Bind to the subtypes
sca_data <- merge(somsig, pat_sub[, c("Patient.ID", "Subtype")], by = "Patient.ID")

# Start the SomaticSignatures pipeline
## Create the range by making an IRange object with the start and end positions
range <- IRanges(start = sca_data$Start_Position, end = sca_data$End_Position)

## Create a VRange object consisting of the Chromosome, Reference and Mutated Bases
### Use the NCBI_build for the study
sca_vr <- VRanges(seqnames = sca_data$Chromosome, ranges = range, ref = sca_data$Reference_Allele,
                  alt = sca_data$MutatedBase, sampleNames = sca_data$Tumor_Sample_Barcode, subtype = sca_data$Subtype,
                  study = sca_data$NCBI_Build, PatID = sca_data$Patient.ID)

## Count number of Somatic variants per study
sort(table(sca_vr$subtype), decreasing = T)

## Find the context of the mutations with respect to a library genome
sca_motifs <- mutationContext(sca_vr, BSgenome.Hsapiens.1000genomes.hs37d5)

# Somatic spectrum across patient types
pdf("./Figures/Mutation_Spectrum/Somatic_spec_Subtype.pdf", width = 6, height = 6)
plotMutationSpectrum(sca_motifs, "subtype") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = cbcols)
dev.off()

# Construct matrix M {motifs × studies}. Normalize to return frequency 
sca_mm <- motifMatrix(sca_motifs, group = "PatID", normalize = T) # Look for PatID
head(round(sca_mm, 4))

# Remove any NA columns
sca_mm[sca_mm == "NaN"] <- NA
sca_mm1 <- as.data.frame(sca_mm)
sca_mm2 <- sca_mm1[ , apply(sca_mm1, 2, function(x) !any(is.na(x)))]

## Back to matrix
sca_mm3 <- as.matrix(sca_mm2)

## Use to estimate the number of signatures 
set.seed(123)
n_sigs = 2:7

gof_nmf <- assessNumberSignatures(sca_mm3, n_sigs, nReplicates = 10)

pdf("./Figures/Mutation_Spectrum/NMF_Sigs.pdf", width = 6, height = 6)
plotNumberSignatures(gof_nmf)
dev.off()

# gof_pca <- assessNumberSignatures(sca_mm3, n_sigs, pcaDecomposition) # these are pants
# plotNumberSignatures(gof_pca)

# n_sigs should be based on gof
# Somatic signatures identified using NMF and PCA
n_sigs <- 4
sigs_nmf <- identifySignatures(sca_mm3, n_sigs, nmfDecomposition)
# sigs_pca <- identifySignatures(sca_mm3, n_sigs, pcaDecomposition)

sigs_nmf@signatures
sigs_nmf@samples
sigs_nmf@fitted

w <- signatures(sigs_nmf)
w_norm <- t(t(w) / colSums(w))
sum(as.numeric(w_norm[,1]))

# Plot the NMF
pdf("./Figures/Mutation_Spectrum/NMF_Sigs_Contrib.pdf", width = 8, height = 8)
SomaticSignatures:: plotSignatures(sigs_nmf) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

# plotObservedSpectrum(sigs_nmf) # too many colours
# plotFittedSpectrum(sigs_nmf) # too many colours
plotSignatureMap(sigs_nmf)

# Order by most contributing
library(reshape2)
w_df <- melt(signatures(sigs_nmf))
w_df1 <- spread(w_df, key = "Var2", value = "value")
w_df1[order(w_df1$S3, decreasing = T)[1:10],]

#### Cosine Similarity ####
# install.packages("lsa", dependencies = T)
library(lsa)
data("signatures21")
resultS1 <- as.data.frame(cosine(sigs_nmf@signatures[,1], signatures21)) %>% rownames_to_column(var = "Signature")
resultS2 <- as.data.frame(cosine(sigs_nmf@signatures[,2], signatures21)) %>% rownames_to_column(var = "Signature")
resultS3 <- as.data.frame(cosine(sigs_nmf@signatures[,3], signatures21)) %>% rownames_to_column(var = "Signature")
resultS4 <- as.data.frame(cosine(sigs_nmf@signatures[,4], signatures21)) %>% rownames_to_column(var = "Signature")

result <- Reduce(merge, list(resultS1, resultS2, resultS3, resultS4))

colnames(result) <- c("Alex_Signature", "S1", "S2", "S3", "S4")
write.csv("./Output/Cosine_results.csv", row.names = F, x = result)

# Compare signatures across patients
try <- as.data.frame(sigs_nmf@samples) %>% rownames_to_column(var = "Patient.ID")
try1 <- merge(try, pat_sub[, c("Patient.ID", "Subtype")], by = "Patient.ID")

Sig_Freq <- as.data.frame(try1 %>%
                            gather(key = "Signature", value = "Frequency", -Patient.ID, -Subtype) %>%
                            dplyr:: group_by(Subtype, Signature) %>%
                            dplyr:: summarise(mean(Frequency)) %>%
                            spread(key = "Signature", value = "mean(Frequency)"))

library(gplots)
mycol <- colorpanel(1000,"blue","white","red")
SF <- data.frame(Sig_Freq[, names(Sig_Freq) != "Subtype"], row.names = Sig_Freq[, names(Sig_Freq) == "Subtype"])

pdf("./Figures/Mutation_Spectrum/Sigs_Heatmap.pdf", width = 80, height = 61)
heatmap.2(as.matrix(SF), scale = "row",
          labRow = rownames(SF), labCol = colnames(SF), 
          col = mycol,
          trace = "none", density.info = "none", 
          margin = c(8,12), lhei = c(2,10),
          main = "")
dev.off()


### Unimportant
# Sig_Freq1 <- as.data.frame(try1 %>% gather(key = "Signature", value = "Frequency", -Patient.ID, -Subtype))
# 
# 
# S3 <- droplevels(subset(Sig_Freq1, Signature == "S3"))
# S4 <- droplevels(subset(Sig_Freq1, Signature == "S4"))
# S1 <- droplevels(subset(Sig_Freq1, Signature == "S1"))
# S2 <- droplevels(subset(Sig_Freq1, Signature == "S2"))
# 
# ggplot(S3, aes(x = Subtype, y = Frequency)) +
#   geom_boxplot(alpha = 0.5, width = 0.2) + 
#   geom_violin(aes(Subtype, fill = Subtype), scale = "width", alpha = 0.8) +
#   scale_fill_manual(values = cbcols) +
#   labs(x = "Subtype", y = "Rank transformed total number of mutations") +
#   theme_bw() +
#   theme(axis.text = element_text(size = 16)) +
#   theme(legend.direction = "horizontal", legend.position = "top") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   stat_compare_means(comparisons = my_comparisons, label = "p.signif")
