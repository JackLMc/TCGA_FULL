# A R script to investigate associations of SNPs with Patient groups
library(tidyverse)
library(UsefulFunctions)

# Looking at the SNPs
MSS_hiCIRC <- read.delim("./Output/associations.annotated_13_02.txt") # A = MSS, B = MSS_hiCIRC

# Split into Autosomes and Mitochondria 
Autosomes <- MSS_hiCIRC[!('%in%'(MSS_hiCIRC$Chromosome, c("X", "Y", "MT"))),] %>% droplevels() # Look at the allele test, phred score around 60 = 10E-6 Pval
# write.csv(x = Autosomes, file = "./Output/Boris_Autosome_Associations_MSS_MSShiCIRC.csv", row.names = F)

Mito <- MSS_hiCIRC[('%in%'(MSS_hiCIRC$Chromosome, c("MT"))),] %>% droplevels() # Look at the genotype test, phred score around 25?
# write.csv(x = Mito, file = "./Output/Boris_Mito_Associations_MSS_MSShiCIRC.csv", row.names = F)




ordered_Mito <- Mito[order(Mito$genotype_rest_phred, decreasing = T),]
head(ordered_Mito)

ordered_Auto <- Autosomes[order(Autosomes$allele_test_phred, decreasing = T),]
head(ordered_Auto)

head(ordered_Mito)


#### END ####