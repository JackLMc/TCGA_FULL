# A script to reduce the matrix genotypes file to only those SNPs that are in Toju's analysis
## Read in required data
SNP <- read.delim("./Data/SNPs/matrix_genotypes.final.txt")

SNP1 <- column_to_rownames(SNP, var = "snpID") %>%
  as.matrix() %>%
  t()
rownames(SNP1) <- gsub("^X", "", rownames(SNP1))



SNPs_to_read <- read.csv("./Data/SNPs/Toju/SNP_Toju.csv")$SNP # these are rs numbers... Need to do this on the final matrix...


SNP2 <- SNP1[, colnames(SNP1) %in% SNPs_to_read]



