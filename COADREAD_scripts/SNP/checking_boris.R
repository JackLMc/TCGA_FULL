getwd()

TB <- read.csv("./Output/file_id_pats_SNP.csv")[, c("file_id", "Subtype")]
Boris <- read.delim("./Data/Boris_clusters.txt")
Boris$file_id <- gsub("-", ".", Boris$file_id)


colnames(Boris) <- c("file_id", "Subtype_B")
dim(Boris)
dim(TB)

see <- merge(Boris, TB, by = "file_id")
all.equal(see$Subtype_B, see$Subtype)




### Are the SNPs the same in the cancer?





