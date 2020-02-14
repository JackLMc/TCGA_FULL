#!/usr/bin/env Rscript
# Run this script as
# Rscript --vanilla ProcessGenotypeAssociation.R clusters.txt matrix_genotypes.final.txt > association.txt

library("stringr")

# Change it if needed, can use multiple labels in one cluster
clustersA=c("MSS","MSS-hiCIRC")
clustersB=c("MSI-H")

args = commandArgs(trailingOnly=TRUE)

# First argument is the file with cluster annotations
cluster_file=args[[1]]
#cluster_file="clusters.txt"

# Second argument is the genotype file 
genotype_file=args[[2]]
#genotype_file="tmp.txt"

sets=read.delim(cluster_file)

# IDs of elements in different clusters
idA=as.character(sets[sets[,2] %in% clustersA,1,drop=T])
idB=as.character(sets[sets[,2] %in% clustersB,1,drop=T])


inputpipe=file(genotype_file,"rt",blocking=T)

# Read the first line, it should contain the header
line1 <- str_split(readLines(inputpipe,n=1),"\t")[[1]]
nc=length(line1)

# Create the index for columns in the 2 clusters:
idxA=(line1 %in% idA)
#table(idxA)

idxB=(line1 %in% idB)
#table(idxB)

# Converting p-value to phred score
phred=function(x){if(x==0){10000} else {-10*log10(x)}}

# Function to calculate the p-values for one SNP (one row)
AssociationTest=function(line){
  # Genotype counts
  lineA=line[idxA]
  cA=sapply(c("0","1","2"),function(x){sum(lineA==x)})
  
  lineB=line[idxB]
  cB=sapply(c("0","1","2"),function(x){sum(lineB==x)})
  
  # Matrix of allele counts, 2x2
  m2=rbind(c(2*cA[[1]]+cA[[2]],2*cA[[3]]+cA[[2]]),
          c(2*cB[[1]]+cB[[2]],2*cB[[3]]+cB[[2]]))
  # Matrix of genotype counts, 2x3
  m3=rbind(cA,cB)
  
  res=NULL
  res$phred_allele=phred(fisher.test(m2)$p.value)
  res$phred_genotype=phred(fisher.test(m3)$p.value)
  res$genotypes=m3

  return(res)
}

# Write the header of the table
header=c(line1[[1]],"A00","A01","A11","B00","B01","B11","allele_test_phred","genotype_rest_phred")
writeLines(paste0(header,collapse = "\t"))

#system.time(
while(length(line <- readLines(inputpipe,n=1)) > 0) {
  line=str_split_fixed(line,"\t",nc)[1,]
  res=AssociationTest(line)
  output=c(line[[1]],c(t(res$genotypes)),round(c(res$phred_allele,res$phred_genotype),digits = 2))
  writeLines(paste0(output,collapse = "\t"))
}
#)

close(inputpipe)
