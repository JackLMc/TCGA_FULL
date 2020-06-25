# A script to filter the GO genesets
library(tidyverse)
library(UsefulFunctions)


# BiocManager::install("qusage")
library(qusage)
# BiocManager::install("PCGSE")
library(PCGSE)


load("./R_Data/Counts_clean.RData")
BioProc <- read.gmt("./Data/Genesets/GSEA/Symbol/c5.bp.v7.1.symbols.gmt")
prcomp.output <- prcomp(t(Counts_cqn), center = T, scale = T)

# try <- read.delim("./Data/Genesets/Investigative_Genesets.txt")
# head(try)

# attempt <- table(try) %>% as.matrix()
# head(attempt)

head(BioProc)
# Create binary file of the genesets
library(tidyverse)
attempt <- BioProc %>% 
  enframe %>% 
  unnest %>% 
  count(name, value) %>% 
  spread(value, n, fill = 0) %>% column_to_rownames(., var = "name") %>%
  as.matrix()

dimnames(attempt) <- NULL

## Execute SGSE using Fisher-transformed correlation coefficients as 
## the gene-level statistics, the standardized mean difference as the 
## gene set statistic and a correlation adjusted two-sided, 
## two-sample t-test for the determination of statistical significance,
## all PCs with non-zero eigenvalues for spectral enrichment and 
## variance weights 
sgse.results = sgse(data = t(Counts_cqn), 
                    prcomp.output = prcomp.output, 
                    gene.sets = attempt,
                    gene.statistic = "z", 
                    transformation = "none",
                    gene.set.statistic = "mean.diff",
                    gene.set.test = "cor.adj.parametric",
                    pc.selection.method = "all",
                    pcgse.weight = "variance")

## Display the PCGSE p-values for the first 5 gene sets for PC 1 
sgse.results$pcgse$p.values[1:5,1]

## Display the SGSE weights for the first 5 PCs 
sgse.results$weights[1:5]   

## Display the SGSE p-values for the first 5 gene sets 
sgse.results$sgse[1:5]

## Execute SGSE again but using RMT scaled variance weights
sgse.results = sgse(data = cqn_Counts, 
                    prcomp.output = prcomp.output, 
                    gene.sets = attempt,
                    gene.statistic="z", 
                    transformation="none",
                    gene.set.statistic="mean.diff",
                    gene.set.test="cor.adj.parametric",
                    pc.selection.method="all",
                    pcgse.weight="rmt.scaled.var")

## Display the SGSE weights for the first 5 PCs 
sgse.results$weights[1:5]   

## Display the SGSE p-values for the first 5 gene sets 
sgse.results$sgse[1:5]                          

## Execute SGSE again using RMT scaled variance weights and  
## all RMT-significant PCs at alpha=.05
sgse.results = sgse(data=cqn_Counts, 
                    prcomp.output=prcomp.output, 
                    gene.sets=attempt,
                    gene.statistic="z", 
                    transformation="none",
                    gene.set.statistic="mean.diff",
                    gene.set.test="cor.adj.parametric",
                    pc.selection.method="rmt",
                    rmt.alpha=.05,
                    pcgse.weight="rmt.scaled.var")

## Display the indexes of the RMT-significant PCs
sgse.results$pc.indexes

## Display the SGSE p-values for the first 5 gene sets 
sgse.results$sgse[1:5]                             



