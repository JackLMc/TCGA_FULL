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

# Create binary file of the genesets
library(tidyverse)
Binary_BP <- BioProc %>% 
  enframe %>% 
  unnest %>% 
  count(name, value) %>% 
  spread(value, n, fill = 0) %>% column_to_rownames(., var = "name") %>%
  as.matrix()

dimnames(Binary_BP) <- NULL

## Execute SGSE using Fisher-transformed correlation coefficients as 
## the gene-level statistics, the standardized mean difference as the 
## gene set statistic and a correlation adjusted two-sided, 
## two-sample t-test for the determination of statistical significance,
## all PCs with non-zero eigenvalues for spectral enrichment and 
## variance weights 
sgse.results = sgse(data = t(Counts_cqn), 
                    prcomp.output = prcomp.output, 
                    gene.sets = Binary_BP,
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
# 
# save.image("./Data/SGSE.RData")
load("./Data/SGSE.RData")

# ## Execute SGSE again but using RMT scaled variance weights
# sgse.results = sgse(data = t(Counts_cqn), 
#                     prcomp.output = prcomp.output, 
#                     gene.sets = Binary_BP,
#                     gene.statistic = "z", 
#                     transformation = "none",
#                     gene.set.statistic = "mean.diff",
#                     gene.set.test = "cor.adj.parametric",
#                     pc.selection.method = "all",
#                     pcgse.weight = "rmt.scaled.var")
# 
# ## Display the SGSE weights for the first 5 PCs 
# sgse.results$weights[1:5]
# 
# ## Display the SGSE p-values for the first 5 gene sets 
# sgse.results$sgse[1:5]                       

## Execute SGSE again using RMT scaled variance weights and  
## all RMT-significant PCs at alpha=.05
sgse.results = sgse(data = t(Counts_cqn), 
                    prcomp.output = prcomp.output, 
                    gene.sets = Binary_BP,
                    gene.statistic = "z", 
                    transformation = "none",
                    gene.set.statistic = "mean.diff",
                    gene.set.test = "cor.adj.parametric",
                    pc.selection.method = "rmt",
                    rmt.alpha = .05,
                    pcgse.weight = "rmt.scaled.var")


# ## Display the indexes of the RMT-significant PCs
# sgse.results$pc.indexes
# 
# ## Display the SGSE p-values for the first 5 gene sets 
# sgse.results$sgse[1:5]                             

# Get the most significant genesets
filtered_genesets <- names(BioProc[sgse.results$sgse <= 0.01]) # highly significant ones
write.table(filtered_genesets, file = "./Data/Genesets/Filtered_set_PCGSE.txt", row.names = F, quote = F)
filtered_genesets[grepl("MITO", filtered_genesets)]









