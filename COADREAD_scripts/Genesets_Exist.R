# A script to ensure that all genes in the Genesets I've chosen are in the Counts_cqn matrix
library(UsefulFunctions)
library(tidyverse)
library(ggpubr)

# Load data
load("./R_Data/Counts_clean.RData")

# The Genesets
## Cell type Genesets
CTGenesets <- read.csv("./Data/Genesets/Raw/Pre_Change_Investigative_Genesets.csv", stringsAsFactors = T)
CTGenesets$HUGO.symbols <- as.character(CTGenesets$HUGO.symbols)
CTGenesets[!'%in%'(CTGenesets$HUGO.symbols, rownames(Counts_cqn)), ]

rownames(Counts_cqn)[grep("CCN4", rownames(Counts_cqn))]

### Castro
CTGenesets$HUGO.symbols[CTGenesets$HUGO.symbols == "NLF1"] <- "C2CD4A"
CTGenesets$HUGO.symbols[CTGenesets$HUGO.symbols == "NLF2"] <- "C2CD4B"
CTGenesets$HUGO.symbols[CTGenesets$HUGO.symbols == "C6orf145"] <- "PXDC1"
CTGenesets$HUGO.symbols[CTGenesets$HUGO.symbols == "C1orf113"] <- "SH3D21"
CTGenesets$HUGO.symbols[CTGenesets$HUGO.symbols == "IL1"] <- "IL1B"

### CAF_Alfie
CTGenesets$HUGO.symbols[CTGenesets$HUGO.symbols == "GPR124"] <- "ADGRA2"
CTGenesets$HUGO.symbols[CTGenesets$HUGO.symbols == "LPPR4"] <- "PLPPR4"
CTGenesets$HUGO.symbols[CTGenesets$HUGO.symbols == "ODZ4"] <- "TENM4"
CTGenesets$HUGO.symbols[CTGenesets$HUGO.symbols == "TXNDC3"] <- "NME8"
CTGenesets$HUGO.symbols[CTGenesets$HUGO.symbols == "WISP1"] <- "CCN4"

write.table("./Data/Genesets/Investigative_Genesets.txt", x = CTGenesets, row.names = F, quote = F, sep = "\t") # Change HUGO.symbols 03/06/2009 to 
CTGenesets <- read.delim("./Data/Genesets/Investigative_Genesets.txt", stringsAsFactors = T)

CTGenesets$HUGO.symbols <- as.character(CTGenesets$HUGO.symbols)
CTGenesets[!'%in%'(CTGenesets$HUGO.symbols, rownames(Counts_cqn)), ]


# Signature Genesets
SigGen <- read.csv("./Data/Genesets/Raw/Pre_Change_Signature_Genesets.csv", stringsAsFactors = T)
SigGen$HUGO.symbols <- as.character(SigGen$HUGO.symbols)
SigGen[!'%in%'(SigGen$HUGO.symbols, rownames(Counts_cqn)), ]



rownames(Counts_cqn)[grep("ICAM", rownames(Counts_cqn))]

### HTG
SigGen$HUGO.symbols[SigGen$HUGO.symbols == "ATP5F1"] <- "ATP5PB"

write.table("./Data/Genesets/Signature_Genesets.txt", x = SigGen, row.names = F, quote = F, sep = "\t") # Change HUGO.symbols 03/06/2009 to 
SigGen <- read.delim("./Data/Genesets/Signature_Genesets.txt", stringsAsFactors = T)

SigGen$HUGO.symbols <- as.character(SigGen$HUGO.symbols)
SigGen[!'%in%'(SigGen$HUGO.symbols, rownames(Counts_cqn)), ]

# GO Terms



