# A script to sort out the Bindea immunome dataset for Cytoscape
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


# Get all the genes that I can...
Counts_clean <- cqn_Counts$counts
ensembl_DB <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
Gene_Map <- getBM(attributes = c("entrezgene_id", "hgnc_symbol"),
                  filters = "hgnc_symbol", values = rownames(Counts_clean), ensembl_DB)
colnames(Gene_Map)[colnames(Gene_Map) == "entrezgene_id"] <- "EntrezGene"


#### Immunome Dataset - Bindea et al., 2013 ####
immunome <- read.csv("./Data/Genesets/immunome.csv", stringsAsFactors = F)
immunome$EntrezGene <- trim(immunome$EntrezGene)
immunome$Symbol <- trim(immunome$Symbol)
immunome <- immunome[!duplicated(immunome), ]

# Gene_Map <- Gene_Map[!is.na(Gene_Map$EntrezGene), ]
# Find the genes that were lost
lost_data <- immunome[!'%in%'(immunome$Symbol, Gene_Map$hgnc_symbol), ]
write.csv("./Output/Geneset_vetting/lost_data.csv", x = lost_data, row.names = F)

#### Probe for alternative genes ####
# Gene_Map[grep("RNU2", Gene_Map$hgnc_symbol, ignore.case = T), ]
lost_data_alt <- read.csv("./Output/Geneset_vetting/lost_data_mismatch_combine.csv")[, c()]

# Replace data in immunome with the alternative gene aliases
immunome1 <- immunome[, c("CellType", "Symbol")]
immunome1$Symbol <- ifelse((immunome1$Symbol %in% lost_data_alt$Symbol), as.character(lost_data_alt$Alternative), as.character(immunome1$Symbol))

## Remove NAs
immunome2 <- immunome1[!is.na(immunome1$Symbol), ]

## Lose 13 genes

# Remove duplicates
immunome3 <- immunome2[!duplicated(immunome2), ]
write.csv("./Output/Genesets/Immunome_vetted.csv", x = immunome3, row.names = F)

# Make it a part of a big dataset
immunome3$Source <- "Immunity, 2013"
immunome3$Type_of_data <- "CellType"

## Th17 geneset
CT <- read.csv("./Data/Genesets/Castro_collated.csv")
need_alt <- CT[!('%in%'(CT$Symbol, Gene_Map$hgnc_symbol)), ] %>% droplevels()

write.csv("./Output/Geneset_vetting/Castro_collated_required.csv", row.names = F, x = need_alt)

# NLF1 = C2CD4A
# NLF2 = C2CD4B
# C6orf145 = PXDC1
# C1orf113 = SH3D21

# Gene_Map[grep("IL1B", Gene_Map$hgnc_symbol), ]
Cast_v <- read.csv("./Output/Geneset_vetting/Castro_collated_vetted.csv")

# Replace data in castro with the alternative gene aliases
CT1 <- CT[, c("CellType", "Symbol")]
CT1$Symbol <- ifelse((CT1$Symbol %in% Cast_v$Symbol), as.character(Cast_v$Alternative), as.character(CT1$Symbol))

## Remove NAs
CT2 <- CT1[!is.na(CT1$Symbol), ]

# 0 lost

CT2$Source <- "PLoS, 2017"
CT2$Type_of_data <- "CellType"












#  Bind all data together
Genesets_of_interest <- rbind(immunome3, CT2)


write.csv("./Output/Genesets/Genesets_of_interest.csv", x = Genesets_of_interest, row.names = F)
