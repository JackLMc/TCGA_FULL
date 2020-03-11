# A script to investigate protein expression between patient subgroups
library(UsefulFunctions)
library(tidyverse)
library(ggpubr)

cbcols <- c("MSS-hiCIRC" = "#999999",
            "MSI-H" = "#56B4E9",
            "MSS" = "#009E73")

my_comparisons <- list(c("MSS-hiCIRC", "MSI-H"),
                       c("MSS-hiCIRC", "MSS"),
                       c("MSI-H", "MSS"))



RPPA <- read.delim("./Data/RPPA/COADREAD.rppa.txt")
RPPA1 <- separate(RPPA, "Composite.Element.REF", c("Gene.Name", "Composite.Element.REF"), sep = "\\|")

# Look at the expression of various components across subtype
RPPA2 <- RPPA1 %>% gather(key = "Patient.ID", value = "Protein_expression", -Gene.Name, -Composite.Element.REF)
RPPA2$Patient.ID <- samptopat(RPPA2$Patient.ID)

## Read in clinical data
pat_sub <- read.csv("./Output/Patient_Subtypes_09_03.csv")
RPPA3 <- merge(RPPA2, pat_sub, by = "Patient.ID")
RPPA3 <- factorthese(RPPA3, c("Gene.Name", "Composite.Element.REF"))

levels(RPPA3$Gene.Name)

GoI <- droplevels(subset(RPPA3, Gene.Name == "TP53"))
nlevels(GoI$Composite.Element.REF) == 1 # Check there is only one probe

ggplot(GoI, aes(x = Subtype, y = Protein_expression)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(aes(Subtype, fill = Subtype),
              scale = "width", alpha = 0.8) +
  scale_fill_manual(values = cbcols) +
  labs(x = "Subtype", y = "Protein Expression") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top") + 
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test")


## Look at stats for all genes - will need to average across probe
data <- data.frame()
c <- 1
for(i in levels(RPPA3$Gene.Name)){
  print(i)
  work <- droplevels(subset(RPPA3, Gene.Name == i))
  data[c, "Gene"] <- i
  data[c, "Number.Probes"] <- nlevels(work$Composite.Element.REF)
  c <- c + 1
}


head(data)

data[data$Number.Probes != 1, ]

