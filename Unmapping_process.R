library(UsefulFunctions)
library(tidyverse)
library(ggpubr)

my_comparisons <- list(c("MSS-hiCIRC", "MSI-H"),
                       c("MSS-hiCIRC", "MSS"),
                       c("MSI-H", "MSS"))
cbcols <- c("MSS-hiCIRC" = "#999999",
            "MSI-H" = "#56B4E9",
            "MSS" = "#009E73")


stats <- read.delim("./Data/statistics_unmapped.txt")
mapping <- read.delim("./Data/GDC_large_mapping_TCGA.txt")


stats <- separate(stats, "sample", c("sample", "extra"), sep = ":")
stats$extra <- as.factor(stats$extra)
stats$SAMPLES <- gsub(".idxstats", "", stats$sample)
stats_mapped <- mapping[mapping$file_name %in% stats$SAMPLES, ]
mapped_no_norm <- droplevels(subset(stats_mapped, cases.0.samples.0.sample_type != "Solid Tissue Normal" & cases.0.samples.0.sample_type != "Blood Derived Normal"))
mapped_no_norm$Patient.ID <- gsub("-", ".", mapped_no_norm$cases.0.submitter_id)

# Combine with patient subtypes
pat_sub <- read.csv("./Output/Patient_Subtypes.csv")
TCGA_pats <- merge(mapped_no_norm, pat_sub, by = "Patient.ID")

TCGA_pats <- TCGA_pats[, c("Patient.ID", "file_name", "Subtype", "CIRC_Genes")]
stats$file_name <- stats$SAMPLES

merged_data <- merge(stats[,c("file_name", "unmapped")], TCGA_pats, by = "file_name")

merged_data_clean <- merged_data[ c("Patient.ID", "unmapped", "Subtype")]


merged_data_clean$Patient.ID <- as.factor(merged_data_clean$Patient.ID)
DF1 <- data.frame(stringsAsFactors = F)
c <- 1
for(i in levels(merged_data_clean$Patient.ID)){
  work <- droplevels(subset(merged_data_clean, Patient.ID == i))
  total <- sum(work$unmapped)
  DF1[c, "Patient.ID"] <- i
  DF1[c, "unmapped"] <- total
  DF1[c, "Subtype"] <- as.character(levels(work$Subtype))
  c <- c + 1
}

DF1$logged <- log10(DF1$unmapped + 1)
ggplot(DF1, aes(x = Subtype, y = logged)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(aes(Subtype, fill = Subtype),
              scale = "width", alpha = 0.8) +
  scale_fill_manual(values = cbcols) +
  labs(x = "MSI Status", y = "Log10(unmapped reads + 1)") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top") + 
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test")

range(DF1$unmapped)
range(droplevels(subset(DF1, Subtype == "MSS-hiCIRC"))$unmapped)
range(droplevels(subset(DF1, Subtype == "MSS"))$unmapped)
range(droplevels(subset(DF1, Subtype == "MSI-H"))$unmapped)



