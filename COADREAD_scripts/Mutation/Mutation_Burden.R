# A script to look at the mutation totals between the groups in the COADREAD TCGA project
## Packages
library(UsefulFunctions)
library(tidyverse)
library(ggpubr)

my_comparisons <- list(c("MSS-hiCIRC", "MSI-H"),
                       c("MSS-hiCIRC", "MSS"),
                       c("MSI-H", "MSS"))
cbcols <- c("MSS-hiCIRC" = "#999999",
            "MSI-H" = "#56B4E9",
            "MSS" = "#009E73")

load("./R_Data/Mutation_clean.RData")

pat_sub <- read.csv("Output/Patient_Subtypes_13_02.csv")[, c("Patient.ID", "Subtype")]

# Number of mutations between groups ----
# Count up
## Including silent
# find nonsynonymous mutations??
Mutation_numbers <- tcga_mut %>%
  dplyr:: group_by(Patient.ID, Variant_Type) %>%
  dplyr:: summarise(length(Variant_Type)) %>%
  spread(key = "Variant_Type", value = "length(Variant_Type)") %>% 
  as.data.frame()
Mutation_numbers[is.na(Mutation_numbers)] <- 0

Mutation_numbers$TOTAL <- rowSums(Mutation_numbers[!'%in%'(names(Mutation_numbers), "Patient.ID")])

mut_clin <- Mutation_numbers %>% 
  merge(., pat_sub[, c("Patient.ID", "Subtype")], by = "Patient.ID") %>%
  gather(., -"Patient.ID", -"Subtype", key = "Variant", value = "Number") %>% droplevels()

mut_clin <- factorthese(mut_clin, c("Patient.ID", "Variant"))

# Plot
## All
for(i in levels(mut_clin$Variant)){
  print(i)
  work <- droplevels(subset(mut_clin, Variant == i))
  work$Rank <- log(work$Number + 1)
  temp_plot <- ggplot(work, aes(x = Subtype, y = Rank)) +
    geom_boxplot(alpha = 0.5, width = 0.2) + 
    geom_violin(aes(Subtype, fill = Subtype), scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) +
    labs(x = "Subtype", y = paste0("ln(", i, " mutations + 1)")) +
    theme_bw() +
    theme(axis.text = element_text(size = 16)) +
    theme(legend.direction = "horizontal", legend.position = "top") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif")
  filen <- paste0(i, ".pdf")
  ggplot2:: ggsave(filen, plot = temp_plot, device = "pdf", path = "./Figures/Mutation/Numbers",
                   height = 6, width = 6)}

# Remove silent
# tcga_mut[, grep("Silent", ignore.case = T, tcga_mut)] %>% head()

Mutation_numbers <- tcga_mut %>%
  dplyr:: group_by(Patient.ID, Consequence) %>%
  dplyr:: summarise(length(Consequence)) %>%
  spread(key = "Consequence", value = "length(Consequence)") %>% 
  as.data.frame()
Mutation_numbers[is.na(Mutation_numbers)] <- 0

Mutation_numbers$TOTAL <- rowSums(Mutation_numbers[!'%in%'(names(Mutation_numbers), "Patient.ID")])

mut_clin <- Mutation_numbers %>% 
  merge(., pat_sub[, c("Patient.ID", "Subtype")], by = "Patient.ID") %>%
  gather(., -"Patient.ID", -"Subtype", key = "Consequence", value = "Number")


# Plot
## All
mut_clin$Consequence <- as.factor(mut_clin$Consequence)
for(i in levels(mut_clin$Consequence)){
  print(i)
  work <- droplevels(subset(mut_clin, Consequence == i))
  work$Rank <- log(work$Number + 1)
  temp_plot <- ggplot(work, aes(x = Subtype, y = Rank)) +
    geom_boxplot(alpha = 0.5, width = 0.2) + 
    geom_violin(aes(Subtype, fill = Subtype), scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) +
    labs(x = "Subtype", y = paste0("ln(", i, " mutations + 1)")) +
    theme_bw() +
    theme(axis.text = element_text(size = 16)) +
    theme(legend.direction = "horizontal", legend.position = "top") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif")
  filen <- paste0(i, ".pdf")
  ggplot2:: ggsave(filen, plot = temp_plot, device = "pdf", path = "./Figures/Mutation/Numbers/Consequence",
                   height = 6, width = 6)}
