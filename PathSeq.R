# PathSeq analysis
library(UsefulFunctions)
library(tidyverse)
library(ggpubr)

my_comparisons <- list(c("MSS-hiCIRC", "MSI-H"),
                       c("MSS-hiCIRC", "MSS"),
                       c("MSI-H", "MSS"))
cbcols <- c("MSS-hiCIRC" = "#999999",
            "MSI-H" = "#56B4E9",
            "MSS" = "#009E73",
            "MSI-L" = "#E69F00")




# Looking for the files
## Mount the drive
thousand.folders <- list.dirs(path = "/Volumes/2018/beggsa-tcgacolorectal/download_rest/bacterial_project/results", full.names = T)
filelist1 <- sapply(thousand.folders, function(x){
  list.files(x, pattern = ".pathseq.txt$", full.names = T)})
filelist = unlist(filelist1)



# Read in files and combine
lists <- lapply(filelist, read.delim, header = T, stringsAsFactors = F)
listsDF <- lists
lists <- listsDF
names(lists) <- filelist
names(lists) <- gsub("/Volumes/2018/beggsa-tcgacolorectal/download_rest/bacterial_project/results/", "", names(lists))
names(lists) <- gsub("pathseq.txt", "bam", names(lists))

# Colnames = c(tax_id, taxonomy, type, name,
# kingdom, score, score_normalized, reads, unambiguous, reference_length)
GDC_convert <- read.delim("./Data/PathSeq/gdc_sample_sheet.2019-04-10.tsv")
GDC_convert$Patient.ID <- gsub("-", ".", GDC_convert$Case.ID)

# colnames(GDC_convert)[colnames(GDC_convert) == "cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id"] <- "file.ID"
GDC_convert_t <- droplevels(subset(GDC_convert, Sample.Type != "Blood Derived Normal" &
                                Sample.Type != "Solid Tissue Normal"))
GDC_convert_n <- droplevels(subset(GDC_convert, Sample.Type == "Blood Derived Normal" |
                                     Sample.Type == "Solid Tissue Normal"))


lists1 <- lists[names(lists) %in% GDC_convert_t$File.Name]
lists1[["TCGA-AA-3846-01A-01W-0995-10_Illumina_gdc_realn.bam"]]

lists1a <- lists1[sapply(lists1, function(x) dim(x)[1]) > 0]
lists2 <- lapply(names(lists1a), 
                  function(n, x){
                    x[[n]]$File.Name <- n
                    return (x[[n]])},
                 lists1)


multi_join <- function(list_of_loaded_data, join_func, ...){
  require("dplyr")
  output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
  return(output)
}

pathseq <- multi_join(lists2, full_join, by = c("tax_id", "taxonomy", "type", "name",
                                                    "kingdom", "score", "score_normalized",
                                                    "reads", "unambiguous", "reference_length", "File.Name"))

pathseq1 <- merge(pathseq, GDC_convert, by = "File.Name") %>% droplevels()
pathseq1$Patient.ID <- as.factor(pathseq1$Patient.ID)
pathseq1$kingdom <- as.factor(pathseq1$kingdom)
pathseq1$name <- as.factor(pathseq1$name)
pathseq1$taxonomy <- as.factor(pathseq1$taxonomy)
pathseq1$Patient.ID <- as.factor(pathseq1$Patient.ID)
pathseq1$type <- as.factor(pathseq1$type)


pathseq2 <- pathseq1[, c("Patient.ID", "name", "score",
                         "score_normalized", "reads", "kingdom", "type", "unambiguous", "reference_length")]

# Fill in the missing data for microbes not found in those samples
microbe_taxa <- pathseq1[, c("kingdom", "type", "name", "reference_length")]
microbe_taxa <- microbe_taxa[!duplicated(microbe_taxa), ]

missing_microbes <- list()
for(i in levels(pathseq2$Patient.ID)){
  work <- droplevels(subset(pathseq2, Patient.ID == i))
  pat_microbes <- as.character(levels(work$name))
  all_microbes <- as.character(levels(pathseq2$name))
  see <- all_microbes[!(all_microbes %in% pat_microbes)]
  missing <- ifelse((see != 0), c(see), "None")
  missing_microbes[[i]] <- missing
  }

missing_microbes1 <- lapply(missing_microbes, as.data.frame, stringsAsFactors = F)
missing_microbes2 <- lapply(missing_microbes1, setNames, c("name"))

missing_microbes3 <- lapply(names(missing_microbes2), 
                 function(n, x){
                   x[[n]]$Patient.ID <- n
                   return (x[[n]])},
                 missing_microbes2)

missing_microbes3 <- mapply(cbind, missing_microbes3, "score" = 0, SIMPLIFY = F)
missing_microbes3 <- mapply(cbind, missing_microbes3, "score_normalized" = 0, SIMPLIFY = F)
missing_microbes3 <- mapply(cbind, missing_microbes3, "reads" = 0, SIMPLIFY = F)
missing_microbes3 <- mapply(cbind, missing_microbes3, "unambiguous" = 0, SIMPLIFY = F)


missing_data <- multi_join(missing_microbes3, full_join, by = c("Patient.ID", "name", "score", "score_normalized",
                                                "reads", "unambiguous"))
missing_data <- merge(missing_data, microbe_taxa, by = "name")
pathseq3 <- rbind(pathseq2, missing_data)

# ## Double-check
# output <- data.frame(stringsAsFactors = F)
# c <- 1
# for(i in levels(pathseq3$Patient.ID)){
#   work <- droplevels(subset(pathseq3, Patient.ID == i))
#   these <- nlevels(work$name)
#   output[c, "Patient.ID"] <- i
#   output[c, "num_of_tax"] <- these
#   c <- c + 1
# }

# Merge with the pat_sub file
pat_sub <- read.csv("./Output/Patient_Subtypes.csv")
pathseq4 <- merge(pathseq3, pat_sub, by = "Patient.ID")
pathseq4 <-  pathseq4[, c("Patient.ID", "kingdom", "type", "name",
                          "reference_length", "score", "score_normalized",
                          "reads", "unambiguous", "CIRC_Genes", "Subtype")] %>% droplevels()


library(reshape2)
 
droplevels(subset(pathseq4, Subtype == "MSS"))$Patient.ID %>% nlevels()
droplevels(subset(pathseq4, Subtype == "MSS-hiCIRC"))$Patient.ID %>% nlevels()
droplevels(subset(pathseq4, Subtype == "MSI-H"))$Patient.ID %>% nlevels()

empty_ <- data.frame("name" = character(),
                     "num_lev" = double(),
                     stringsAsFactors = F)
c <- 1
for(i in levels(pathseq4$name)){
  print(i)
  working <- droplevels(subset(pathseq4, name == i))
  empty_[c, "name"] <- i
  empty_[c, "num_lev"] <- nlevels(working$Patient.ID)
  c <- c + 1
  }

empty_[empty_$num_lev != 529, ]


try <- pathseq4[grepl("Clostridium", pathseq4$name, ignore.case = T), ] %>% droplevels()

Ecol <- droplevels(subset(pathseq4, name == "Buty"))


# Coprobacillus_sp._D6, Bifidobacterium_animalis_subsp._lactis (not exact to one in paper), Clostridium_sp._ASF356, Clostridium_saccharobutylicum (not exact)
# Eubacterium_sp._3_1_31 (up MSI-h), Erysipelotrichaceae_bacterium_21_3 (up in MSI-H, slightly), Firmicutes_bacterium_ASF500,
# Clostridium_hathewayi 12489931 (not in), Ruminococcus_gnavus_AGR2154 (MSI-H), Ruminococcus bromii (not in) or Ruminococcus obeum
# Subdoligranulum_sp._4_3_54A2FAA (highly expressed in MSI-H), Bifidobacterium_breve (not eact), Clostridium symbiosium (not in)
# Bacteroides_dorei (not exact, but up in MSI-H)



# Combine stats into a list, not working because some bugs are 0
maximum_nums <- data.frame(stringsAsFactors = F)
c <- 1
for(i in levels(pathseq4$name)){
  print(i)
  work <- droplevels(subset(pathseq4, name == i))
  max_num <- max(work$score)
  maximum_nums[c, "name"] <- i
  maximum_nums[c, "max"] <- max_num
  c <- c + 1
}

remove_these <- droplevels(subset(maximum_nums, max == 0))$name
pathseq5 <- pathseq4[!('%in%'(pathseq4$name, remove_these)), ] %>% droplevels()

# save.image("./PathSeq/PathSeq.RData")
load("./PathSeq/PathSeq.RData")

library(ggpubr)
stat_list <- list()
c <- 1
for(i in levels(pathseq5$name)){
  name <- basename(i)
  cat("Processing", i, "\n")
  workingon <- droplevels(subset(pathseq5, name == i))
  # assign your ggplot call to the i"th position in the list
  x <- compare_means(score ~ Subtype, data = workingon, method = "wilcox.test")
  y <- as.data.frame(x)
  stat_list[[i]]  <- cbind(i, y)
  c <- c + 1}


# Bind and remove row names
z <- do.call(rbind, stat_list)
rownames(z) <- c()

# Separate
CellDensity_Megapixel_Stat <- z 

# Write out the statistics
write.csv("Output/CellDensity_Megapixel_Stat.csv", x = CellDensity_Megapixel_Stat, row.names = F)

  
Ecol$Rank <- rank(Ecol$score)

ggplot(Ecol, aes(x = Subtype, y = Rank)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(aes(Subtype, fill = Subtype),
              scale = "width", alpha = 0.8) +
  scale_fill_manual(values = cbcols) +
  labs(x = "Subtype", y = "Normalised bact") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top") + 
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test")

# Peptostreptococcus, fusobacterium, Parvimonas, Lachnospiraceae all up in MSI-H patients...

