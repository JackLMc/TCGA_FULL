getwd()
setwd("/Volumes/2018/beggsa-tcgacolorectal/download_rest/harvest/")


filelist1 <- list.files(full.names = F)
filelist = unlist(filelist1)


GDC_convert <- read.delim("~/OneDrive/UoB/PhD/Projects/6_TCGA_Full/TCGA_FULL/Data/PathSeq/gdc_sample_sheet.2019-04-10.tsv")
my_files <- GDC_convert[GDC_convert$File.Name %in% filelist, ]


hmmm <- filelist[!('%in%'(filelist, GDC_convert$File.Name))]
hmmm[!grepl("bai", hmmm)]





library(tidyverse)
files_with_gaps <- my_files[grepl("gap", my_files$File.Name), ] %>% droplevels()

head(files_with_gaps)





droplevels(subset(GDC_convert, Sample.ID == "TCGA-AG-A015-01A"))


