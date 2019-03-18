library(UsefulFunctions)
library(tidyverse)
library(ggpubr)
library(MSIseq)

# BiocManager:: install("RWeka")
library(RWeka)

my_comparisons <- list(c("MSS-hiCIRC", "MSI-H"),
                       c("MSS-hiCIRC", "MSS"),
                       c("MSI-H", "MSS"))
cbcols <- c("MSS-hiCIRC" = "#999999",
            "MSI-H" = "#56B4E9",
            "MSS" = "#009E73",
            "MSI-L" = "#E69F00")

url <-"http://steverozen.net/data/Hg19repeats.rda"
file <- basename(url)
download.file(url, file)
load("Hg19repeats.rda")

load("~/Downloads/NGStestdata.rda")
load("~/Downloads/NGStestseqLen.rda")

tail(NGStestseqLen)
tail(NGStestdata)


test.mutationNum <- Compute.input.variables(NGStestdata, repeats = Hg19repeats, seq.len = NGStestseqLen)

## find overlaps function
.intersects<-function(spl, ref){
  
  spl$Chrom = ifelse(substring(spl$Chrom,1,3)!='chr', paste('chr', spl$Chrom, sep=''), spl$Chrom)
  spl.out <- data.frame()
  for (chr in unique(spl$Chrom)){
    spl.iter <- spl[spl$Chrom==chr,]
    ref.iter <- ref[ref$Chrom==chr,]
    spl.range <- IRanges(spl.iter$Start_Position, spl.iter$End_Position)
    ref.range <- IRanges(ref.iter$Start_Position, ref.iter$End_Position)
    overlaps <- as.data.frame(findOverlaps(spl.range, ref.range, minoverlap=1))
    ins.filter <- as.data.frame(findOverlaps(spl.range, ref.range, minoverlap=2))
    MSI_status <- rep(0, nrow(spl.iter))
    for (i in which(spl.iter$Variant_Type!='INS')){
      MSI_status[i] = sum(overlaps$queryHits==i)
    }
    for (i in which(spl.iter$Variant_Type=='INS')){
      MSI_status[i] = sum(ins.filter$queryHits==i)
    }
    spl.iter.new <- cbind(spl.iter, MSI_status)
    spl.out <- rbind(spl.out, spl.iter.new)
  }
  return (spl.out) 
}

## function of getting mutation numbers
Compute.input.variables<-function(data, repeats, uniform.seq.len=38, seq.len = NULL){
  ## check if data has all the columns
  if(!all(c("Chrom", "Start_Position", "End_Position", "Tumor_Sample_Barcode", "Variant_Type")%in%colnames(data))){
    stop("Wrong column names in maf file.\nPlease check if the following columns are in your data:\nChrom, Start_Position, End_Position, Tumor_Sample_Barcode, Variant_Type.")    
  }
  ## check if Variant_Type has the right format;
  if(!all(unique(data$Variant_Type)%in%c("SNP", "INS", "DEL"))){
    stop("Wrong variant types. Only SNP, INS and DEL are accepted.")  
  }
  sorted.data<-.intersects(data, repeats)
  
  ## get the mutation numbers for each sample
  uni.sample<-unique(sorted.data$Tumor_Sample_Barcode)
  result.mat<-matrix(0,length(uni.sample),9) 
  rownames(result.mat)<-uni.sample
  colnames(result.mat)<-c("T.sns", "S.sns", "T.ind", "S.ind", "T", "S", "ratio.sns", "ratio.ind", "ratio")
  ## T: total mutation number
  ## R: mutation number in simple repeats
  
  for (i in seq_along(uni.sample)){
    result.mat[i,1]<-sum(sorted.data$Tumor_Sample_Barcode==rownames(result.mat)[i] & sorted.data$Variant_Type=="SNP")
    result.mat[i,2]<-sum(sorted.data$MSI_status[sorted.data$Tumor_Sample_Barcode==rownames(result.mat)[i] & sorted.data$Variant_Type=="SNP"])
    result.mat[i,3]<-sum(sorted.data$Tumor_Sample_Barcode==rownames(result.mat)[i] & sorted.data$Variant_Type!="SNP")
    result.mat[i,4]<-sum(sorted.data$MSI_status[sorted.data$Tumor_Sample_Barcode==rownames(result.mat)[i] & sorted.data$Variant_Type!="SNP"])
    result.mat[i,5]<-result.mat[i,1]+result.mat[i,3]
    result.mat[i,6]<-result.mat[i,2]+result.mat[i,4]
    result.mat[i,7]<-result.mat[i,2]/result.mat[i,1]
    result.mat[i,8]<-result.mat[i,4]/result.mat[i,3]
    result.mat[i,9]<-result.mat[i,6]/result.mat[i,5]
  }
  
  ## divide the absolute mutation number by length
  ## if there is no seq.len, use uniform.seq.len
  if (is.null(seq.len)) {
    result.mat[,1:6] = result.mat[,1:6]/uniform.seq.len
  } else { ## use seq.len
    ## seq.len contains 2 columns: Tumor_Sample_Barcode & Sequence_Length (Mb)
    if(!all(c("Tumor_Sample_Barcode", "Sequence_Length")%in%colnames(seq.len))){
      stop("Wrong column names in seq.len.")
    }
    ## check if the sample size is the same in maf and capture file
    if (nrow(result.mat)!=nrow(seq.len)|all(rownames(result.mat)%in%seq.len$Tumor_Sample_Barcode)==FALSE){
      stop("Samples do not match between seq.len and data.")
    }
    
    result.mat[,1:6] = result.mat[,1:6]/seq.len$Sequence_Length[match(rownames(result.mat), seq.len$Tumor_Sample_Barcode)]
  }
  mutationNum = as.data.frame(result.mat)
  mutationNum[is.na(mutationNum)] = 0
  mutationNum}








if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("maftools")

library(maftools)

muta <- read.maf("~/Downloads/mc3.v0.2.8.PUBLIC.maf")
pat_sub <- read.csv("./Data/Important/patient_subtypes.csv")

library(UsefulFunctions)
muta@data$Patient.ID <- samptopat(muta@data$Tumor_Sample_Barcode)
muta@data$Patient.ID <- gsub("-", ".", muta@data$Patient.ID)
muta@data <- muta@data[muta@data$Patient.ID %in% pat_sub$Patient.ID, ]

muta1 <- muta
pat_sub1 <- droplevels(subset(pat_sub, Subtype == "MSS-hiCIRC"))
nlevels(pat_sub1$Patient.ID)

muta1@data <- muta1@data[muta1@data$Patient.ID %in% pat_sub1$Patient.ID,]


lollipopPlot(muta1, gene = "APC")
geneCloud(muta1, minMut = 8)

