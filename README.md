# README File for the interrogation of the TCGA COADREAD Dataset

## Delineating a hiCIRC cohort
### Generating the Clinical dataframe
Run the first part of the FPKM.R script to get the Patient_list.txt document 
Run the Clinical.R script, to source the Clinical data for the 614 patients. This data has the microsatellite status for the patients of interest, this data was collated from two sources:
- TCGA_public_clinical.csv ( downloaded from cBioPortal)
- Pan-cancer immunogenomic perspective on the tumor microenvironment based on PD-L1 and CD8 T cell infiltration

### FPKM CIRC enrichment
This part of the script begins at a commented sections titled " # START CIRC Enrichment ------ "
It begins by investigating the CIRC enrichment score in the two microsatellite status groups - MSI-H and MSS, a strong significant difference is found.
However, the MSS cohort is highly variable, as shown by both the fligner test and the aymptotic tests of differences in variance.
Save the image to create "FPKMs.RData"

Attempts in clustering use both kmeans and Phenograph clustering to find clusters of MSS patients which have a high CIRC enrichment score, similar to that of MSI-H
Those patients which have a high CIRC score, are written in the Patient_Subtypes.csv file (in Output/ folder)

## Interrogating the cohort based on FPKM values
Starts at "START Interrogation of the cohort ----"

### Interrogation of Genesets
- Interrogate Ping and CellType/Signature enrichment scores as both Pearson correlation and enrichment across the subtypes
- Interrogate GO datasets for genesets of interest (Response to Reactive Oxygen Species)

### Interrogation of individual genes
- Requires merging of FPKM and pat_sub into MA dataframe
- For correlation of Genes of Interest with the CIRC Score, remove MSI-H first 
- For comparisons across Subtype, they can stay...

- Part 1: list of genes of interest across the cohorts
- Part 2: Bespoke genes.
-- Add genes that are of interest from Part 2 into Part 1

### Last part of the script investigates whatever Gary wants me to investigate...




## Proving that MHC Class II is the correct way to delineate these patients.
Script called "Class II impact" 
- Dependency is "FPKMs.RData" 
- Label the data based on the CIRC

### Perform logistic regression (simple)
Generalised linear model picks out 4 HLA class II genes. and a McFadden index of 0.33
Model just based on the Class II genes in the signature had a McFadden index of 0.23

#### ROC analysis
Using just the Class II genes in the CIRC signature, and a naive Bayes classifier 97% of MSS patients are correctly assigned to MSS-hiCIRC or MSS
- This is an average of 100 random training and test sets



To be continued...

