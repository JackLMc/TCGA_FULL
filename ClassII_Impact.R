library(UsefulFunctions)
library(tidyverse)

cbcols <- c("MSS-hiCIRC" = "#999999",
            "MSI-H" = "#56B4E9",
            "MSS" = "#009E73")

my_comparisons <- list(c("MSS-hiCIRC", "MSI-H"),
                       c("MSS-hiCIRC", "MSS"),
                       c("MSI-H", "MSS"))

#### Read in data and process ####
load("FPKMs.RData")

# CIRC Geneset
CIRC_IG <- read.csv("./Exploratory_Data/Genesets/CIRC.csv")
CIRC_IG$Hugo_Symbol <- as.factor(CIRC_IG$Hugo_Symbol)

# Clinical
patient_subtypes <- read.csv("./Output/Patient_Subtypes.csv")
Clin_614 <- read.csv("./Output/Clinical_Data_614.csv") # REDO THIS

#### Label CIRC genes ####
CIRC_genes <- droplevels(subset(CIRC_IG, CIRC == T)) %>%
  takegenelevels()

FPKM2$CIRC <- ifelse(
  (FPKM2$SYMBOL %in% CIRC_genes), T, F)
doCIRC <- droplevels(subset(FPKM2, CIRC == T))

#### Process to get in format for PCA ####
# Remove uneeded stuff
pca1 <- doCIRC %>% gather(contains("TCGA"), key = "Patient.ID", value = "FPKM") %>%
  dplyr:: select(matches("Patient.ID"), matches("SYMBOL"), matches("FPKM")) %>%
  spread(key = "SYMBOL", value = "FPKM") %>% 
  merge(patient_subtypes, by = "Patient.ID")# Merge with cleaned clinical

pca3 <- droplevels(subset(pca1, Subtype == "MSS" | Subtype == "MSS-hiCIRC"))

# make the Patient ID the row names
pca4 <- data.frame(pca3[, names(pca3) != "Patient.ID"], row.names = pca3[, names(pca3) == "Patient.ID"])

#### Complete PCA ####
## Check rows with NA and columns
prin_comp <- prcomp(pca4[, !'%in%'(names(pca4), c("Subtype", "CIRC_Genes"))], scale. = T)
Subtype <- pca4[, "Subtype"]
library(ggbiplot)

# pdf("./Figures/Contribution_ClassII/PCA_of_CIRC_Gene_List_MSS.pdf", height = 6, width = 6)
# ggbiplot(prin_comp, obs.scale = 1, var.scale = 1, groups = Subtype,
#          ellipse = T, circle = T, var.axes = F) +
#   theme_bw() +
#   theme(legend.direction = "horizontal", legend.position = "top") +
#   #ylim(-5, 5) +
#   scale_colour_manual(values = cbcols) +
#   theme(axis.text = element_text(size = 16)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# dev.off()


# Rtsne
library(Rtsne)
set.seed(1)
tsne_out <- Rtsne(as.matrix(pca4[, !'%in%'(names(pca4), c("Subtype", "CIRC_Genes"))]))
tsne_dimensions <- as.data.frame(tsne_out$Y)
colnames(tsne_dimensions) <- c("Dim1", "Dim2")

## tSNE plot
# pdf("./Figures/Clustering/MSS_tsNE.pdf")
ggplot(tsne_dimensions, aes(x = Dim1, y = Dim2, colour = Subtype)) +
  geom_point(size = 4, alpha = 0.8, pch = 20) +
  scale_colour_manual(values = c("MSS" = "#009E73", "MSS-hiCIRC" = "#999999")) +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.direction = "horizontal", legend.position = "top")
# dev.off()


#### end ####

#### Investigate PCA ####
# How many PC's are important?
library(factoextra)
pdf("./Figures/Contribution_ClassII/Screeplot_CIRC_List.pdf",
    width = 6, height = 6)
fviz_screeplot(prin_comp, ncp = 10, choice = "eigenvalue", ggtheme = theme_bw(), main = "") +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 18)) +
  geom_hline(yintercept = 1, linetype = "dashed")
dev.off()

# Find the Eigenvalues
eig <- (prin_comp$sdev)^2

## Variances in percentage
variance <- eig*100/sum(eig)

## Cumulative variances
cumvar <- cumsum(variance)

## Store the variances
var <- get_pca_var(prin_comp)

## Find the correlation between variables and principal components
loadings <- prin_comp$rotation

### Find orientation of loadings
# loads <- as.data.frame(loadings)
# loads[order(loads$PC1, decreasing = T)[1:10],]
# loads1 <- tibble:: rownames_to_column(loads, "Parameter")
# loads_contrib <- droplevels(subset(loads1, PC1 > 0.1))
# loads_contrib$Parameter <- as.factor(loads_contrib$Parameter)
# levels(loads_contrib$Parameter)

sdev <- prin_comp$sdev

### Find the correlataions
var.coord <- t(apply(loadings, 1, var_cor_func, sdev))

## Calculate the Cos2 (square of the coordinates)
var.cos2 <- var.coord^2

## Contribution of variables to the Components ((cos2) * 100) / (total cos2 of the component)
comp.cos2 <- apply(var.cos2, 2, sum)
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
contrib.var <- as.data.frame(var.contrib)

## Find the most contributing variable
Df1 <- tibble:: rownames_to_column(contrib.var, var = "Gene")
Contribute_to_CIRC <- Df1[order(Df1$PC1, decreasing = T),]
# write.csv("./Output/Contribute_to_CIRC.csv", x = Contribute_to_CIRC, row.names = F)

## Calculate the cumulative contributions to most important PC
PC <- paste("PC", seq(1:28), sep = "")
variance_for_PC <- as.data.frame(cbind(PC, variance))

Df2 <- Df1 %>%
  gather(contains("PC"), key = "PC", value = "variance") %>%
  spread(key = "Gene", value = "variance") %>% column_to_rownames(., var = "PC")
Df2a <- Df2[, !'%in%'(names(Df2), c("PC", "variance"))]/100
Df2b <- rownames_to_column(Df2a, var = "PC")

Df3 <- merge(variance_for_PC, Df2b, by = "PC")
Df3$variance <- as.numeric(as.character(Df3$variance))

## Calculate percentage to scale each PC not to be 100%
Df4 <- Df3 %>% dplyr:: mutate_each(funs(.*variance), -matches("PC|variance"))
all.equal(rowSums(Df4[-c(1, 2)]), Df3$variance) # check that the cumulative addition of all genes = total variance explained by that PC


## Take only the eigenvalues above 1 (explains more variance than a single gene alone)
eig_above_1 <- droplevels(subset(Df4, PC == "PC1" | PC == "PC2" |
                                   PC == "PC3" | PC == "PC4"))


sum(eig_above_1$variance) # Total variance that's being covered by these 4 PCs
all.equal(sum(eig_above_1$variance), sum(colSums(eig_above_1[-c(1,2)]))) # Check cumulative still = to total


# Add up the cumulative for each gene to the top 4
cumulative_contrib <- as.data.frame(colSums(eig_above_1[,-c(1,2)])) %>% rownames_to_column(., var = "Gene")
names(cumulative_contrib) <- c("Gene", "Contribution")

## Reorder to most contribution
cumulative_contrib$Gene <- reorder(cumulative_contrib$Gene,-cumulative_contrib$Contribution) 

pdf("./Figures/Contribution_ClassII/Cumulative_top4_PC.pdf", height = 6, width = 6)
ggplot(cumulative_contrib, aes(y = Contribution, x = Gene)) +
  geom_bar(stat = "identity", fill = "#4783B5") + 
  theme_bw() +
  theme(legend.direction = "horizontal", legend.position = "top") +
  scale_fill_manual(values = "blue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("Cumulative PC1-4 variance (%)") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3.5)) 
dev.off()

## How much does Class II contribute to this?
all_class_2 <- dplyr:: filter(cumulative_contrib, grepl("HLA", Gene))
sum(all_class_2$Contribution)
stained_class_2 <- dplyr:: filter(all_class_2, grepl("DR|DP|DQ", Gene))
sum(stained_class_2$Contribution)


#### Plot the PCs and colour if a class II molecule ####
# PCA
CIRC_Genetypes <- read.csv("./Exploratory_Data/Genesets/CIRC_Genetypes.csv")

ready_plot <- Contribute_to_CIRC %>% dplyr:: select(matches("Gene"), matches("PC1$"), matches("PC2$")) %>% 
  merge(CIRC_Genetypes, by = "Gene")
ready_plot$Gene_type <- as.factor(ready_plot$Gene_type)

## Colours
cbPalette <- c("Adhesion" = "#999999", 
               "Chemokine" = "#E69F00", 
               "Class_2" = "#56B4E9", 
               "Costimulatory" =  "#009E73",
               "Immune_Checkpoint" = "#CC79A7",
               "TCR_Related" = "#0072B2",
               "Th1_Function" = "#D55E00")


## Plotting
library(ggrepel)
pdf("./Figures/Contribution_ClassII/Scatterplot_of_PC_contributions.pdf", height = 7, width = 7)
ggplot(ready_plot, aes(PC1, PC2, color = Gene_type, label = ready_plot$Gene)) +
  labs(x = "Contribution to PC1",
       y = "Contribution to PC2") +
  theme_bw() +
  scale_color_manual(values = cbPalette) + 
  geom_label_repel() + 
  geom_point() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

# Cumulative contribution plot (genes)
fill <- c("#5F9EA0", "#E1B378","#0072B2","#CC79A7")
ready_plot1 <- ready_plot %>% 
  mutate(Total = PC1 + PC2)

ready_plot1$Gene <- reorder(ready_plot1$Gene,-ready_plot1$Total) 
ready_plot2 <- ready_plot1 %>% gather(contains("PC"), key = "Parameter", value = "Contribution")

ready_plot2$Parameter <- factor(ready_plot2$Parameter, levels = c("PC3", "PC2", "PC1"))

pdf("./Figures/Contribution_ClassII/Contribution_PC1_PC2_Individual.pdf", height = 6, width = 6)
ggplot() + geom_bar(aes(y = Contribution, x = Gene, fill = Parameter),
                    data = ready_plot2, stat = "identity") + 
  theme_bw() +
  theme(legend.direction = "horizontal", legend.position = "top") +
  scale_fill_manual(values = fill) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


# Cumulative contribution plot (Genetype groups)
## Gather PCs
ready_plot$Gene <- as.factor(ready_plot$Gene)
No_genes <- data.frame(Gene_type = character(),
                       Counts = double(),
                       stringsAsFactors = F)

c <- 1
for(i in levels(ready_plot$Gene_type)){
  print(i)
  work <- droplevels(subset(ready_plot, Gene_type == i))
  cnts <- nlevels(work$Gene)
  No_genes[c, "Gene_type"] <- i
  No_genes[c, "Counts"] <- cnts
  c <- c + 1
}

ready_plot3 <- ready_plot %>% 
  gather(contains("PC"), key = "Parameter", value = "Contribution") %>%
  dcast(Gene_type ~ Parameter, sum, value.var = "Contribution") %>% 
  merge(., No_genes, by = "Gene_type") %>%
  mutate(PC1 = PC1/Counts) %>% mutate(PC2 = PC2/Counts) %>%
  mutate(Total = PC1 + PC2) 

## Reorder
ready_plot3$Gene_type <- factor(ready_plot3$Gene_type,
                                levels = ready_plot3$Gene_type[order(-ready_plot3$Total)],
                                ordered = T)

## Get in format 
ready_plot4 <- ready_plot3 %>% dplyr:: select(-matches("Total")) %>% 
  gather(contains("PC"), key = "Parameter", value = "Contribution")

## Plot
### Reorder PC values so PC1 is on the bottom 
fill <- c("#5F9EA0", "#E1B378")
ready_plot4$Parameter <- factor(ready_plot4$Parameter, levels = c("PC2", "PC1"))

pdf("./Figures/Contribution_ClassII/Contribution_PC1_PC2_Genesets.pdf", height = 6, width = 6)
ggplot() + geom_bar(aes(y = Contribution, x = Gene_type, fill = Parameter), 
                    data = ready_plot4, stat = "identity") +
  xlab("Gene Type") +
  ylab("Normalised Contribution") +
  theme_bw() +
  scale_fill_manual(values = fill) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.direction = "horizontal", legend.position = "top") +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

#### ROC Curve ####
# Logistic regression - MSS versus hiCIRC MSS
## Clean data
pca1 <- doCIRC %>% gather(contains("TCGA"), key = "Patient.ID", value = "FPKM") %>%
  dplyr:: select(matches("Patient.ID"), matches("SYMBOL"), matches("FPKM")) %>%
  spread(key = "SYMBOL", value = "FPKM") %>% 
  merge(patient_subtypes, by = "Patient.ID") # Merge with cleaned clinical


## hiCIRC versus MSS
my_data <- pca1 %>%
  filter(Subtype == "MSS" | Subtype == "MSS-hiCIRC") %>%
  droplevels() %>%
  #select(matches("HLA|MSI|Patient")) %>%
  na.omit()


duplicated(my_data$Patient.ID)
contrasts(my_data$Subtype)

## Remove uneeded columns
my_data <- data.frame(my_data[, names(my_data) != "Patient.ID"],
                      row.names = my_data[, names(my_data) == "Patient.ID"])
# my_data <- my_data %>% dplyr:: select(matches("HLA|Subtype"))
my_data <- my_data[!('%in%'(colnames(my_data), c("CIRC_Genes")))]

## Straight forward model
model <- glm(Subtype ~.,family = binomial(link = "logit"), data = my_data, control = list(maxit = 50))
summary(model)
library(MASS)
stepAIC(model)
bet_model <- glm(formula = Subtype ~ CCL5 + CXCL10 + CXCL9 + HAVCR2 + HLA.DOA + 
                   HLA.DPA1 + HLA.DRA + HLA.DRB5 + ICAM1 + LAG3 + PDCD1LG2 + 
                   STAT1, family = binomial(link = "logit"), data = my_data, control = list(maxit = 50))
stepAIC(bet_model)
bet_model1 <- glm(formula = Subtype ~ CXCL10 + HAVCR2 + HLA.DPA1 + HLA.DRA + 
      HLA.DRB5 + ICAM1 + LAG3 + PDCD1LG2, family = binomial(link = "logit"), 
    data = my_data, control = list(maxit = 50))
stepAIC(bet_model1)
summary(bet_model1)
test_model<- glm(formula = Subtype ~  HLA.DPA1 + HLA.DPB1 + HLA.DQA1 + HLA.DQA2 + HLA.DRA + HLA.DRB5, family = binomial(link = "logit"), 
                 data = my_data)
stepAIC(test_model)
summary(test_model)
# While no exact equivalent to the R2 of linear regression exists, the McFadden R2 index can be used to assess the model fit.
library(pscl)
pR2(bet_model1)

#Get ORs and confidence intervals
exp(cbind(OR = coef(model), confint(model)))

#Get Chi-squared ANOVA P values between your groups
anova(model, test = "Chisq")

#Perform cross validation (CV) analysis
#The delta values should not greatly differ
#K=number of samples, i.e., leave-one-out CV.
library(boot)
cv.glm(my_data, model, K=nrow(my_data))$delta

#### Complex model for ROC Curve ####
# attempt
# 1 = MSS, 2 = MSS-hiCIRC
library(e1071)
library(ROCR)
lvls <- levels(my_data$Subtype)
smp_size <- floor(0.75 * nrow(my_data))

#### Trying to take an average of 100 seeds ####
aucs = c()
# df_test <- data.frame()
# g <- ggplot(df_test) + geom_line() + xlim(0, 1) + ylim(0, 1) +
#   xlab("False Positive Rate") +
#   ylab("True Positive Rate") +
#   theme_bw() +
#   theme(axis.text = element_text(size = 16)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# dev.off()

pdf("./Figures/Contribution_ClassII/ROC_hiCIRC_MSS_CIRC_CLASS2.pdf", height = 6, width = 6)
plot(x = NA, y = NA, xlim = c(0, 1), ylim = c(0, 1),
     ylab = "True Positive Rate",
     xlab = "False Positive Rate",
     bty = "n")


# 1 = MSS, 2 = hiCIRC

AUC_df = data.frame()
coordinate_list = list()

c <- 1
for(seed in 1:100){
  set.seed(seed)
  type.id <- 1
  testidx <- sample(seq_len(nrow(my_data)), size = smp_size)
  training <- my_data[testidx, ]
  testing <- my_data[-testidx, ]
  
  type <- as.factor(training$Subtype == lvls[type.id])
  
  nbmodel <- naiveBayes(type ~ ., data = dplyr:: select(training, matches("HLA"))) # How well does Class 2 predict CIRC
  nbprediction <- predict(nbmodel, dplyr:: select(testing, matches("HLA")), type = "raw")
  
  score <- nbprediction[, "TRUE"]
  actual.class <- testing$Subtype == lvls[type.id]
  
  pred <- prediction(score, actual.class)
  nbperf <- performance(pred, "tpr", "fpr")
  
  roc.x <- unlist(nbperf@x.values)
  roc.y <- unlist(nbperf@y.values)
  
  coordinate_list[[seed]] <- as.data.frame(cbind(X_coord = roc.x, Y_coord = roc.y))
  
  lines(roc.y ~ roc.x, col = alpha(cbPalette[type.id + 1], 0.2), lwd = 2)
  nbauc <- performance(pred, "auc")
  nbauc <- unlist(slot(nbauc, "y.values"))
  
  AUC_df[c, "Seed_number"] <- seed
  AUC_df[c, "ROC"] <- nbauc
  c <- c + 1
}

legend("bottomright", legend = c("Average AUC - 0.97"),
       col = c("#E69F00"), lty = 1, cex = 1.2, lwd = 2)

dev.off()
mean(AUC_df$ROC)
se <- function(x) sqrt(var(x)/length(x))
se(AUC_df$ROC)
NAME <- paste0("Seed_", 1:length(coordinate_list))
names(coordinate_list) <- NAME

#### Any Difference in Clinical Outcome? ####
my_data4 <- my_data %>% rownames_to_column(var = "Patient.ID") %>% dplyr:: select(-matches("HLA"))
my_data5 <- merge(my_data4, Clin_614, by = "Patient.ID")
my_data6 <- droplevels(subset(my_data5, OS_STATUS != "NC"))

my_data6$died <- ifelse((my_data6$OS_STATUS == "LIVING"), F, T)

# Survival
os.surv <- Surv(my_data6$OS_MONTHS, my_data6$died)
os.surv.fitted <- survfit(os.surv ~ my_data6$Subtype)
plot(os.surv.fitted,
     col = c("red","blue", "black", "green"),
     xlab = "Overall Survival (Months)",
     ylab = "Probability of Survival")

legend("bottomright",
       inset = 0.05,
       c("MSS-hiCIRC","MSI-H", "MSI-L", "MSS"),
       fill = c("red","blue", "black", "green"))

os.surv.fitted

#### end ####

#### Normalise results for stage? ####

HLA <- FPKM2[grepl("HLA", FPKM2$SYMBOL), ] %>%
  droplevels()



HLA$SYMBOL

