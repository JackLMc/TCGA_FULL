# A script to evaluate the impact of Class II

library(UsefulFunctions)
library(tidyverse)
library(ggpubr)

cbcols <- c("MSS-hiCIRC" = "#999999",
            "MSI-H" = "#56B4E9",
            "MSS" = "#009E73")

my_comparisons <- list(c("MSS-hiCIRC", "MSI-H"),
                       c("MSS-hiCIRC", "MSS"),
                       c("MSI-H", "MSS"))

# Read in data and process ----
load("./R_Data/Counts_clean.RData")

# Don't se a seed, as lower down the seed is set again and again for the loop

## CIRC Geneset
CIRC_IG <- read.csv("./Exploratory_Data/Genesets/CIRC.csv")
CIRC_IG$SYMBOL <- as.factor(CIRC_IG$SYMBOL)

## Clinical
pat_sub <- read.csv("./Output/Patient_Subtypes_09_03.csv")
Clin_542 <- read.csv("./Output/Clinical_Data_542.csv")

## Label CIRC genes 
CIRC_genes <- droplevels(subset(CIRC_IG, CIRC == T)) %>%
  takegenelevels()

LongCQN <- Counts_cqn %>% as.data.frame() %>%
  rownames_to_column(., var = "Gene") %>%
  gather(contains("TCGA"), key = "Patient.ID", value = "Enrich") %>%
  spread(., key = "Gene", value = "Enrich") %>%
  gather(-Patient.ID, key = "Gene", value = "CQN")

LongCQN <- Counts_cqn[rownames(Counts_cqn) %in% CIRC_genes, ]

# ROC Curve ----
# Logistic regression - MSS versus hiCIRC MSS
## Clean data
pca1 <- LongCQN %>% as.data.frame() %>%
  rownames_to_column(., var = "Gene") %>%
  gather(contains("TCGA"), key = "Patient.ID", value = "CQN") %>%
  spread(key = "Gene", value = "CQN") %>% 
  merge(pat_sub, by = "Patient.ID")  # Merge with cleaned clinical

my_data <- pca1 %>%
  filter(Subtype == "MSS" | Subtype == "MSS-hiCIRC") %>%
  droplevels() %>%
  #select(matches("HLA|MSI|Patient")) %>%
  na.omit() # hiCIRC versus MSS

# my_data[duplicated(my_data$Patient.ID)]
contrasts(my_data$Subtype)

my_data <- data.frame(my_data[, names(my_data) != "Patient.ID"],
                      row.names = my_data[, names(my_data) == "Patient.ID"])
# my_data <- my_data %>% dplyr:: select(matches("HLA|Subtype"))
my_data <- my_data[!('%in%'(colnames(my_data), c("CIRC_Genes")))]

## Straight forward model
model <- glm(Subtype ~.,family = binomial(link = "logit"), data = my_data, control = list(maxit = 50))

summary(model)
library(MASS)
stepAIC(model)
bet_model <- glm(formula = Subtype ~ CCL5 + CD4 + CXCL9 + GNLY + HLA.DOA + 
                   HLA.DPB1 + HLA.DRA + ICAM1 + ICOS + IFNG + IL18RAP + LAG3, 
                 family = binomial(link = "logit"), data = my_data, control = list(maxit = 50))
stepAIC(bet_model)
summary(bet_model)


test_model <- glm(formula = Subtype ~  HLA.DPA1 + HLA.DPB1 + HLA.DQA1 + HLA.DQA2 + HLA.DRA + HLA.DRB5, family = binomial(link = "logit"), 
                 data = my_data, control = list(maxit = 50))
stepAIC(test_model)
summary(test_model)

# While no exact equivalent to the R2 of linear regression exists, the McFadden R2 index can be used to assess the model fit.
library(pscl)
pR2(bet_model)
pR2(test_model)

#Get ORs and confidence intervals
exp(cbind(OR = coef(model), confint(model)))

#Get Chi-squared ANOVA P values between your groups
anova(bet_model, test = "Chisq")

# Perform cross validation (CV) analysis
# The delta values should not greatly differ
# K=number of samples, i.e., leave-one-out CV.
library(boot)
cv.glm(my_data, bet_model, K=nrow(my_data))$delta



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
  
  lines(roc.y ~ roc.x, col = alpha("#E69F00", 0.2), lwd = 2)
  nbauc <- performance(pred, "auc")
  nbauc <- unlist(slot(nbauc, "y.values"))
  
  AUC_df[c, "Seed_number"] <- seed
  AUC_df[c, "ROC"] <- nbauc
  c <- c + 1
}

legend("bottomright", legend = c("Average AUC - 0.94"),
       col = c("#E69F00"), lty = 1, cex = 1.2, lwd = 2)

dev.off()
mean(AUC_df$ROC)
se <- function(x) sqrt(var(x)/length(x))
se(AUC_df$ROC)
NAME <- paste0("Seed_", 1:length(coordinate_list))
names(coordinate_list) <- NAME


#### END ####


