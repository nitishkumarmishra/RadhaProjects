### This program will select top feature and build model by using multivariate lm
library(ROCR); library(MASS)
library(dplyr)
setwd("F:/OneDrive - University of Nebraska Medical Center/Radha/Autism/NewAutism/")
data <- read.csv("7-Autism Blood Individual 230 CpG AUC Aug 10 2018 to Nitish.txt", header = TRUE, row.names = 2, sep = "\t")
select <- grep("AVG_Beta", names(data), value = TRUE)
data <- data[, select]
colnames(data) <- gsub("X", "", colnames(data))
colnames(data) <- gsub(".AVG_Beta", "", colnames(data))
#colnames(Meth) <- gsub("PanCancer", "", colnames(Meth))
data1 <- as.data.frame(t(data))
data1$Status <- as.factor(ifelse(grepl("Affected", rownames(data1)), "Autism", "Normal"))# we need class as factor for SVM 


col_idx <- grep("Status", names(data1))
data1 <- data1[, c(col_idx, (1:ncol(data1))[-col_idx])] ## Move last (Status) at first
## for msvmRFE we need class in first colum and in for of factor

set.seed(12345)
library(e1071)
source('msvmRFE.R')

svmRFE(data1, k=10, halve.above=100)


nfold = 10
nrows = nrow(data1)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds
folds = lapply(1:nfold, function(x) which(folds == x))
folds

# Perform feature ranking on all training sets
results = lapply(folds, svmRFE.wrap, data1, k=10, halve.above=100)
length(results)
results

top.features = WriteFeatures(results, data1, save=F)
head(top.features, 17)$FeatureName ## We use different number of feature and finally find 17 CpG can achieve AUC 1
#################################################################


data <- read.csv("7-Autism Blood Individual 230 CpG AUC Aug 10 2018 to Nitish.txt", header = TRUE, row.names = 2, sep = "\t")
select <- grep("AVG_Beta", names(data), value = TRUE)
data <- data[, select]
colnames(data) <- gsub("X", "", colnames(data))
colnames(data) <- gsub(".AVG_Beta", "", colnames(data))
#colnames(Meth) <- gsub("PanCancer", "", colnames(Meth))
data1 <- as.data.frame(t(data)) ## for lm we need Status 1 or 0
data1$Status <- ifelse(grepl("Affected", rownames(data1)), 1, 0)# we need class as factor for SVM 
#model.all <- glm(Status ~ cg20129082+cg08590939+cg20187719+cg10989317+cg15028160+cg04663916+cg21341586+cg15371711+cg05806645+cg05317207+ch.2.1116759R+cg02648941+cg07317062+cg10218876+cg05468028+cg04674383+cg26748794, data = data1)
model.all <- glm(Status ~ cg20129082+cg08590939+cg20187719, data = data1)
step <- stepAIC(model.all, direction = "both")
step$anova ## Select the list of features based on The Akaike information criterion (AIC)
model.stepAIC <- glm(Status ~ cg08590939 + cg20187719 + cg10989317 + cg15028160 + cg04663916 + cg21341586 + cg15371711 + cg05806645 + cg05317207 + ch.2.1116759R + cg07317062 + cg10218876 + cg04674383, data=data1)

#step <- stepAIC(model.all, direction = "both")
#step$anova ## Select the list of features based on The Akaike information criterion (AIC)
#model.stepAIC <- glm(Status ~ cg10979903 + cg26814276 + cg18073151 + cg11125369 + cg25999722 + cg22943986 + cg12853563 + cg25729826 + cg21241839, data=data1)

se_auc <- function(auc, n1, n2) { #n1 is positive, n2 is negative
  q1 <- auc / (2 - auc)
  q2 <- (2 * auc ^ 2) / (1 + auc)
  sqrt((auc * (1 - auc) + (n1 - 1) * (q1 - auc ^ 2) + (n2 - 1) * (q2 - auc ^ 2)) / (n1 * n2))
}

ci_auc <- function(auc, n1, n2, level = 0.95) {
  ci <- auc + c(-1, 1) * qnorm((1 + level) / 2) * se_auc(auc, n1, n2)
  ci[1] <- ifelse(ci[1] < 0, 0, ci[1])
  ci[2] <- ifelse(ci[2] > 1, 1, ci[2])
  ci }
save_roc_plot <- function(roc_obj) {}

rocs <- performance(prediction(fitted(model.all), data1$Status), "tpr", "fpr")
aucs <- performance(prediction(fitted(model.all), data1$Status), "auc")
roc <- rocs
auc <- aucs@y.values[[1]]
attributes(roc)$alpha.values[[1]][1:4] <- 1
#attributes(roc)$alpha.values[[1]][1:5] <- 1 ## Alpha value >1 then replace it with 1. Here first five are more than 1 
#attributes(roc)$alpha.values[[1]][21:25] <- 0 ## Replace alpha vale < 1 with zero, here last five have bnegative values
ci <- ci_auc(auc, 12, 12) # Here we have 12 affcted and 12 control samples
ci_chr <- paste0(round(ci, 2), collapse = ", ") 
auc_df <- data.frame(name = c("Mulivariate"),
                     AUC = numeric(length(1)),
                     CI_lower = numeric(length(1)),
                     CI_upper = numeric(length(1)), stringsAsFactors = F)
auc_df[1, c("AUC", "CI_lower", "CI_upper")] <- c(auc, ci[1], ci[2])
auc_df

save_roc_plot <- function(roc_obj) {}
## PDF image quality is better, so I will save figure in PDF
#file_name = paste0("ROC_mulivariate", "CpG", ".png")    ### change directory name here (the one you created above.
#png(filename = file_name, width = 1024, height = 768)
file_name = paste0("ROC_mulivariate 3 ", "CpG", ".pdf")    ### change directory name here (the one you created above.
pdf(file_name)
#plot(roc, main = name, col = "red")
plot(roc, main = "Multiple linear regression", colorize=T)  
abline(0, 1, lty = 2)
### Adjust position based on figure (pdf or png)
text(0.80, 0.13, paste0("AUC = ", round(unlist(auc), 2)), cex = 1.4)
text(0.80, 0.08, paste0("95% CI [", ci_chr, "]"), cex = 1.2)
rect(0.65, 0.05, 0.95, 0.16)
dev.off()



