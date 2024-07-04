### This program will select top feature and build model by using multivariate lm
library(ROCR); library(MASS)
library(dplyr)
setwd("F:/OneDrive - University of Nebraska Medical Center/Radha/Cerebral Palcy/Latest/")
TableInput <- fread(input = "CP Analysis 341 Targets for AUC ROC Nov 23 2018 to Nitish regression.txt")
select <- grep("AVG_Beta|TargetID", names(TableInput), value = TRUE)
Meth <- TableInput[, select, with=FALSE]
colnames(Meth) <- gsub(".AVG_Beta", "", colnames(Meth))
Meth <- as.data.frame(Meth); Meth <- na.omit(Meth); rownames(Meth) <- Meth$TargetID; Meth$TargetID <- NULL
#colnames(data) <- gsub("Pre-eclampsia_", "", colnames(data))
data1 <- as.data.frame(t(Meth))
data1$Status <- as.factor(ifelse(grepl("Case", rownames(data1)), "Cerebral", "Normal"))# we need class as factor for SVM 


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


# data <- read.table("2-HS Individual-multivariate linear regression analysis-Aug 3 2018-To Nitish.txt", header = TRUE, row.names = 2)
# data[,c("AUC","Index", "CI_lower", "CI_upper")] <- NULL
# colnames(data) <- gsub("X", "", colnames(data))
# colnames(data) <- gsub(".AVG_Beta", "", colnames(data))
# data1 <- as.data.frame(t(data))
data1$Status <- ifelse(grepl("Case", rownames(data1)), 1, 0) ## for lm we need Status 1 or 0
model.all <- glm(Status ~ cg12425861 + cg13187827 + cg12204727 + cg24455365 + cg20415053 + cg08894153 + cg19499452 + cg03586379 + cg26707202 + cg20339553 + cg04761648 + cg05332869 + cg00167275 + cg08052428 + cg15277906 + cg08634464 + cg23000734, data = data1)
#model.all <- glm(Status ~ cg12425861 + cg13187827 + cg12204727 + cg24455365 + cg20415053 + cg08894153 + cg19499452 + cg03586379 + cg26707202 + cg20339553 + cg04761648 + cg05332869 + cg00167275 + cg08052428 + cg15277906, data = data1)
step <- stepAIC(model.all, direction = "both")
step$anova ## Select the list of features based on The Akaike information criterion (AIC)
model.stepAIC <- glm(Status ~ cg26707202+cg05332869+cg15277906+cg04761648+cg08052428+cg13187827+cg08634464, data=data1)

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
attributes(roc)$alpha.values[[1]][1:10] <- 1
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
file_name = paste0("ROC_mulivariate Top 17 ", "CpG", ".pdf")    ### change directory name here (the one you created above.
pdf(file_name)
#plot(roc, main = name, col = "red")
plot(roc, main = "Multiple linear regression", colorize=T)  
abline(0, 1, lty = 2)
### Adjust position based on figure (pdf or png)
text(0.80, 0.13, paste0("AUC = ", round(unlist(auc), 2)), cex = 1.4)
text(0.80, 0.08, paste0("95% CI [", ci_chr, "]"), cex = 1.2)
rect(0.65, 0.05, 0.95, 0.16)
dev.off()



