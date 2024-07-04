library(ROCR); library(MASS)
setwd("F:/OneDrive - University of Nebraska Medical Center/Radha/TOF Placenta/")

data <- read.table("TOF_Placenta IndividualData 58CpGs to Nitish for linear Regression July 31 2018.txt", header = TRUE, row.names = 2, sep = "\t")
data[,c("Index","ProbeID_A","ProbeID_B", "Cases.Average", "Controls.Average", "Diff_Avg.cases.and.controls", "CHR","UCSC_REFGENE_NAME", "ILMNID"  )] <- NULL
#colnames(data) <- gsub("X", "", colnames(data))
#colnames(data) <- gsub(".Pre.eclampsia_", "", colnames(data))
data <- as.data.frame(t(data))
data$Status <- ifelse(grepl("TOFCA", rownames(data)), 1, 0)

model.all <- glm(Status ~ cg00033213+cg00435526+cg02502145+cg02511231+cg03449867+cg03667593+cg03893271+cg04657470+cg05111616+cg05238741+cg05273049+cg06459104+cg06960824+cg06968724+cg07326438+cg08024264+cg08732526+cg09084244+cg09279736+cg09968723+cg09980477+cg10306192+cg10701801+cg10909185+cg11035303+cg11251367+cg11377136+cg12466610+cg12615916+cg13569207+cg13665593+cg13767940+cg14345128+cg15829294+cg16178625+cg16884400+cg17369694+cg17386473+cg17707870+cg17770035+cg18208707+cg18478319+cg19393008+cg19496978+cg19863210+cg20459037+cg20973720+cg21435568+cg21514997+cg21532325+cg22152407+cg23165899+cg24072924+cg25046571+cg26311944+cg26349606+cg27226920+cg27392792, data=data)
step <- stepAIC(model.all, direction = "both")
step$anova ## Select the list of features based on The Akaike information criterion (AIC)
model.stepAIC <- glm(Status ~ cg00033213+cg00435526+cg02502145+cg02511231 + cg03449867+cg03667593+cg03893271+cg04657470+cg05111616+cg05238741+cg05273049+cg06459104+cg06960824+cg06968724+cg07326438+cg08024264+cg08732526, data=data)

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

rocs <- performance(prediction(fitted(model.stepAIC), data$Status), "tpr", "fpr")
aucs <- performance(prediction(fitted(model.stepAIC), data$Status), "auc")
roc <- rocs
auc <- aucs@y.values[[1]]
#attributes(roc)$alpha.values[[1]][1] <- 1
attributes(roc)$alpha.values[[1]][1] <- 1 ## Alpha value >1 then replace it with 1. Here first five are more than 1 
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
file_name = paste0("ROC_mulivariate_stepAIC", "CpG", ".pdf")    ### change directory name here (the one you created above.
pdf(file_name)
#plot(roc, main = name, col = "red")
plot(roc, main = "Multiple linear regression", colorize=T)  
abline(0, 1, lty = 2)
### Adjust position based on figure (pdf or png)
text(0.80, 0.13, paste0("AUC = ", round(unlist(auc), 2)), cex = 1.4)
text(0.80, 0.08, paste0("95% CI [", ci_chr, "]"), cex = 1.2)
rect(0.65, 0.05, 0.95, 0.16)
dev.off()



#### Without feature selection ########
model.all <- glm(Status ~ cg00033213+cg00435526+cg02502145+cg02511231+cg03449867+cg03667593+cg03893271+cg04657470+cg05111616+cg05238741+cg05273049+cg06459104+cg06960824+cg06968724+cg07326438+cg08024264+cg08732526+cg09084244+cg09279736+cg09968723+cg09980477+cg10306192+cg10701801+cg10909185+cg11035303+cg11251367+cg11377136+cg12466610+cg12615916+cg13569207+cg13665593+cg13767940+cg14345128+cg15829294+cg16178625+cg16884400+cg17369694+cg17386473+cg17707870+cg17770035+cg18208707+cg18478319+cg19393008+cg19496978+cg19863210+cg20459037+cg20973720+cg21435568+cg21514997+cg21532325+cg22152407+cg23165899+cg24072924+cg25046571+cg26311944+cg26349606+cg27226920+cg27392792, data=data)

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

rocs <- performance(prediction(fitted(model.all), data$Status), "tpr", "fpr")
aucs <- performance(prediction(fitted(model.all), data$Status), "auc")
roc <- rocs
auc <- aucs@y.values[[1]]
#attributes(roc)$alpha.values[[1]][1] <- 1
attributes(roc)$alpha.values[[1]][1] <- 1 ## Alpha value >1 then replace it with 1. Here first five are more than 1 
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
file_name = paste0("ROC_mulivariate All", "CpG", ".pdf")    ### change directory name here (the one you created above.
pdf(file_name)
#plot(roc, main = name, col = "red")
plot(roc, main = "Multiple linear regression", colorize=T)  
abline(0, 1, lty = 2)
### Adjust position based on figure (pdf or png)
text(0.80, 0.13, paste0("AUC = ", round(unlist(auc), 2)), cex = 1.4)
text(0.80, 0.08, paste0("95% CI [", ci_chr, "]"), cex = 1.2)
rect(0.65, 0.05, 0.95, 0.16)
dev.off()


performance(prediction(fitted(model.stepAIC), data$Status), "spec", "sens")
performance(prediction(fitted(model.all), data$Status), "spec", "sens")

######################################################
save.image("Multiple_Linear_analysis.RData")
