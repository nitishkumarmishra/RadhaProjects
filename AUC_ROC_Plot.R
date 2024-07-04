library(ChAMP)
library(dplyr)
library(tibble)

setwd("D:/OneDrive - University of Nebraska Medical Center/NAS/")
load("D:/OneDrive - University of Nebraska Medical Center/NAS/NAS_withXY.RData")

########### DMP analysis by using ChAMP ############
pd <- ifelse(grepl("noOPIOD", colnames(myNorm.PBC)), "Group3", ifelse(grepl("OPIOD_No_NAS", colnames(myNorm.PBC)), "Group2", "Group1"))
myDMP <- champ.DMP(beta = myNorm.PBC,pheno=pd, arraytype = "EPIC", adjPVal = 0.05)

myDMP$Group3_to_Group1 <- myDMP$Group3_to_Group1 %>%
  rownames_to_column("rowname") %>%
  mutate(deltaBeta=logFC,FoldChange=Group3_AVG/Group1_AVG) %>%
  mutate(logFC=log2(FoldChange)) %>%
  column_to_rownames("rowname")

myDMP$Group3_to_Group2 <- myDMP$Group3_to_Group2 %>%
  rownames_to_column("rowname") %>%
  mutate(deltaBeta=logFC,FoldChange=Group3_AVG/Group2_AVG) %>%
  mutate(logFC=log2(FoldChange)) %>%
  column_to_rownames("rowname")



myDMP$Group1_to_Group2 <- myDMP$Group1_to_Group2 %>%
  rownames_to_column("rowname") %>%
  mutate(deltaBeta=logFC,deltaBeta=deltaBeta*-1) %>%
  mutate(FoldChange=Group1_AVG/Group2_AVG) %>%
  mutate(logFC=log2(FoldChange)) %>%
  column_to_rownames("rowname")


pd.1 <- gsub("Group3", "Group2", pd)
myDMP.1 <- champ.DMP(beta = myNorm.PBC,pheno=pd.1, arraytype = "EPIC", adjPVal = 0.05)
myDMP.1$Group2_to_Group1 <- myDMP.1$Group2_to_Group1 %>%
  rownames_to_column("rowname") %>%
  mutate(FoldChange=Group1_AVG/Group2_AVG) %>%
  mutate(logFC=log2(FoldChange)) %>%
  column_to_rownames("rowname")

pd.2 <- gsub("Group2", "Group1", pd)
myDMP.2 <- champ.DMP(beta = myNorm.PBC,pheno=pd.2, arraytype = "EPIC", adjPVal = 0.05)
myDMP.2$Group3_to_Group1 <- myDMP.2$Group3_to_Group1 %>%
  rownames_to_column("rowname") %>%
  mutate(FoldChange=Group1_AVG/Group3_AVG) %>%
  mutate(logFC=log2(FoldChange)) %>%
  column_to_rownames("rowname")
####################################################################
########################## roc calculation #########################
####################################################################
library(dplyr)
library(reshape2)
library(ROCR)


m1 <-myNorm.PBC[rownames(myDMP$Group3_to_Group1), grep("Group2", invert = TRUE, pd)]

ind_Exp <- as.data.frame(m1)
ind_Exp$Index <- rownames(ind_Exp)
ind_Exp$TargetID <- rownames(ind_Exp)
ind_Exp <- ind_Exp %>% filter(!is.na(Index)) %>% droplevels()
rownames(ind_Exp) <- ind_Exp$TargetID
#ind_Exp <- ind_Exp[rownames(diffMeth),]
dat <- ind_Exp %>% dplyr::select(Index, TargetID, contains("NAS")) %>% melt(id = c("Index", "TargetID")) %>% mutate(group = ifelse(grepl("-Aff", variable), 1, 0)) %>% split(ind_Exp$TargetID)
#### Now code for R function ########### 
#It will write two functio se_auc ci_auc
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

############ R code for AUC ###########
aucs <- lapply(dat, function(x) performance(prediction(x[, "value"], x[, "group"]), "auc")@y.values)
rocs <- lapply(dat, function(x) performance(prediction(x[, "value"], x[, "group"]), "tpr", "fpr"))
# Create directory (if it doesn't exist)
if (!dir.exists("DiffMeth_ROC_Plot")) dir.create("DiffMeth_ROC_Plot")   ### change directory name here
# Loop through ROCs to make plots.
ns <- unname(table(dat[[1]]$group))
auc_df <- data.frame(name = character(length(dat)), 
                     AUC = numeric(length(dat)), 
                     CI_lower = numeric(length(dat)), 
                     CI_upper = numeric(length(dat)), stringsAsFactors = F)
for (i in 1:length(rocs)) {
  roc <- rocs[[i]]
  auc <- aucs[[i]]
  name <- names(aucs[i])
  attributes(roc)$alpha.values[[1]][1] <- 1
  # Flip results if AUC < 0.5
  ## Why he is flipping AUC values
  if (auc < 0.5) {
    auc <- 1 - unlist(auc)
    temp <- roc@x.values
    roc@x.values <- roc@y.values
    roc@y.values <- temp
  }
  ci <- ci_auc(unlist(auc), ns[1], ns[2])
  ci_chr <- paste0(round(ci, 2), collapse = ", ") 
  auc_df[i, "name"] <- name
  auc_df[i, c("AUC", "CI_lower", "CI_upper")] <- c(auc, ci[1], ci[2])
  if (auc > 0.5) {
    file_name = paste0("DiffMeth_ROC_Plot/roc_", name, ".png")    ### change directory name here (the one you created above.
    png(file_name)
    #plot(roc, main = name, col = "red")
    plot(roc, main = name, colorize=T)  
    abline(0, 1, lty = 2)
    ## We can change position by changing 0.65, 0.20 etc.
    ## These range are based on ROC plot X and Y axis value [1,1] range
    ## AUC start from X aix=0.8, Y axis = 0.13, similarly CI start from position 0.8 on X and 0.08 on Y
    ## rect(0.62, 0.03, 0.98, 0.18); It start on X=0.62 end X=0.98 (go from 0.66 to 0.98 on X)
    ## On Y axis start Y= 0.03 go to 0.18
    text(0.80, 0.13, paste0("AUC = ", round(unlist(auc), 2)), cex = 1.4)
    text(0.80, 0.08, paste0("95% CI [", ci_chr, "]"), cex = 1.2)
    rect(0.62, 0.03, 0.98, 0.18)
    dev.off()
  }
}
auc_df_Group3_to_Group1 <- arrange(auc_df, desc(AUC)) %>% filter(AUC > 0.05)
write.csv(auc_df_Group3_to_Group1, "AUC_table_DiffMeth_ROC_Plot_Group3_to_Group1.csv", row.names = F)   ### change filename to save AUC values.
write.csv(auc_df_Group3_to_Group1, "AUC_table_DiffMeth_ROC_Plot_Group3_to_Group1.txt", row.names = F)   ### change filename to save AUC values.

#################################################################
#################################################################

m2 <-myNorm.PBC[rownames(myDMP$Group3_to_Group2), grep("Group1", invert = TRUE, pd)]

ind_Exp <- as.data.frame(m2)
ind_Exp$Index <- rownames(ind_Exp)
ind_Exp$TargetID <- rownames(ind_Exp)
ind_Exp <- ind_Exp %>% filter(!is.na(Index)) %>% droplevels()
rownames(ind_Exp) <- ind_Exp$TargetID
#ind_Exp <- ind_Exp[rownames(diffMeth),]
dat <- ind_Exp %>% dplyr::select(Index, TargetID, contains("NAS")) %>% melt(id = c("Index", "TargetID")) %>% mutate(group = ifelse(grepl("noOPIOD", variable), 0, 1)) %>% split(ind_Exp$TargetID)
#### Now code for R function ########### 
#It will write two functio se_auc ci_auc
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

############ R code for AUC ###########
aucs <- lapply(dat, function(x) performance(prediction(x[, "value"], x[, "group"]), "auc")@y.values)
rocs <- lapply(dat, function(x) performance(prediction(x[, "value"], x[, "group"]), "tpr", "fpr"))
# Create directory (if it doesn't exist)
if (!dir.exists("DiffMeth_ROC_Plot_Group3_to_Group2")) dir.create("DiffMeth_ROC_Plot_Group3_to_Group2")   ### change directory name here
# Loop through ROCs to make plots.
ns <- unname(table(dat[[1]]$group))
auc_df <- data.frame(name = character(length(dat)), 
                     AUC = numeric(length(dat)), 
                     CI_lower = numeric(length(dat)), 
                     CI_upper = numeric(length(dat)), stringsAsFactors = F)
for (i in 1:length(rocs)) {
  roc <- rocs[[i]]
  auc <- aucs[[i]]
  name <- names(aucs[i])
  attributes(roc)$alpha.values[[1]][1] <- 1
  # Flip results if AUC < 0.5
  ## Why he is flipping AUC values
  if (auc < 0.5) {
    auc <- 1 - unlist(auc)
    temp <- roc@x.values
    roc@x.values <- roc@y.values
    roc@y.values <- temp
  }
  ci <- ci_auc(unlist(auc), ns[1], ns[2])
  ci_chr <- paste0(round(ci, 2), collapse = ", ") 
  auc_df[i, "name"] <- name
  auc_df[i, c("AUC", "CI_lower", "CI_upper")] <- c(auc, ci[1], ci[2])
  if (auc > 0.5) {
    file_name = paste0("DiffMeth_ROC_Plot_Group3_to_Group2/roc_", name, ".png")    ### change directory name here (the one you created above.
    png(file_name)
    #plot(roc, main = name, col = "red")
    plot(roc, main = name, colorize=T)  
    abline(0, 1, lty = 2)
    ## We can change position by changing 0.65, 0.20 etc.
    ## These range are based on ROC plot X and Y axis value [1,1] range
    ## AUC start from X aix=0.8, Y axis = 0.13, similarly CI start from position 0.8 on X and 0.08 on Y
    ## rect(0.62, 0.03, 0.98, 0.18); It start on X=0.62 end X=0.98 (go from 0.66 to 0.98 on X)
    ## On Y axis start Y= 0.03 go to 0.18
    text(0.80, 0.13, paste0("AUC = ", round(unlist(auc), 2)), cex = 1.4)
    text(0.80, 0.08, paste0("95% CI [", ci_chr, "]"), cex = 1.2)
    rect(0.62, 0.03, 0.98, 0.18)
    dev.off()
  }
}
auc_df_Group3_to_Group2 <- arrange(auc_df, desc(AUC)) %>% filter(AUC > 0.05)
write.csv(auc_df_Group3_to_Group2, "AUC_table_DiffMeth_ROC_Plot_Group3_to_Group2.csv", row.names = F)   ### change filename to save AUC values.
write.csv(auc_df_Group3_to_Group2, "AUC_table_DiffMeth_ROC_Plot_Group3_to_Group2.txt", row.names = F)   ### change filename to save AUC values.

#################################################################
#################################################################

m3 <-myNorm.PBC[rownames(myDMP$Group1_to_Group2), grep("Group3", invert = TRUE, pd)]

ind_Exp <- as.data.frame(m3)
ind_Exp$Index <- rownames(ind_Exp)
ind_Exp$TargetID <- rownames(ind_Exp)
ind_Exp <- ind_Exp %>% filter(!is.na(Index)) %>% droplevels()
rownames(ind_Exp) <- ind_Exp$TargetID
#ind_Exp <- ind_Exp[rownames(diffMeth),]
dat <- ind_Exp %>% dplyr::select(Index, TargetID, contains("NAS")) %>% melt(id = c("Index", "TargetID")) %>% mutate(group = ifelse(grepl("OPIOD_NAS-Aff", variable), 1, 0)) %>% split(ind_Exp$TargetID)
#### Now code for R function ########### 
#It will write two functio se_auc ci_auc
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

############ R code for AUC ###########
aucs <- lapply(dat, function(x) performance(prediction(x[, "value"], x[, "group"]), "auc")@y.values)
rocs <- lapply(dat, function(x) performance(prediction(x[, "value"], x[, "group"]), "tpr", "fpr"))
# Create directory (if it doesn't exist)
if (!dir.exists("DiffMeth_ROC_Plot_Group1_to_Group2")) dir.create("DiffMeth_ROC_Plot_Group1_to_Group2")   ### change directory name here
# Loop through ROCs to make plots.
ns <- unname(table(dat[[1]]$group))
auc_df <- data.frame(name = character(length(dat)), 
                     AUC = numeric(length(dat)), 
                     CI_lower = numeric(length(dat)), 
                     CI_upper = numeric(length(dat)), stringsAsFactors = F)
for (i in 1:length(rocs)) {
  roc <- rocs[[i]]
  auc <- aucs[[i]]
  name <- names(aucs[i])
  attributes(roc)$alpha.values[[1]][1] <- 1
  # Flip results if AUC < 0.5
  ## Why he is flipping AUC values
  if (auc < 0.5) {
    auc <- 1 - unlist(auc)
    temp <- roc@x.values
    roc@x.values <- roc@y.values
    roc@y.values <- temp
  }
  ci <- ci_auc(unlist(auc), ns[1], ns[2])
  ci_chr <- paste0(round(ci, 2), collapse = ", ") 
  auc_df[i, "name"] <- name
  auc_df[i, c("AUC", "CI_lower", "CI_upper")] <- c(auc, ci[1], ci[2])
  if (auc > 0.5) {
    file_name = paste0("DiffMeth_ROC_Plot_Group1_to_Group2/roc_", name, ".png")    ### change directory name here (the one you created above.
    png(file_name)
    #plot(roc, main = name, col = "red")
    plot(roc, main = name, colorize=T)  
    abline(0, 1, lty = 2)
    text(0.80, 0.13, paste0("AUC = ", round(unlist(auc), 2)), cex = 1.4)
    text(0.80, 0.08, paste0("95% CI [", ci_chr, "]"), cex = 1.2)
    rect(0.62, 0.03, 0.98, 0.18)
    dev.off()
  }
}
auc_df_Group1_to_Group2 <- arrange(auc_df, desc(AUC)) %>% filter(AUC > 0.05)
write.csv(auc_df_Group1_to_Group2, "AUC_table_DiffMeth_ROC_Plot_Group1_to_Group2.csv", row.names = F)   ### change filename to save AUC values.
write.csv(auc_df_Group1_to_Group2, "AUC_table_DiffMeth_ROC_Plot_Group1_to_Group2.txt", row.names = F)   ### change filename to save AUC values.

#################################################################
#################################################################
m4 <-myNorm.PBC[rownames(myDMP.1$Group2_to_Group1), grep("Group3", invert = TRUE, pd)]
ind_Exp <- as.data.frame(m4)
ind_Exp$Index <- rownames(ind_Exp)
ind_Exp$TargetID <- rownames(ind_Exp)
ind_Exp <- ind_Exp %>% filter(!is.na(Index)) %>% droplevels()
rownames(ind_Exp) <- ind_Exp$TargetID
#ind_Exp <- ind_Exp[rownames(diffMeth),]
dat <- ind_Exp %>% dplyr::select(Index, TargetID, contains("NAS")) %>% melt(id = c("Index", "TargetID")) %>% mutate(group = ifelse(grepl("OPIOD_NAS-Aff", variable), 1, 0)) %>% split(ind_Exp$TargetID)
#### Now code for R function ########### 
#It will write two functio se_auc ci_auc
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

############ R code for AUC ###########
aucs <- lapply(dat, function(x) performance(prediction(x[, "value"], x[, "group"]), "auc")@y.values)
rocs <- lapply(dat, function(x) performance(prediction(x[, "value"], x[, "group"]), "tpr", "fpr"))
# Create directory (if it doesn't exist)
if (!dir.exists("DiffMeth_ROC_Plot_Group2_to_Group1")) dir.create("DiffMeth_ROC_Plot_Group2_to_Group1")   ### change directory name here
# Loop through ROCs to make plots.
ns <- unname(table(dat[[1]]$group))
auc_df <- data.frame(name = character(length(dat)), 
                     AUC = numeric(length(dat)), 
                     CI_lower = numeric(length(dat)), 
                     CI_upper = numeric(length(dat)), stringsAsFactors = F)
for (i in 1:length(rocs)) {
  roc <- rocs[[i]]
  auc <- aucs[[i]]
  name <- names(aucs[i])
  attributes(roc)$alpha.values[[1]][1] <- 1
  # Flip results if AUC < 0.5
  ## Why he is flipping AUC values
  if (auc < 0.5) {
    auc <- 1 - unlist(auc)
    temp <- roc@x.values
    roc@x.values <- roc@y.values
    roc@y.values <- temp
  }
  ci <- ci_auc(unlist(auc), ns[1], ns[2])
  ci_chr <- paste0(round(ci, 2), collapse = ", ") 
  auc_df[i, "name"] <- name
  auc_df[i, c("AUC", "CI_lower", "CI_upper")] <- c(auc, ci[1], ci[2])
  if (auc > 0.5) {
    file_name = paste0("DiffMeth_ROC_Plot_Group2_to_Group1/roc_", name, ".png")    ### change directory name here (the one you created above.
    png(file_name)
    #plot(roc, main = name, col = "red")
    plot(roc, main = name, colorize=T)  
    abline(0, 1, lty = 2)
    text(0.80, 0.13, paste0("AUC = ", round(unlist(auc), 2)), cex = 1.4)
    text(0.80, 0.08, paste0("95% CI [", ci_chr, "]"), cex = 1.2)
    rect(0.62, 0.03, 0.98, 0.18)
    dev.off()
  }
}
auc_df_Group2_to_Group1 <- arrange(auc_df, desc(AUC)) %>% filter(AUC > 0.05)
write.csv(auc_df_Group2_to_Group1, "AUC_table_DiffMeth_ROC_Plot_Group2_to_Group1.csv", row.names = F)   ### change filename to save AUC values.
write.csv(auc_df_Group2_to_Group1, "AUC_table_DiffMeth_ROC_Plot_Group2_to_Group1.txt", row.names = F)   ### change filename to save AUC values.

#################################################################
#################################################################
## Group1+Group2 Vs. Group3
#pd.2 <- gsub("Group2", "Group1", pd)

m5 <-myNorm.PBC[rownames(myDMP.2$Group3_to_Group1), ]
ind_Exp <- as.data.frame(m5)
ind_Exp$Index <- rownames(ind_Exp)
ind_Exp$TargetID <- rownames(ind_Exp)
ind_Exp <- ind_Exp %>% filter(!is.na(Index)) %>% droplevels()
rownames(ind_Exp) <- ind_Exp$TargetID
#ind_Exp <- ind_Exp[rownames(diffMeth),]
dat <- ind_Exp %>% dplyr::select(Index, TargetID, contains("NAS")) %>% melt(id = c("Index", "TargetID")) %>% mutate(group = ifelse(grepl("noOPIOD", variable), 0, 1)) %>% split(ind_Exp$TargetID)
#### Now code for R function ########### 
#It will write two functio se_auc ci_auc
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

############ R code for AUC ###########
aucs <- lapply(dat, function(x) performance(prediction(x[, "value"], x[, "group"]), "auc")@y.values)
rocs <- lapply(dat, function(x) performance(prediction(x[, "value"], x[, "group"]), "tpr", "fpr"))
# Create directory (if it doesn't exist)
if (!dir.exists("DiffMeth_ROC_Plot_Group1_and_2_to_Group3")) dir.create("DiffMeth_ROC_Plot_Group1_and_2_to_Group3")   ### change directory name here
# Loop through ROCs to make plots.
ns <- unname(table(dat[[1]]$group))
auc_df <- data.frame(name = character(length(dat)), 
                     AUC = numeric(length(dat)), 
                     CI_lower = numeric(length(dat)), 
                     CI_upper = numeric(length(dat)), stringsAsFactors = F)
for (i in 1:length(rocs)) {
  roc <- rocs[[i]]
  auc <- aucs[[i]]
  name <- names(aucs[i])
  attributes(roc)$alpha.values[[1]][1] <- 1
  # Flip results if AUC < 0.5
  ## Why he is flipping AUC values
  if (auc < 0.5) {
    auc <- 1 - unlist(auc)
    temp <- roc@x.values
    roc@x.values <- roc@y.values
    roc@y.values <- temp
  }
  ci <- ci_auc(unlist(auc), ns[1], ns[2])
  ci_chr <- paste0(round(ci, 2), collapse = ", ") 
  auc_df[i, "name"] <- name
  auc_df[i, c("AUC", "CI_lower", "CI_upper")] <- c(auc, ci[1], ci[2])
  if (auc > 0.5) {
    file_name = paste0("DiffMeth_ROC_Plot_Group1_and_2_to_Group3/roc_", name, ".png")    ### change directory name here (the one you created above.
    png(file_name)
    #plot(roc, main = name, col = "red")
    plot(roc, main = name, colorize=T)  
    abline(0, 1, lty = 2)
    text(0.80, 0.13, paste0("AUC = ", round(unlist(auc), 2)), cex = 1.4)
    text(0.80, 0.08, paste0("95% CI [", ci_chr, "]"), cex = 1.2)
    rect(0.62, 0.03, 0.98, 0.18)
    dev.off()
  }
}
auc_df_Group1_and_2_to_Group3 <- arrange(auc_df, desc(AUC)) %>% filter(AUC > 0.05)
write.csv(auc_df_Group1_and_2_to_Group3, "AUC_table_DiffMeth_ROC_Plot_Group1_and_2_to_Group3.csv", row.names = F)   ### change filename to save AUC values.
write.csv(auc_df_Group1_and_2_to_Group3, "AUC_table_DiffMeth_ROC_Plot_Group1_and_2_to_Group3.txt", row.names = F)   ### change filename to save AUC values.

#################################################################
######### Merge the AUC and DMP file and save it ################
Group3_to_Group1 <- merge(myDMP$Group3_to_Group1, auc_df_Group3_to_Group1, by.x="row.names", by.y="name")
rownames(Group3_to_Group1) <- Group3_to_Group1$Row.names; Group3_to_Group1$Row.names=NULL
write.csv(Group3_to_Group1, "DiffMeth_AUC_Group3_to_Group1.csv",row.names = TRUE)


Group3_to_Group2 <- merge(myDMP$Group3_to_Group2, auc_df_Group3_to_Group2, by.x="row.names", by.y="name")
rownames(Group3_to_Group2) <- Group3_to_Group2$Row.names; Group3_to_Group2$Row.names=NULL
write.csv(Group3_to_Group2, "DiffMeth_AUC_Group3_to_Group2.csv",row.names = TRUE)

Group1_to_Group2 <- merge(myDMP$Group1_to_Group2, auc_df_Group1_to_Group2, by.x="row.names", by.y="name")
rownames(Group1_to_Group2) <- Group1_to_Group2$Row.names; Group1_to_Group2$Row.names=NULL
write.csv(Group1_to_Group2, "DiffMeth_AUC_Group1_to_Group2.csv",row.names = TRUE)

Group2_to_Group1 <- merge(myDMP.1$Group2_to_Group1, auc_df_Group2_to_Group1, by.x="row.names", by.y="name")
rownames(Group2_to_Group1) <- Group2_to_Group1$Row.names; Group2_to_Group1$Row.names=NULL
write.csv(Group2_to_Group1, "DiffMeth_AUC_Group2_to_Group1.csv",row.names = TRUE)

Group1_and_2_to_Group3 <- merge(myDMP.2$Group3_to_Group1, auc_df_Group1_and_2_to_Group3, by.x="row.names", by.y="name")
rownames(Group1_and_2_to_Group3) <- Group1_and_2_to_Group3$Row.names; Group1_and_2_to_Group3$Row.names=NULL
write.csv(Group1_and_2_to_Group3, "DiffMeth_AUC_Group1_and_2_to_Group3.csv",row.names = TRUE)

#################################################################
#################################################################
save.image("DiffMethAUC_ROC_Plot.RData")
#################################################################