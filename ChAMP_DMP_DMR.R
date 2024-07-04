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

######## DMR analysis by using Bumphunter #########
pd.Beta <- pd[grepl("Group1|Group2", pd)]
Beta <- myNorm.PBC[,grepl("Group1|Group2", pd)]
myDMR.Group1_to_Group2 <- champ.DMR(beta = Beta,pheno=pd.Beta, arraytype = "EPIC", method = "Bumphunter")


pd.Beta <- pd[grepl("Group1|Group3", pd)]
Beta <- myNorm.PBC[,grepl("Group1|Group3", pd)]
myDMR.Group1_to_Group3 <- champ.DMR(beta = Beta,pheno=pd.Beta, arraytype = "EPIC", method = "Bumphunter")


pd.Beta <- pd[grepl("Group2|Group3", pd)]
Beta <- myNorm.PBC[,grepl("Group2|Group3", pd)]
myDMR.Group2_to_Group3 <- champ.DMR(beta = Beta,pheno=pd.Beta, arraytype = "EPIC", method = "Bumphunter")


pd.1 <- gsub("Group3", "Group2", pd)
myDMR.Group1_to_Group2_and_3 <- champ.DMR(beta = myNorm.PBC,pheno=pd.1, arraytype = "EPIC", method = "Bumphunter")

############## GSEA analysis of CpGs ##############
###################################
pd.Beta <- pd[grepl("Group1|Group2", pd)]
Beta <- myNorm.PBC[,grepl("Group1|Group2", pd)]
myGSEA.Group1_to_Group2 <- champ.ebGSEA(beta=Beta,pheno=pd.Beta,arraytype="EPIC")


pd.Beta <- pd[grepl("Group1|Group3", pd)]
Beta <- myNorm.PBC[,grepl("Group1|Group3", pd)]
myGSEA.Group1_to_Group3 <- champ.ebGSEA(beta=Beta,pheno=pd.Beta,arraytype="EPIC")


pd.Beta <- pd[grepl("Group2|Group3", pd)]
Beta <- myNorm.PBC[,grepl("Group2|Group3", pd)]
myGSEA.Group2_to_Group3 <- champ.ebGSEA(beta=myNorm.PBC,pheno=pd.Beta,arraytype="EPIC")


pd.1 <- gsub("Group3", "Group2", pd)
myGSEA.Group1_to_Group2_and_3 <- champ.ebGSEA(beta=Beta,pheno=pd.1,arraytype="EPIC")

##############################################################
##############################################################
save.image("ChAMP_DMP_DMR.RData")
