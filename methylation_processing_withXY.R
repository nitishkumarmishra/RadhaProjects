setwd("C:/Users/Siddesh.Southekal/OneDrive - University of Nebraska Medical Center/nitish_radhakrishna/NAS")

library(ChAMP)
library(doParallel)


myLoad <- champ.load(directory = getwd(),
                     method="ChAMP",
                     methValue="B", # Beta values
                     autoimpute=TRUE,
                     filterDetP=TRUE,
                     ProbeCutoff=0.25,
                     SampleCutoff=0.25,
                     detPcut=0.05,
                     filterNoCG=TRUE,
                     arraytype="EPIC",
                     filterXY = FALSE)

setwd("C:/Users/Siddesh.Southekal/OneDrive - University of Nebraska Medical Center/nitish_radhakrishna/NAS")


myNorm.BMIQ <- champ.norm(arraytype="EPIC",plotBMIQ=TRUE, resultsDir="./BMIQ_withXY/")  # normalization

setwd("C:/Users/Siddesh.Southekal/OneDrive - University of Nebraska Medical Center/nitish_radhakrishna/NAS")

myNorm.PBC <- champ.norm(arraytype="EPIC",method = "PBC", resultsDir="./PBC_withXY/")  # normalization

setwd("C:/Users/Siddesh.Southekal/OneDrive - University of Nebraska Medical Center/nitish_radhakrishna/NAS")


save.image("NAS_withXY.RData")






library(data.table)
library(ggplot2)
#Meth <- as.data.frame(Meth); Meth <- na.omit(Meth); rownames(Meth) <- Meth$TargetID; Meth$TargetID <- NULL
pca <- prcomp(t(myNorm.PBC), scale. = TRUE)
pd <- ifelse(grepl("noOPIOD", colnames(myNorm.PBC)), "noOPIOD", ifelse(grepl("OPIOD_No_NAS", colnames(myNorm.PBC)), "OPIOD_No_NAS", "OPIOD_NAS"))
#color <- as.factor(c(rep("noOPIOD", 32),rep("OPIOD_NAS", 32),rep("OPIOD_No_NAS", 32)))
color <- as.factor(pd)
PCi<-data.frame(pca$x,Sample=color)


ggplot(PCi,aes(x=PC1,y=PC2, col=Sample))+
  geom_point(size=4,alpha=0.9)+ #Size and alpha just for fun
  scale_color_manual(values = c("blue4","red4","green4"))+ #your colors here
  theme_classic()

ggsave("2D_PCA_withXY.pdf", dpi = 600, width = 5, height = 5)

################### 3D PCA plot ######################
groups <- levels(color)
#group.color <- c("noOPIOD"="blue4","OPIOD_NAS"="red4","OPIOD_No_NAS"="green4")
#colors <- c(rep("red4", 23), rep("blue4", 24))
#pca$pcolor <- colors


pca$pcolor <- pd

pca$pcolor <- ifelse(grepl("noOPIOD", pca$pcolor), "blue4", ifelse(grepl("OPIOD_No_NAS", pca$pcolor), "green4", "red4"))

s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=1.2, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), angle = 45)
s3d.coords <- s3d$xyz.convert(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100)
legend("top", inset = 0.05, bty = "n", cex = 1.2, groups, col = c("blue4", "red4","green4"), pch = 19)
dev.print(pdf, '3D_PCA_withXY.pdf', width = 10, height = 10)



########### 3D-PCA plot with sample name #############
s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=0.01, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), angle = 45)

s3d.coords <- s3d$xyz.convert(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80)
legend("top", inset = 0.05, bty = "n", cex = 1.2, groups, col = c("blue4", "red4","green4"), pch = 19)
text(s3d$xyz.convert(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80),labels = rownames(pca$x), cex= 0.6, col = pca$pcolor)
dev.print(pdf, '3D_PCA_withXY_withName.pdf', width = 12, height = 12)
