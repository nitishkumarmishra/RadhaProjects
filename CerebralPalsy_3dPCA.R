library(data.table)
library(ggplot2)
setwd("F:/OneDrive - University of Nebraska Medical Center/Radha/Cerebral Palcy/Latest/")
TableInput <- fread(input = "CP_Individual_analysis Feb 17 2017 Original.txt")
select <- grep("AVG_Beta|TargetID", names(TableInput), value = TRUE)
Meth <- TableInput[, select, with=FALSE]
colnames(Meth) <- gsub(".AVG_Beta", "", colnames(Meth))
Meth <- as.data.frame(Meth); Meth <- na.omit(Meth); rownames(Meth) <- Meth$TargetID; Meth$TargetID <- NULL
colnames(Meth) <- gsub("Autism_Control", "CP_Control", colnames(Meth))
colnames(Meth) <- gsub("-B_CP_Control", "_CP_Control", colnames(Meth))


pca <- prcomp(t(Meth), scale. = TRUE)
color <- as.factor(c(rep("Cerebral", 23), rep("Normal", 21)))
PCi<-data.frame(pca$x,Sample=color)

ggplot(PCi,aes(x=PC1,y=PC2, col=Sample))+
  geom_point(size=4,alpha=0.9)+ #Size and alpha just for fun
  scale_color_manual(values = c("blue4","red4"))+ #your colors here
  theme_classic()
ggsave("Cerebral.pdf", dpi = 600, width = 5, height = 5)

################### 3D PCA plot ######################
groups <- levels(color)
colors <- c(rep("red4", 23), rep("blue4", 21))
pca$pcolor <- colors


s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=1.2, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), angle = 45)
s3d.coords <- s3d$xyz.convert(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100)
legend("top", inset = 0.05, bty = "n", cex = 1.2, groups, col = c("blue4", "red4"), pch = 19)
dev.print(pdf, 'Cerebral 3dPCA.pdf', width = 10, height = 10)


########### 3D-PCA plot with sample name #############
s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=0.01, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), angle = 45)

s3d.coords <- s3d$xyz.convert(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100)
legend("top", inset = 0.05, bty = "n", cex = 1.2, groups, col = c("blue4", "red4"), pch = 19)
text(s3d$xyz.convert(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100),labels = rownames(pca$x), cex= 0.5, col = pca$pcolor)
dev.print(pdf, 'Cerebral  With Name.pdf', width = 12, height = 12)

#######################################################
save.image("Cerebral_3dPCA.RData")
