library(data.table)
library(ggplot2)

setwd("F:/OneDrive - University of Nebraska Medical Center/Radha/HLH/")

TableInput <- fread(input = "1-HLHS Individual data Oct 31 2017 Original.txt")
select <- grep("AVG_Beta", names(TableInput), value = TRUE)
Meth <- TableInput[, select, with=FALSE]
Meth <- as.data.frame(Meth)
rownames(Meth) <- TableInput$TargetID
colnames(Meth) <- gsub(".AVG_Beta", "", colnames(Meth))

Meth <- na.omit(Meth)

######### PCA Plot #########
pca <- prcomp(t(Meth), scale. = TRUE)
color <- as.factor(c(rep("HLHS", 23), rep("Normal", 24)))
PCi<-data.frame(pca$x,Sample=color)

ggplot(PCi,aes(x=PC1,y=PC2, col=Sample))+
  geom_point(size=4,alpha=0.9)+ #Size and alpha just for fun
  scale_color_manual(values = c("red4","blue4"))+ #your colors here
  theme_classic()
ggsave("PCA HLHS .pdf", dpi = 600, width = 5, height = 5)




### 3D PCA plot ###########
groups <- levels(color)
colors <- c(rep("red4", 23), rep("blue4", 24))
pca$pcolor <- colors

s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/10, pca$x[, 2]/10, pca$x[, 3]/10, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=1.2, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2) *  100, 1), "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2) * 100, 1), "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2) * 100, 1), "%", sep = ""), angle = 40)

s3d.coords <- s3d$xyz.convert(pca$x[, 1]/10, pca$x[, 2]/10, pca$x[, 3]/10)
legend("bottom", inset = 0.05, bty = "n", cex = 1.2, title = "Group", groups, col = c("blue4", "red4"), pch = 19)

dev.print(pdf, 'PCA HLHS 3D.pdf', width = 8, height = 8)

save.image("HLHS_PCA_Plot.RData")

#########################################################################################


