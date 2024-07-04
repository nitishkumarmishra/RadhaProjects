library(data.table)
library(ggplot2)

setwd("F:/OneDrive - University of Nebraska Medical Center/Radha/NSCLP/")

TableInput <- fread(input = "NSCLPIndividual_data_May 212017 Original newfile.txt")
select <- grep("AVG_Beta", names(TableInput), value = TRUE)
Meth <- TableInput[, select, with=FALSE]
colnames(Meth) <- gsub(".AVG_Beta", "", colnames(Meth))
Meth <- na.omit(Meth)

pca <- prcomp(t(Meth), scale. = TRUE)
color <- as.factor(c(rep("NSCLP", 21), rep("Normal", 21)))
PCi<-data.frame(pca$x,Sample=color)

ggplot(PCi,aes(x=PC1,y=PC2, col=Sample))+
  geom_point(size=4,alpha=0.9)+ #Size and alpha just for fun
  scale_color_manual(values = c("red4","blue4"))+ #your colors here
  theme_classic()
ggsave("PCA NSCLP .pdf", dpi = 600, width = 5, height = 5)




### 3D PCA plot ###########
groups <- levels(color)
colors <- c(rep("red4", 21), rep("blue4", 21))
pca$pcolor <- colors

s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/10, pca$x[, 2]/10, pca$x[, 3]/10, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=1.2, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2) *  100, 1), "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2) * 100, 1), "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2) * 100, 1), "%", sep = ""), angle = 210)

s3d.coords <- s3d$xyz.convert(pca$x[, 1]/10, pca$x[, 2]/10, pca$x[, 3]/10)
legend("top", inset = 0.05, bty = "n", cex = 1.2, title = "Group", groups, col = c("blue4", "red4"), pch = 19)

dev.print(pdf, 'PCA NSCLP 3D.pdf', width = 8, height = 8)

#####################
save.image("NSCLP_PCA_Plot.RData")
