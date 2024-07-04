library(data.table)
library(ggplot2)

setwd("F:/OneDrive - University of Nebraska Medical Center/Radha")

TableInput <- fread(input = "TableControl.txt  DEc 30 analysis al individually.txt")
select <- grep("AVG_Beta", names(TableInput), value = TRUE)
Meth <- TableInput[, select, with=FALSE]
colnames(Meth) <- gsub(".AVG_Beta", "", colnames(Meth))
Meth <- na.omit(Meth)
pca <- prcomp(t(Meth), scale. = TRUE)
color <- as.factor(c(rep("TOF", 23), rep("Normal", 24)))

PCi<-data.frame(pca$x,Sample=color)

ggplot(PCi,aes(x=PC1,y=PC2, col=Sample))+
  geom_point(size=4,alpha=0.9)+ #Size and alpha just for fun
  scale_color_manual(values = c("blue4","red4"))+ #your colors here
  theme_classic()
ggsave("PCA TOF .pdf", dpi = 600, width = 5, height = 5)




### 3D PCA plot ###########
groups <- levels(color)
colors <- c(rep("red4", 23), rep("blue4", 24))
pca$pcolor <- colors
#rownames(pca$x)strsplit(rownames(pca$x),"\\-"),'[[',2
#sapply(strsplit(rownames(pca$x),"\\-"),'[[',2)
s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/10, pca$x[, 2]/10, pca$x[, 3]/10, color = pca$pcolor, pch = 20, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=0.01, scale.y = 0.75, 
                             xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2) *  100, 1), "%", sep = ""), 
                             ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2) * 100, 1), "%", sep = ""), 
                             zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2) * 100, 1), "%", sep = ""), angle = 40)

s3d.coords <- s3d$xyz.convert(pca$x[, 1]/10, pca$x[, 2]/10, pca$x[, 3]/10)
legend("top", inset = 0.05, bty = "n", cex = 1.5, title = "Group", groups, col = c("blue4", "red4"), pch = 20)
text(s3d$xyz.convert(pca$x[, 1]/10, pca$x[, 2]/10, pca$x[, 3]/10),labels = sapply(strsplit(rownames(pca$x),"\\-"),'[[',2), cex= 0.6, col = colors, font = 15)
dev.print(pdf, 'PCA_TOF With Name.pdf', width = 10, height = 10)


#####################
save.image("TOF_PCA.RData")
