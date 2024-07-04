library(data.table)
library(ggplot2)
library(scatterplot3d)
setwd("F:/OneDrive - University of Nebraska Medical Center/Radha/NSCLP")
NSCLP <- fread(input = "NSCLPIndividual_data_May 212017 Original newfile.txt")
select <- grep("AVG_Beta", names(NSCLP), value = TRUE)
Meth <- NSCLP[, select, with=FALSE]
#Meth.t <- t(Meth)
Meth<- na.omit(Meth)


colnames(Meth) <- gsub("AVG_", "", colnames(Meth))
#colnames(Meth) <- gsub("NSCLP_n", "N", colnames(Meth))
colnames(Meth) <- gsub(".Beta", "", colnames(Meth))
colnames(Meth) <- gsub("Affected", "Aff", colnames(Meth))
colnames(Meth) <- gsub("normal", "norm", colnames(Meth))
#colnames(Meth) <- gsub("Affected\\d+.","", colnames(Meth))
#colnames(Meth) <- gsub("NCLP|NSCP", "NSCLP", colnames(Meth))


pca <- prcomp(t(Meth))
color <- as.factor(c(rep("NSLP", 14), rep("NCLP", 7), rep("Normal", 21)))

PCi<-data.frame(pca$x,Sample=color)

# s3d <- scatterplot3d(pca$x[, 1], pca$x[, 2], pca$x[, 
#                                                    3], color = pcolor, pch = 19, type = "h", lty.hplot = 0, 
#                      cex.axis = 1.5, cex.lab = 1.5, cex.symbols=1.5, scale.y = 0.75, xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2) * 
#                                                                                                                       100, 1), "%", sep = ""), ylab = paste("Comp 2: ", 
#                                                                                                                                                             round(pca$sdev[2]^2/sum(pca$sdev^2) * 100, 1), 
#                                                                                                                                                             "%", sep = ""), zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2) * 
#                                                                                                                                                                                                              100, 1), "%", sep = ""), angle = angle, ...)
# s3d.coords <- s3d$xyz.convert(pca$x[, 1], pca$x[, 2], 
#                               pca$x[, 3])
# legend("topleft", inset = 0.05, bty = "n", cex = 1.2, title = "Group", 
#        groups, col = c(levels(as.factor(color))), pch = 19)


ggplot(PCi,aes(x=PC1,y=PC2,col=Sample))+
  geom_point(size=5,alpha=0.9)+ #Size and alpha just for fun
  scale_color_manual(values = c("red4","black", "blue4"))+ #your colors here
  theme_classic()
ggsave("PCA.pdf", dpi = 600, width = 5, height = 5)
save.image("PCALplot.RData")
# colors <- c("red4", "blue4", "black")
# colors <- colors[as.numeric(PCi$Sample)]
# scatterplot3d(PCi[,1:3], pch = 16, color=colors)
# 
# 
# 
# s3d <- scatterplot3d(PCi[,1:3], pch = 16, color=colors, angle = 120)
# #legend(s3d$xyz.convert(7.5, 3, 4.5), legend = levels(PCi$Sample),
#        #col =  c("red4", "blue4", "black"), pch = 16)

# plot(pca$x[,1], pca$x[,2], xlab = 'PCA1', ylab = 'PCA2', main = 'PCA plot')
# text(pca$x[,1], pca$x[,2], labels = row.names(pca$x))