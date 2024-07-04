setwd("C:/Users/siddesh.southekal/University of Nebraska Medical Center/Mishra, Nitish K - NAS/HTSeqCount")

load("Group1_vs_Group3.RData")

library(ggplot2)
library(ggrepel)

#2D pca
pca <- prcomp(t(countdata), scale. = TRUE)
pd <- condition
#color <- as.factor(c(rep("FSH-21", 4), rep("FSH-24", 4), rep("Ctrl",4)))
color <- factor(pd, levels = levels(condition))
PCi<-data.frame(pca$x,Sample=color)


ggplot(PCi,aes(x=PC1,y=PC2, col=Sample))+
  geom_point(size=4,alpha=0.9)+ #Size and alpha just for fun
  scale_color_manual(values = c("blue4","red4"))+ #your colors here
  theme_classic() +
  xlab(paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2), 3)*100, "%", sep = "")) +
  ylab(paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""))

ggsave("2D_PCA_Group1_Vs_Group3.pdf", dpi = 600, width = 5, height = 5)



################### 3D PCA plot ######################
#groups <- levels(color)
groups <- levels(condition)
colors <- c(rep("blue4", 32),rep("red4", 32))
pca$pcolor <- colors


s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=1.2, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), angle = 45)
s3d.coords <- s3d$xyz.convert(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100)
legend("topright", inset = 0.02, bty = "n", cex = 1.2, groups, col = c("blue4","red4"), pch = 19)


dev.print(pdf, '3D_PCA_Group1_vs_Group3.pdf', width = 10, height = 10)



########### 3D-PCA plot with sample name #############
s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=0.01, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), angle = 45)

s3d.coords <- s3d$xyz.convert(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80)
legend("topright", inset = 0.05, bty = "n", cex = 1.2, groups, col = c("blue4","red4"), pch = 19)
text(s3d$xyz.convert(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80),labels = rownames(pca$x), cex= 1, col = pca$pcolor)
dev.print(pdf, '3D_PCA_with_Name_Group3_vs_Group1.pdf', width = 12, height = 12)


#######

## Volcano Plot
load("gencode.v33.Rdata")
data.table::setnames(gencode33,"gene_id","Gene")
resdata <- dplyr::left_join(resdata, gencode33, by=c("Gene"))

x.cut=1;y.cut=0.05

Significance <- ifelse(resdata$log2FoldChange >= x.cut & resdata$padj < y.cut, "Upregulated", ifelse(resdata$log2FoldChange <= -x.cut & resdata$padj < y.cut, "Downregulated", "Not significant"))

ggplot(resdata, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significance), size = 3) +
  scale_color_manual(values = c( "green4","gray41" ,"red4")) +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept = -x.cut), colour = "black", linetype = "dashed") + 
  geom_vline(aes(xintercept = x.cut),  colour = "black", linetype = "dashed") + 
  geom_hline(aes(yintercept = -1 * log10(y.cut)), colour = "black", linetype = "dashed") +
  
  geom_text_repel(data = subset(resdata, (padj < 0.027 & abs(log2FoldChange) >=2.6)),
                  aes(label = gene_name)) +
  
  theme(legend.position="right")+
  ggtitle("Volcano plot Group 1 Vs Group 3")+
  xlab("log2(FoldChange)") + ylab("-log10(padj)") +
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5))+
  theme(legend.key = element_rect(colour = "white", fill = "white"))+ ### remove box from legends
  theme(legend.title=element_text(colour="black", size = 10, face = "bold"))+
  theme(legend.text = element_text(colour="black", size = 10, face = "bold"))+
  theme(axis.title = element_text(color="black", face="bold", size=14)) 


ggsave(filename = "VolcanoPlot_Group1_vs_3.pdf", width = 12, height = 8, dpi = 800)






ggplot(resdata, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significance), size = 3) +
  scale_color_manual(values = c( "green4","gray41" ,"red4")) +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept = -x.cut), colour = "black", linetype = "dashed") + 
  geom_vline(aes(xintercept = x.cut),  colour = "black", linetype = "dashed") + 
  geom_hline(aes(yintercept = -1 * log10(y.cut)), colour = "black", linetype = "dashed") +
  
  theme(legend.position="right")+
  ggtitle("Volcano plot Group 1 Vs Group 3")+
  xlab("log2(FoldChange)") + ylab("-log10(padj)") +
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5))+
  theme(legend.key = element_rect(colour = "white", fill = "white"))+ ### remove box from legends
  theme(legend.title=element_text(colour="black", size = 10, face = "bold"))+
  theme(legend.text = element_text(colour="black", size = 10, face = "bold"))+
  theme(axis.title = element_text(color="black", face="bold", size=14)) 


ggsave(filename = "VolcanoPlot_Group1_vs_3_NoName.pdf", width = 12, height = 8, dpi = 800)


