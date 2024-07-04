library(data.table)
library(ggplot2)

setwd("F:/OneDrive - University of Nebraska Medical Center/Radha/Alzheimer/")

TableInput <- fread(input = "ALZEH_1.txt")


select <- grep("AVG_Beta", names(TableInput), value = TRUE)
Meth <- TableInput[, select, with=FALSE]
colnames(Meth) <- gsub(".AVG_Beta", "", colnames(Meth))
Meth <- na.omit(Meth)

pca <- prcomp(t(Meth), scale. = TRUE)
color <- as.factor(c(rep("Alzheimer", 24), rep("Normal", 24)))
PCi<-data.frame(pca$x,Sample=color)

ggplot(PCi,aes(x=PC1,y=PC2, col=Sample))+
  geom_point(size=4,alpha=0.9)+ #Size and alpha just for fun
  scale_color_manual(values = c("red4","blue4"))+ #your colors here
  theme_classic()
ggsave("PCA Alzheimer .pdf", dpi = 600, width = 5, height = 5)




### 3D PCA plot ###########
groups <- levels(color)
colors <- c(rep("red4", 24), rep("blue4", 24))
pca$pcolor <- colors

s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/10, pca$x[, 2]/10, pca$x[, 3]/10, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=1.2, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2) *  100, 1), "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2) * 100, 1), "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2) * 100, 1), "%", sep = ""), angle = 30)

s3d.coords <- s3d$xyz.convert(pca$x[, 1]/10, pca$x[, 2]/10, pca$x[, 3]/10)
legend("topright", inset = 0.05, bty = "n", cex = 1.2, groups, col = c("blue4", "red4"), pch = 19)

dev.print(pdf, 'PCA Alzheimer 3D.pdf', width = 8, height = 8)

#####################
rm(TableInput)
save.image("Alzheimer_PCA_Plot.RData")




TableInput <- read.csv(file = "4-Alzheimers_48_Samples_New_523 targets-Individual_June 20 2018_to Nithish.txt", header = TRUE, row.names = 2, sep = "\t")

select <- grep("AVG_Beta", names(TableInput), value = TRUE)
Meth <- TableInput[, select]
colnames(Meth) <- gsub(".AVG_Beta", "", colnames(Meth))
Meth <- na.omit(Meth) ## Now it's 422 line

pca <- prcomp(t(Meth), scale. = TRUE)
color <- as.factor(c(rep("Alzheimer", 24), rep("Normal", 24)))
PCi<-data.frame(pca$x,Sample=color)

ggplot(PCi,aes(x=PC1,y=PC2, col=Sample))+
  geom_point(size=4,alpha=0.9)+ #Size and alpha just for fun
  scale_color_manual(values = c("red4","blue4"))+ #your colors here
  theme_classic()
ggsave("PCA Alzheimer 523 CpG.pdf", dpi = 600, width = 5, height = 5)


### 3D PCA plot ###########
groups <- levels(color)
colors <- c(rep("red4", 24), rep("blue4", 24))
pca$pcolor <- colors

s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/10, pca$x[, 2]/10, pca$x[, 3]/10, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=1.2, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2) *  100, 1), "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2) * 100, 1), "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2) * 100, 1), "%", sep = ""), angle = 60)

s3d.coords <- s3d$xyz.convert(pca$x[, 1]/10, pca$x[, 2]/10, pca$x[, 3]/10)
legend("topright", inset = 0.05, bty = "n", cex = 1.2, title = "Group", groups, col = c("blue4", "red4"), pch = 19)

dev.print(pdf, 'PCA Alzheimer 3D 422 CpG.pdf', width = 8, height = 8)
#######################################
s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/10, pca$x[, 2]/10, pca$x[, 3]/10, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=0.01, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), angle = 55)

s3d.coords <- s3d$xyz.convert(pca$x[, 1]/10, pca$x[, 2]/10, pca$x[, 3]/10)
legend("top", inset = 0.05, bty = "n", cex = 1.2, groups, col = c("blue4", "red4"), pch = 19)
text(s3d$xyz.convert(pca$x[, 1]/10, pca$x[, 2]/10, pca$x[, 3]/10),labels = substr(rownames(pca$x), 10, 17), cex= 0.6, col = pca$pcolor)
dev.print(pdf, 'PCA Alzheimer 3D 422 CpG With Name.pdf', width = 12, height = 12)



########################################
save.image("Alzheimer_PCA_Plot_523CpGs.RData")



library(dplyr)
TableInput <- read.csv(file = "4-Alzheimers_48_Samples_New_523 targets-Individual_June 20 2018_to Nithish.txt", header = TRUE, row.names = 2, sep = "\t")
Radha1 <- read.csv("0-Alzheimers_48-Samples_525 targets-Table 1_June 18 2018_Sent to Nithish.txt", header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
Radha1$delta <- abs(Radha1$MethylationCases - Radha1$MethylationControl)/100
tmp <- TableInput[, c("UCSC_REFGENE_GROUP", "ILMNID")]
tmp1 <- merge(Radha1, tmp, by = "row.names")
rownames(tmp1) <- tmp1$Row.names
tmp1 <- tmp1[rownames(Meth),]
annotation <- tmp1$UCSC_REFGENE_GROUP
annotation <- as.character(annotation)
annotation[which(nchar(annotation)<2)]<-NA ### When no attotation
annotation <- strsplit(annotation,"\\;")
annotation <- sapply(annotation, "[[", 1)
logFC <- log2(tmp1$`Fold change`)
direction <- ifelse( logFC>0,"Hyper","Hypo")
#logP <- -log10(tmp1$FDR.p.Val)
delta <- abs(tmp1$MethylationCases-tmp1$MethylationControl)/100
AUC <- tmp1$AUC
Beta = t(scale(t(Meth)))

library(ComplexHeatmap)
library(circlize)
set.seed(121)
ha = HeatmapAnnotation(df = data.frame(Samples = c(rep("Alzheimer", 24), rep("Control", 24))), col = list(Samples = c("Alzheimer" = "red4", "Control"="blue4")))
#column_tree = hclust(dist(t(meth)), method = "ward.D2")
ht_list = Heatmap(Beta, clustering_method_columns = "ward.D2", name = "Methylation", col = colorRamp2(c(-3.0, -2,-1, 0, 1, 2, 4,6), c("green4", "green", "aquamarine","white","red1","red2","firebrick","red4")),
                  cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 4, show_row_names = FALSE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(AUC, name = "AUC", col = colorRamp2(c(0.6, 0.8, 1.0), c("blue4", "purple", "red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(logFC, name = "log2(FC)", col = colorRamp2(c(2, 1, 0, -1, -2), c("red4","red","purple", "blue", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(direction, name = "Direction", col = c("Hyper" = "red4", "Hypo" = "blue4"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(delta, name = "Delta Beta", col = colorRamp2(c(0.25, 0.12, 0), c("red4","red", "blue4")), column_names_gp = gpar(fontsize = 8))+
  #Heatmap(logP, name = "-log10 (FDR)", col = colorRamp2(c(0, 10, 20), c("red4","red", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(annotation, name = "Annotation", col = c("Body" = "aquamarine", "TSS200" = "red4", "TSS1500" = "violetred", "1stExon" = "black", "5'UTR" = "blue", "3'UTR" = "green", ExonBnd= "purple"), column_names_gp = gpar(fontsize = 8))
#Heatmap(relation, name = "Probe Relationship", column_names_gp = gpar(fontsize = 8))
pdf("Alzheimer.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

#######################################################################
library(dplyr)
Radha1 <- read.csv("0-Alzheimers_48-Samples_525 targets-Table 1_June 18 2018_Sent to Nithish.txt", header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
Radha1$delta <- abs(Radha1$MethylationCases - Radha1$MethylationControl)/100
tmp <- TableInput[, c("UCSC_REFGENE_GROUP", "ILMNID")]
tmp1 <- merge(Radha1, tmp, by = "row.names")
rownames(tmp1) <- tmp1$Row.names
tmp1 <- tmp1[rownames(Meth),]
tmp1 <- tmp1[order(tmp1$delta, decreasing = TRUE),]
tmp2 <- tmp1[1:250,]

annotation <- tmp2$UCSC_REFGENE_GROUP
annotation <- as.character(annotation)
annotation[which(nchar(annotation)<2)]<-NA ### When no attotation
annotation <- strsplit(annotation,"\\;")
annotation <- sapply(annotation, "[[", 1)
logFC <- log2(tmp2$`Fold change`)
direction <- ifelse( logFC>0,"Hyper","Hypo")
#logP <- -log10(tmp1$FDR.p.Val)
delta <- tmp2$delta
AUC <- tmp2$AUC
Beta <- Meth[rownames(tmp2),]
Beta = t(scale(t(Beta)))
colnames(Beta) <- gsub("Alzheimers.", "", colnames(Beta))

library(ComplexHeatmap)
library(circlize)
set.seed(121)
ha = HeatmapAnnotation(df = data.frame(Samples = c(rep("Alzheimer", 24), rep("Control", 24))), col = list(Samples = c("Alzheimer" = "red4", "Control"="blue4")))
#column_tree = hclust(dist(t(meth)), method = "ward.D2")
ht_list = Heatmap(Beta, clustering_method_columns = "ward.D2", name = "Methylation", col = colorRamp2(c(-4.0, -2,-1.0, 0,1, 2, 4, 6), c("green4","green3","aquamarine","white","firebrick","red2","red3","red4")),
                  cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8, fontface = "bold"), km = 5, show_row_names = FALSE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(AUC, name = "AUC", col = colorRamp2(c(0.6, 0.8, 1.0), c("blue4", "purple", "red4")), column_names_gp = gpar(fontsize = 8, fontface = "bold"))+
  Heatmap(logFC, name = "log2(FC)", col = colorRamp2(c(2, 1, 0, -1, -2), c("red4","red","purple", "blue", "blue4")), column_names_gp = gpar(fontsize = 8, fontface = "bold"))+
  Heatmap(direction, name = "Direction", col = c("Hyper" = "red4", "Hypo" = "blue4"), column_names_gp = gpar(fontsize = 8, fontface = "bold"))+
  Heatmap(delta, name = "Delta Beta", col = colorRamp2(c(0.25, 0.12, 0), c("red4","red", "blue4")), column_names_gp = gpar(fontsize = 8, fontface = "bold"))+
  #Heatmap(logP, name = "-log10 (FDR)", col = colorRamp2(c(0, 10, 20), c("red4","red", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(annotation, name = "Annotation", col = c("Body" = "aquamarine", "TSS200" = "red4", "TSS1500" = "violetred", "1stExon" = "black", "5'UTR" = "blue", "3'UTR" = "green", ExonBnd= "purple"), column_names_gp = gpar(fontsize = 8, fontface = "bold"))
#Heatmap(relation, name = "Probe Relationship", column_names_gp = gpar(fontsize = 8))
pdf("Alzheimer Top 250 CpGs .pdf", width = 8, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 10, fontface = "bold"), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
save.image("Alzheimer_Heatmap.RData")
