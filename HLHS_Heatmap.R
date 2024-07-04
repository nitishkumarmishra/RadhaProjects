library(ComplexHeatmap)
library(circlize)
setwd("F:/OneDrive - University of Nebraska Medical Center/Radha/HLH/")
Beta <- read.csv(file = "11-HLHS Individual data June 21 2018 512 for heat map and PCA to Nithish.txt", header = TRUE, row.names = 1, sep = "\t")
colnames(Beta) <- gsub(".AVG_Beta", "", colnames(Beta))

Beta1 <- Beta[,1:47]
Beta1 <- t(scale(t(Beta1)))

Radha <- read.csv(file = "1L-HLH-512 markers June 21 2018 group data to Nithish.txt", header = TRUE, row.names = 1, sep = "\t")
Radha$delta <- abs(Radha$MethylationCases - Radha$MethylationControl)/100


annotation <- Radha$UCSC_REFGENE_GROUP
annotation <- as.character(annotation)
annotation[which(nchar(annotation)<2)]<-NA ### When no attotation
annotation <- strsplit(annotation,"\\;")
annotation <- sapply(annotation, "[[", 1)
logFC <- log2(Radha$Fold.change)
direction <- ifelse( logFC>0,"Hyper","Hypo")
AUC <- Radha$AUC
delta <- Radha$delta

library(ComplexHeatmap)
library(circlize)
set.seed(121)
ha = HeatmapAnnotation(df = data.frame(Samples = c(rep("HLHS", 23), rep("Control", 24))), col = list(Samples = c("HLHS" = "red4", "Control"="blue4")))
#column_tree = hclust(dist(t(meth)), method = "ward.D2")
ht_list = Heatmap(Beta1, clustering_method_columns = "ward.D2", name = "Methylation", col = colorRamp2(c(-3.0, -1, 0, 1, 3, 7), c("green4", "green", "gray","red1","firebrick","red4")),
                  cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 3, show_row_names = FALSE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(AUC, name = "AUC", col = colorRamp2(c(0.6, 0.8, 1.0), c("blue4", "purple", "red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(logFC, name = "log2(FC)", col = colorRamp2(c(2, 1, 0, -1, -2), c("red4","red","gray", "blue", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(direction, name = "Direction", col = c("Hyper" = "red4", "Hypo" = "blue4"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(delta, name = "Delta Beta", col = colorRamp2(c(0.25, 0.12, 0), c("red4","red", "blue4")), column_names_gp = gpar(fontsize = 8))+
  #Heatmap(logP, name = "-log10 (FDR)", col = colorRamp2(c(0, 10, 20), c("red4","red", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(annotation, name = "Annotation", col = c("Body" = "aquamarine", "TSS200" = "red4", "TSS1500" = "violetred", "1stExon" = "black", "5'UTR" = "blue", "3'UTR" = "green", ExonBnd= "purple"), column_names_gp = gpar(fontsize = 8))
#Heatmap(relation, name = "Probe Relationship", column_names_gp = gpar(fontsize = 8))
pdf("HLHS Heatmap.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

##########################################

tmp <- Radha[order(Radha$delta, decreasing = TRUE),]
tmp2 <- tmp[1:100,]

Beta2 <- Beta1[rownames(tmp2),]

annotation <- tmp2$UCSC_REFGENE_GROUP
annotation <- as.character(annotation)
annotation[which(nchar(annotation)<2)]<-NA ### When no attotation
annotation <- strsplit(annotation,"\\;")
annotation <- sapply(annotation, "[[", 1)
logFC <- log2(tmp2$Fold.change)
direction <- ifelse( logFC>0,"Hyper","Hypo")
AUC <- tmp2$AUC
delta <- tmp2$delta

set.seed(121)
ha = HeatmapAnnotation(df = data.frame(Samples = c(rep("HLHS", 23), rep("Control", 24))), col = list(Samples = c("HLHS" = "red4", "Control"="blue4")))
#column_tree = hclust(dist(t(meth)), method = "ward.D2")
ht_list = Heatmap(Beta2, clustering_method_columns = "ward.D2", name = "Methylation", col = colorRamp2(c(-3.0, -1, 0, 1, 3, 7), c("green4", "green", "gray","red1","firebrick","red4")),
                  cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 3, show_row_names = FALSE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(AUC, name = "AUC", col = colorRamp2(c(0.6, 0.8, 1.0), c("blue4", "purple", "red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(logFC, name = "log2(FC)", col = colorRamp2(c(2, 1, 0, -1, -2), c("red4","red","gray", "blue", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(direction, name = "Direction", col = c("Hyper" = "red4", "Hypo" = "blue4"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(delta, name = "Delta Beta", col = colorRamp2(c(0.25, 0.12, 0), c("red4","red", "blue4")), column_names_gp = gpar(fontsize = 8))+
  #Heatmap(logP, name = "-log10 (FDR)", col = colorRamp2(c(0, 10, 20), c("red4","red", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(annotation, name = "Annotation", col = c("Body" = "aquamarine", "TSS200" = "red4", "TSS1500" = "violetred", "1stExon" = "black", "5'UTR" = "blue", "3'UTR" = "green", ExonBnd= "purple"), column_names_gp = gpar(fontsize = 8))
#Heatmap(relation, name = "Probe Relationship", column_names_gp = gpar(fontsize = 8))
pdf("HLHS Heatmap Top 100.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

save.image("HLHS_Heatmap.RData")
