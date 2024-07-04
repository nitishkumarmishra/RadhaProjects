library(ComplexHeatmap)
library(circlize)
methylation <- read.csv("0-HS Individual data 304 targets July 27 2016 to Nitish.txt", header = TRUE, sep = "\t")
results <- read.csv("00-Final to Nitish wih 304 targets July 27 2016.txt", header = TRUE, sep = "\t")
MergeData <- merge(methylation, results, by="TargetID")
rownames(MergeData) <- MergeData$TargetID
Meth <- MergeData[,grep("AVG_Beta", colnames(MergeData))]
colnames(Meth) <- gsub(".AVG", "_AVG", colnames(Meth))
tmp <- strsplit(colnames(Meth),"\\.")
colnames(Meth) <- sapply(tmp, "[[", 2)
##########################

case <- Meth[,grep("HS", colnames(Meth))]
case$CaseMean <- rowMeans(case)
Control <- Meth[,grep("CTR", colnames(Meth))]
Control$ControlMean <- rowMeans(Control)
Meth <- cbind(case, Control)
Meth$deltaBeta <- Meth$CaseMean - Meth$ControlMean
########################################################################
annotation <- MergeData$UCSC_REFGENE_GROUP
annotation <- as.character(annotation)
annotation[which(nchar(annotation)<2)]<-NA ### When no attotation
annotation <- strsplit(annotation,"\\;")
annotation <- sapply(annotation, "[[", 1)
direction <- ifelse(MergeData$Fold..chance. >1,"Hyper","Hypo")
foldchange <- MergeData$Fold..chance.
relation <- MergeData$
AUC <- MergeData$AUC
logP <- MergeData$LOG10p
delta <- Meth$deltaBeta
#########################################################################
Meth1 <- Meth; Meth1$CaseMean <- NULL; Meth1$ControlMean <- NULL; Meth1$deltaBeta <- NULL
Meth1 <- as.matrix(Meth1)
#########################################################################
set.seed(1234)
ha = HeatmapAnnotation(df = data.frame(type = c(rep("HS", 24), rep("Control", 24))), col = list(type = c("HS" = "red4", "Control"="blue4")))
#column_tree = hclust(dist(t(meth)), method = "ward.D2")
ht_list = Heatmap(Meth1, clustering_method_columns = "ward.D2", name = "Methylation", col = colorRamp2(c(0, 0.15, 0.30, 0.45, 0.6, 0.75, 1), c("green4","lightgreen", "aquamarine4", "violet","violetred3","firebrick","red4")),
                  cluster_columns = FALSE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 3, show_row_names = FALSE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(AUC, name = "AUC", col = colorRamp2(c(0.75, 0.88, 1.0), c("blue4", "purple", "red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(foldchange, name = "Fold Change", col = colorRamp2(c(3.5, 2.5, 2), c("red","purple", "blue")), column_names_gp = gpar(fontsize = 8))+
  #Heatmap(direction, name = "Methylation Direction", col = c("Hyper" = "red", "Hypo" = "blue"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(delta, name = "Delta beta", col = colorRamp2(c(0, 0.1, 0.2), c("red4","gray", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(logP, name = "log P-value", col = colorRamp2(c(0, 20, 40), c("red4","gray", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(annotation, name = "Probe Annotation", column_names_gp = gpar(fontsize = 8))
  #Heatmap(relation, name = "Probe Relationship", column_names_gp = gpar(fontsize = 8))
pdf("HS 304 AUC.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom")
dev.off()
save.image("HSdata.RData")
