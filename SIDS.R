library(ComplexHeatmap)
library(circlize)

methylation <- read.csv("1-SIDS_SG_Individual_Sept 11 2016 Heatmap sent to Nitish.txt", header = TRUE, sep = "\t", row.names = 1)
results <- read.csv("2 SIDS_SG Group36targets Sept 11 2016 to Nitish heat map.txt", header = TRUE, sep = "\t")
colnames(results)[1] <- "TargetID"
MergeData <- merge(methylation, results, by="TargetID")
rownames(MergeData) <- MergeData$TargetID
Meth <- MergeData[,grep("AVG_Beta", colnames(MergeData))]
colnames(Meth) <- sapply(strsplit(colnames(Meth),"\\."), "[[", 2)


case <- Meth[,grep("Affect", colnames(Meth))]
case$CaseMean <- rowMeans(case)
Control <- Meth[,grep("Nor", colnames(Meth))]
Control$ControlMean <- rowMeans(Control)
Meth <- cbind(case, Control)
Meth$deltaBeta <- Meth$CaseMean - Meth$ControlMean

annotation <- MergeData$UCSC_REFGENE_GROUP
annotation <- as.character(annotation)
annotation[which(nchar(annotation)<2)]<-NA ### When no attotation
annotation <- strsplit(annotation,"\\;")
annotation <- sapply(annotation, "[[", 1)
direction <- ifelse(MergeData$Fold.chance >1,"Hyper","Hypo")
foldchange <- MergeData$Fold.chance
relation <- MergeData$RELATION_TO_UCSC_CPG_ISLAND
AUC <- MergeData$AUC
logP <- MergeData$LOG10p
delta <- Meth$deltaBeta

Meth1 <- Meth; Meth1$CaseMean <- NULL; Meth1$ControlMean <- NULL; Meth1$deltaBeta <- NULL
Meth1 <- as.matrix(Meth1)
#########################################################################
set.seed(1110)
ha = HeatmapAnnotation(df = data.frame(type = c(rep("SIDS", 17), rep("Control", 7))), col = list(type = c("SIDS" = "red4", "Control"="blue4")))
#column_tree = hclust(dist(t(meth)), method = "ward.D2")
ht_list = Heatmap(Meth1, clustering_method_columns = "ward.D2", name = "Methylation", col = colorRamp2(c(0, 0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 1), c("green4","lightseagreen","aquamarine4", "violet","violetred1","violetred3","firebrick","red4")),
                  cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 3, show_row_names = TRUE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(AUC, name = "AUC", col = colorRamp2(c(0.75, 0.88, 1.0), c("blue4", "purple", "red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(foldchange, name = "Fold Change", col = colorRamp2(c(3.5, 2.5, 2), c("red","purple", "blue")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(direction, name = "Methylation Direction", col = c("Hyper" = "red", "Hypo" = "blue"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(delta, name = "Delta beta", col = colorRamp2(c(0, 0.1, 0.2), c("red4","gray", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(logP, name = "log P-value", col = colorRamp2(c(-20, 20, 40), c("red4","gray", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(annotation, name = "Probe Annotation", column_names_gp = gpar(fontsize = 8))+
  Heatmap(relation, name = "Probe Relationship", column_names_gp = gpar(fontsize = 8))
pdf("SIDS_36.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom")
dev.off()

################################
Meth2 <- Meth1[which(abs(Meth$deltaBeta) > 0.2),]
AUC1 <- AUC[which(abs(Meth$deltaBeta) > 0.2)]
foldchange1 <- foldchange[which(abs(Meth$deltaBeta) > 0.2)]
direction1 <- direction[which(abs(Meth$deltaBeta) > 0.2)]
delta1 <- delta[which(abs(Meth$deltaBeta) > 0.2)]
logP1 <- logP[which(abs(Meth$deltaBeta) > 0.2)]
annotation1 <- annotation[which(abs(Meth$deltaBeta) > 0.2)]
relation1 <- relation[which(abs(Meth$deltaBeta) > 0.2)]
set.seed(1110)
ha = HeatmapAnnotation(df = data.frame(type = c(rep("SIDS", 17), rep("Control", 7))), col = list(type = c("SIDS" = "red4", "Control"="blue4")))
#column_tree = hclust(dist(t(meth)), method = "ward.D2")
ht_list = Heatmap(Meth2, clustering_method_columns = "ward.D2", name = "Methylation", col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("green4", "aquamarine4", "brown2", "brown4", "red4")),
                  cluster_columns = FALSE, cluster_rows = TRUE, km=2, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), show_row_names = TRUE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(AUC1, name = "AUC", col = colorRamp2(c(0.6, 0.8, 1), c("green4","red", "red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(foldchange1, name = "Fold Change", col = colorRamp2(c(3.5, 2.5, 2), c("red","purple", "blue")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(direction1, name = "Meth. Direction", col = c("Hyper" = "red", "Hypo" = "blue"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(delta1, name = "Delta beta", col = colorRamp2(c(-0.25, 0.1, 0.35), c("red4","gray", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(logP1, name = "log P-value", col = colorRamp2(c(-20, 20, 50), c("red4","gray", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(annotation1, name = "Probe Annot.", column_names_gp = gpar(fontsize = 8))
  #Heatmap(relation1, name = "Probe Relationship", column_names_gp = gpar(fontsize = 8))
pdf("SIDS_beta2.pdf", width = 8, height = 8)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom")
dev.off()

save.image("SIDS.RData")
