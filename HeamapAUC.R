library(ComplexHeatmap)
library(circlize)
annotation <- TOF.Beta.fold$UCSC_REFGENE_GROUP
annotation <- as.character(annotation)
annotation[which(nchar(annotation)<2)]<-NA ### When no attotation
annotation <- strsplit(annotation,"\\;")
annotation <- sapply(annotation, "[[", 1)
direction <- ifelse(TOF.Beta.fold.AUC$Fold.chance >1,"Hyper","Hypo")
foldchange <- TOF.Beta.fold.AUC$Fold.chance
relation <- TOF.Beta.fold.AUC$RELATION_TO_UCSC_CPG_ISLAND
ha = HeatmapAnnotation(df = data.frame(type = c(rep("TOF", 23), rep("Control", 24))), col = list(type = c("TOF" = "red4", "Control"="blue4")))
#column_tree = hclust(dist(t(meth)), method = "ward.D2")
ht_list = Heatmap(meth.AUC, clustering_method_columns = "ward.D2", name = "Methylation", col = colorRamp2(c(0, 0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 1), c("green4","lightseagreen","aquamarine4", "violet","violetred3","red", "firebrick4","red4")),
                  cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 4, show_row_names = TRUE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(AUC, name = "AUC", col = colorRamp2(c(0.75, 0.8, 0.9, 1.0), c("blue4","purple","red", "red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(foldchange, name = "Fold Change", col = colorRamp2(c(2.5, 1.0, 0.4), c("red","purple", "blue")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(direction, name = "Methylation Direction", col = c("Hyper" = "red", "Hypo" = "blue"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(annotation, name = "Probe Annotation", column_names_gp = gpar(fontsize = 8))+
  Heatmap(relation, name = "Probe Relationship", column_names_gp = gpar(fontsize = 8))
pdf("TOF75-Methylation-AUC.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom")
dev.off()