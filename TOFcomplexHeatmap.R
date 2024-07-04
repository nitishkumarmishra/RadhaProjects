library(ComplexHeatmap)
library(circlize)
a <- strsplit(colnames(meth), split = "-")
tmp <- sapply(a, `[`,2)
colnames(meth) <- tmp
ha = HeatmapAnnotation(df = data.frame(type = c(rep("TOF", length(Non.AS.Samples)), rep("Control", length(AS.Samples)))), col = list(type = c("TOF" = "red4", "Control"="blue")))
#column_tree = hclust(dist(t(meth)), method = "ward.D2")
ht_list = Heatmap(meth, name = "Methylation", clustering_method_columns = "ward.D2", col = colorRamp2(c(0, 0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 1), c("green4","lightseagreen","aquamarine4","violet","violetred1", "violetred3","firebrick","red4")),
                  cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 2, show_row_names = TRUE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(foldChange, name = "Fold Change", col = colorRamp2(c(2.5, 1.0, 0.4), c("red","purple", "blue")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(direction, name = "Methylation Direction", col = c("Hyper" = "red", "Hypo" = "blue"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(relation, name = "Probe Relationship", column_names_gp = gpar(fontsize = 8))+
  Heatmap(annot, name = "Probe Annotation", column_names_gp = gpar(fontsize = 8))
pdf("TOF-Methylation-New.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom")
dev.off()