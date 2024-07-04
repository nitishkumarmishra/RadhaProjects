library(ComplexHeatmap)
library(circlize)
ha = HeatmapAnnotation(df = data.frame(type = c(rep("disease", length(AS.Samples)), rep("normal", length(Non.AS.Samples)))), col = list(type = c("normal" = "blue", "disease"="red")))
column_tree = hclust(dist(t(meth)), method = "ward.D2")
ht_list = Heatmap(meth[, column_tree$order], name = "Methylation", col = colorRamp2(c(0, 0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 1), c("green4","lightseagreen","aquamarine4", "aquamarine1","violet","violetred3","firebrick","red4")),
                  cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 3, show_row_names = TRUE, column_title = "Methylation", column_title_gp = gpar(fontsize = 8)) +
  Heatmap(finalFoldChange, name = "Fold Change", col = colorRamp2(c(2.5, 1.0, 0.4), c("red","purple", "blue")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(direction, name = "Methylation Direction", col = c("Hyper" = "red", "Hypo" = "blue"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(finalRelation, name = "Probe Relationship", column_names_gp = gpar(fontsize = 8))+
  Heatmap(finalUCSCgroup, name = "Probe Annotation", column_names_gp = gpar(fontsize = 8))
pdf("AS-Methylation.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE)
dev.off()