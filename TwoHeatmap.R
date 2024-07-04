library(ComplexHeatmap)
library(circlize)
set.seed(12348)
ha = HeatmapAnnotation(df = data.frame(type = c(rep("Autism", 23))), col = list(type = c("Autism" = "red4")))
ha1 = HeatmapAnnotation(df = data.frame(type = c(rep("Control", 21))), col = list(type = c("Control"="blue4")))
#column_tree = hclust(dist(t(meth)), method = "ward.D2")
ht_list = Heatmap(meth1[,22:44], clustering_method_columns = "ward.D2", name = "Autism Methylation", col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("green4","aquamarine4", "violet", "firebrick3","red4")),
                  cluster_columns = FALSE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 3, show_row_names = TRUE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(foldchange, name = "Fold Change", col = colorRamp2(c(3.5, 2.0, 0.25), c("red","purple", "blue")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(direction, name = "Methylation Direction", col = c("Hyper" = "red", "Hypo" = "blue"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(annotation, name = "Probe Annotation", column_names_gp = gpar(fontsize = 8))+
  Heatmap(relation, name = "Probe Relationship", column_names_gp = gpar(fontsize = 8))+
  Heatmap(meth1[,1:21],name = " Normal Methylation", cluster_rows = FALSE, col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("green4","aquamarine4", "violet", "firebrick3","red4")), 
          cluster_columns = FALSE, show_column_names = TRUE, top_annotation = ha1, column_names_gp = gpar(fontsize = 8), show_row_names = TRUE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8))
pdf("Autism-Methylation1.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom")
dev.off()