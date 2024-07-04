ha = HeatmapAnnotation(df = data.frame(type = c(rep("HLH", 23), rep("Control", 24))), col = list(type = c("HLH" = "red4", "Control"="blue4")))
ht_list = Heatmap(Beta.505.case, clustering_method_columns = "ward.D2", name = "Methylation", col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("blue4","lightseagreen","white","firebrick","red4")),
                  cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 4, show_row_names = FALSE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(Beta.505.control, show_row_names = FALSE)+
  Heatmap(AUC, name = "AUC", col = colorRamp2(c(0.75, 0.88, 1.0), c("blue4","gray","red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(foldchange, name = "Fold Change", col = colorRamp2(c(2.5, 1.0, 0.4), c("red","purple", "blue")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(direction, name = "Methylation Direction", col = c("Hyper" = "red", "Hypo" = "blue"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(annotation, name = "Probe Annotation", column_names_gp = gpar(fontsize = 8))+
  Heatmap(log10P, name = "-log10P", column_names_gp = gpar(fontsize = 8))+
  Heatmap(chr, name = "CHR", column_names_gp = gpar(fontsize = 8))+
  Heatmap(relation, name = "Probe Relationship", column_names_gp = gpar(fontsize = 8))
pdf("HLH-Methylation-AUC.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom")
dev.off()