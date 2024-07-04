ha = HeatmapAnnotation(df = data.frame(type = c(rep("HS", length(HS.ID)), rep("Control", length(CTR.ID)))), col = list(type = c("HS" = "red4", "Control"="green4")))
#column_tree = hclust(dist(t(meth1)), method = "ward.D2")
ht_list = Heatmap(meth.tss, name = "Methylation", clustering_method_columns = "ward.D2", col = colorRamp2(c(0, 0.10, 0.2, 0.25, 0.50, 0.75, 1), c("darkgreen", "green", "tomato", "violetred", "red2", "red3","red4")),
                  cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 4, show_row_names =FALSE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(foldchange.TSS, name = "Fold Change", col = colorRamp2(c(3.5, 1.0, 0.3), c("red","purple", "blue")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(direction.TSS, name = "Methylation Direction", col = c("Hyper" = "red", "Hypo" = "blue"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(annotation.TSS, name = "Probe Annotation", column_names_gp = gpar(fontsize = 8))+
  Heatmap(relation.TSS, name = "Probe Relationship", column_names_gp = gpar(fontsize = 8))
pdf("HS-Methylation-TSS.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom")
dev.off()
