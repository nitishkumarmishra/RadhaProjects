ha = HeatmapAnnotation(df = data.frame(type = c(rep("HS", length(HS.ID)), rep("Control", length(CTR.ID)))), col = list(type = c("HS" = "red4", "Control"="blue4")))
#column_tree = hclust(dist(t(meth)), method = "ward.D2")
ht_list = Heatmap(meth.AUC.0.75, clustering_method_columns = "ward.D2", name = "Methylation", col = colorRamp2(c(0, 0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 1), c("green4","aquamarine1", "violetred1", "violetred2","violetred3","firebrick", "red3", "red4")),
                  cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 3, show_row_names = FALSE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(foldchange, name = "Fold Change", col = colorRamp2(c(3.1, 2.5, 2.0), c("red","purple", "blue")), column_names_gp = gpar(fontsize = 8))+
  #Heatmap(direction, name = "Meth. Direction", col = c("Hyper" = "red", "Hypo" = "blue"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(percentChange, name = "% Meth. Change", column_names_gp = gpar(fontsize = 8), col = colorRamp2(c(10, 12.5, 15, 17.5, 20, 22.5, 25), c("green4","aquamarine4", "violet","violetred1","violetred3","firebrick","red4")))+
  Heatmap(log10FDR, name = "-log10(FDR)", column_names_gp = gpar(fontsize = 8))+
  Heatmap(AUC, name = "AUC", column_names_gp = gpar(fontsize = 8), col = colorRamp2(c(1.0, 0.9, 0.8, 0.7), c("red4","firebrick","purple", "blue")))+
  Heatmap(annotation, name = "Probe Annot.", column_names_gp = gpar(fontsize = 8))+
  Heatmap(relation, name = "Probe Relation", column_names_gp = gpar(fontsize = 8))
pdf("HS-AUC-0.75.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom")
dev.off()