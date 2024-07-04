library(ComplexHeatmap)
library(circlize)

annotation <- SIDS.Meth$UCSC_REFGENE_GROUP
annotation <- as.character(annotation)
annotation[which(nchar(annotation)<2)]<-NA ### When no attotation
annotation <- strsplit(annotation,"\\;")
annotation <- sapply(annotation, "[[", 1)
direction <- ifelse(SIDS.FC$Fold.chance >1,"Hyper","Hypo")
foldchange <- SIDS.FC$Fold.chance
relation <- SIDS.Meth$RELATION_TO_UCSC_CPG_ISLAND
AUC <- SIDS.FC$AUC

Meth <- SIDS.Meth[,grep("AVG_Beta", colnames(SIDS.Meth))]
colnames(Meth) <- gsub(".AVG_Beta", "", colnames(Meth))
colnames(Meth) <- sapply(strsplit(colnames(Meth),"\\."),"[[",2)

set.seed(1100)
ha = HeatmapAnnotation(df = data.frame(type = c(rep("SIDS", length(grep("Affect", colnames(Meth)))), rep("Control", length(grep("Nor", colnames(Meth)))))), col = list(type = c("Control" = "blue4", "SIDS"="red")))
#column_tree = hclust(dist(t(meth)), method = "ward.D2")
ht_list = Heatmap(Meth, name = "Methylation", col = colorRamp2(c(0, 0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 1), c("green4","lightseagreen","aquamarine4", "violet","violetred1","violetred3","firebrick","red4")),
                  cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 3, show_row_names = TRUE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(AUC, name = "AUC", col = colorRamp2(c(0.75, 0.8, 0.9, 1.0), c("blue4","purple","red", "red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(foldchange, name = "Fold Change", col = colorRamp2(c(2.5, 1.0, 0.4), c("red","purple", "blue")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(direction, name = "Methylation Direction", col = c("Hyper" = "red", "Hypo" = "blue"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(annotation, name = "Probe Annotation", column_names_gp = gpar(fontsize = 8), col = c("1stExon" = "blue4", "TSS200"="red", "TSS1500"="red4", "Body"="violetred1", "3'UTR"="aquamarine4", "5'UTR"="green4"))+
  Heatmap(relation, name = "Probe Relationship", column_names_gp = gpar(fontsize =8))
pdf("SIDS-Methylation.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom")
dev.off()
#col = c("1stExon" = "blue4", "TSS200"="red", "TSS1500"="red4", "Body"="gray", "3'UTR"="green", "5'UTR"="green4")