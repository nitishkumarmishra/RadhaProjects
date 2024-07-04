library(ComplexHeatmap)
library(circlize)
##############
TOF.75 <- read.csv("10-TOF -Individula-75-markers.txt", header = TRUE, row.names = 1, sep = "\t")
TOF.75$X <- NULL
TOF.75$X.1 <- NULL
TOF.75$X.2 <- NULL
col <- colnames(TOF.75)
col <- gsub("\\.", "-", col)
colnames(TOF.75) <- col
rownames(TOF.75) <- TOF.75$TargetID
TOF.Beta <- TOF.75[,grep("AVG_Beta", colnames(TOF.75))]
Probes.75 <- rownames(TOF.Beta)
fold.10TOF <- read.csv("10A-TOF-fold-change-without-XYRSgaps-1.txt", header = TRUE, sep = "\t")
rownames(fold.10TOF) <- fold.10TOF$TargetID
TOF.Beta.fold <- fold.10TOF[Probes.75,]
annotation <- TOF.Beta.fold$UCSC_REFGENE_GROUP
annotation <- as.character(annotation)
annotation <- strsplit(annotation,"\\;")
annotation <- sapply(annotation, "[[", 1)
direction <- ifelse(TOF.Beta.fold$Fold.chance >1,"Hyper","Hypo")
foldchange <- TOF.Beta.fold$Fold.chance
relation <- TOF.Beta.fold$RELATION_TO_UCSC_CPG_ISLAND
#####################
meth <- TOF.Beta
col <- colnames(meth)
col <- gsub("-TOF-AVG_Beta","",col)
col <- gsub("-AVG_Beta", "", col)
a <- strsplit(colnames(meth), split = "-")
tmp <- sapply(a, `[`,2)
colnames(meth) <- tmp
####################
ha = HeatmapAnnotation(df = data.frame(type = c(rep("TOF", 23), rep("Control", 24))), col = list(type = c("TOF" = "red4", "Control"="blue4")))
#column_tree = hclust(dist(t(meth)), method = "ward.D2")
ht_list = Heatmap(meth, clustering_method_columns = "ward.D2", name = "Methylation", col = colorRamp2(c(0, 0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 1), c("green4","lightseagreen","aquamarine4", "violet","violetred3","red", "firebrick4","red4")),
                  cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 4, show_row_names = TRUE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(foldchange, name = "Fold Change", col = colorRamp2(c(2.5, 1.0, 0.4), c("red","purple", "blue")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(direction, name = "Methylation Direction", col = c("Hyper" = "red", "Hypo" = "blue"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(annotation, name = "Probe Annotation", column_names_gp = gpar(fontsize = 8))+
  Heatmap(relation, name = "Probe Relationship", column_names_gp = gpar(fontsize = 8))
pdf("TOF75-Methylation.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom")
dev.off()