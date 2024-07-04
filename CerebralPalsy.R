library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

setwd("F:/OneDrive - University of Nebraska Medical Center/Radha/Cerebral Palcy/Latest/")
TableInput <- fread(input = "CP Analysis 341 Targets for AUC ROC Nov 23 2018 to Nitish regression.txt")
select <- grep("AVG_Beta|TargetID", names(TableInput), value = TRUE)
Meth <- TableInput[, select, with=FALSE]
colnames(Meth) <- gsub(".AVG_Beta", "", colnames(Meth))
Meth <- as.data.frame(Meth); Meth <- na.omit(Meth); rownames(Meth) <- Meth$TargetID; Meth$TargetID <- NULL
colnames(Meth) <- gsub("Autism_Control", "CP_Control", colnames(Meth))

tmp <- as.data.frame(TableInput[,c( "TargetID", "UCSC_REFGENE_GROUP", "RELATION_TO_UCSC_CPG_ISLAND")])
rownames(tmp) <- tmp$TargetID; tmp$TargetID <- NULL


Delta <- read.csv("6-CP_group_analysis 240 targest-Sep 27 2018.txt", header = TRUE, sep = "\t", row.names = 2)
Delta <- merge(Delta, tmp, by="row.names")
rownames(Delta) <- Delta$Row.names
Delta$PercentageDifference <- abs(Delta$MethylationCases - Delta$MethylationControl)/ 100

#Delta <- Delta[TableInput$TargetID,]
Delta$logFC <- log2(Delta$Fold.chance)
#Delta <- Delta[na.omit(Delta$PercentageDifference),]

Beta = t(scale(t(Meth)))
Beta <- Beta[rownames(Delta),]
#Beta = t(apply(Meth, 1, scale)) ## Both line scale in same way


annotation <- Delta$UCSC_REFGENE_GROUP
annotation <- as.character(annotation)
annotation[which(nchar(annotation)<2)]<-NA ### When no attotation
annotation <- strsplit(annotation,"\\;")
annotation <- sapply(annotation, "[[", 1)


direction <- ifelse(Delta$logFC >0,"Hyper","Hypo")
foldchange <- Delta$logFC
AUC <- Delta$AUC

logP <- -log10(Delta$FDR.p.Val)
delta <- Delta$PercentageDifference

########################
set.seed(121)
ha = HeatmapAnnotation(df = data.frame(Samples = c(rep("Cerebral", 23), rep("Control", 21))), col = list(Samples = c("Cerebral" = "red4", "Control"="blue4")))
#column_tree = hclust(dist(t(meth)), method = "ward.D2")
ht_list = Heatmap(Beta, clustering_method_columns = "ward.D2", name = "Methylation", col = colorRamp2(c(-2, -1, 0, 2, 4, 6), c("green4","green","white","red","firebrick","red4")),
                  cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 5, show_row_names = FALSE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(AUC, name = "AUC", col = colorRamp2(c(0.5, 0.75, 1.0), c("blue", "purple", "red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(foldchange, name = "log2(FC)", col = colorRamp2(c(2, 1, 0, -1, -3), c("red4","red","white", "blue", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(direction, name = "Direction", col = c("Hyper" = "red4", "Hypo" = "blue4"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(delta, name = "Delta Beta", col = colorRamp2(c(0.3, 0.15, 0), c("red4","red", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(logP, name = "-log10 (FDR)", col = colorRamp2(c(0, 15, 30), c("blue4","red","red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(annotation, name = "Annotation", col = c("Body" = "aquamarine", "TSS200" = "red4", "TSS1500" = "violetred", "1stExon" = "black", "5'UTR" = "blue", "3'UTR" = "green"), column_names_gp = gpar(fontsize = 8))
#Heatmap(relation, name = "Probe Relationship", column_names_gp = gpar(fontsize = 8))
pdf("Cerebral Scaled 240 CpGs.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()


###############################
###### Select top 20 CpG based on delta beta values ###########
draw_rownames <- function(.data) .data %>% do(mutate(.,rownames=rownames(.)))
tmp <- Delta %>% draw_rownames() %>% arrange(desc(PercentageDifference)) %>% slice(1:50)

Beta <- Meth[tmp$rownames,]
Beta = t(scale(t(Beta)))

#TableInput1 <- TableInput[TableInput$TargetID%in%rownames(Beta),]
Delta1 <- Delta[tmp$Row.names,]
Meth1 <- Beta[tmp$Row.names,]

annotation <- Delta1$UCSC_REFGENE_GROUP
annotation <- as.character(annotation)
annotation[which(nchar(annotation)<2)]<-NA ### When no attotation
annotation <- strsplit(annotation,"\\;")
annotation <- sapply(annotation, "[[", 1)
direction <- ifelse(Delta1$logFC >0,"Hyper","Hypo")
foldchange <- Delta1$logFC
AUC <- Delta1$AUC
logP <- -log10(Delta1$FDR.p.Val)
delta <- Delta1$PercentageDifference

##########################################
set.seed(121)
ha = HeatmapAnnotation(df = data.frame(Samples = c(rep("Cerebral", 23), rep("Control", 21))), col = list(Samples = c("Cerebral" = "red4", "Control"="blue4")))
#column_tree = hclust(dist(t(meth)), method = "ward.D2")
ht_list = Heatmap(Meth1, clustering_method_columns = "ward.D2", name = "Methylation", col = colorRamp2(c(-2, -1, 0, 2, 4, 6), c("green4","green","white","red","firebrick","red4")),
                  cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 3, show_row_names = FALSE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(AUC, name = "AUC", col = colorRamp2(c(0.5, 0.75, 1.0), c("blue", "purple", "red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(foldchange, name = "log2(FC)", col = colorRamp2(c(2, 1, 0, -1, -3), c("red4","red","white", "blue", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(direction, name = "Direction", col = c("Hyper" = "red4", "Hypo" = "blue4"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(delta, name = "Delta Beta", col = colorRamp2(c(0.3, 0.15, 0), c("red4","red", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(logP, name = "-log10 (FDR)", col = colorRamp2(c(0, 15, 30), c("blue4","red","red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(annotation, name = "Annotation", col = c("Body" = "aquamarine", "TSS200" = "red4", "TSS1500" = "violetred", "1stExon" = "black", "5'UTR" = "blue", "3'UTR" = "green"), column_names_gp = gpar(fontsize = 8))
#Heatmap(relation, name = "Probe Relationship", column_names_gp = gpar(fontsize = 8))
pdf("Cerebral Scaled 50 CpGs.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
##########################################
save.image("Cerebral_Heatmap.RData")
