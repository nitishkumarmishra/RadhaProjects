library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

setwd("F:/OneDrive - University of Nebraska Medical Center/Radha/Autism/NewAutism/")
TableInput <- fread(input = "7-Autism Blood Individual 230 CpG AUC Aug 10 2018 to Nitish.txt")
select <- grep("AVG_Beta|TargetID", names(TableInput), value = TRUE)
Meth <- TableInput[, select, with=FALSE]
colnames(Meth) <- gsub(".AVG_Beta", "", colnames(Meth))
#colnames(Meth) <- gsub("PanCancer", "", colnames(Meth))
Meth <- as.data.frame(Meth); Meth <- na.omit(Meth); rownames(Meth) <- Meth$TargetID; Meth$TargetID <- NULL

#TableInput <- read.csv("1-Pancreatic cancer_Individual data_263 targets_June 29 2018-To Nithish.txt", header = TRUE, sep = "\t", row.names = 1)
#select <- grep("AVG_Beta", names(TableInput), value = TRUE)
#Meth <- TableInput[, select]
#colnames(Meth) <- gsub(".AVG_Beta", "", colnames(Meth))
#colnames(Meth) <- gsub("PanCancer", "", colnames(Meth))

Delta <- read.csv("6-Autism Blood Group 230 CpG AUC Aug 10 2018 Nitish with AUC.txt", header = TRUE, sep = "\t", row.names = 1)
Delta$PercentageDifference <- abs(Delta$MethylationCases-Delta$MethylationControl)
Delta <- Delta[TableInput$TargetID,]
Beta = t(scale(t(Meth)))

annotation <- Delta$UCSC_REFGENE_GROUP
annotation <- as.character(annotation)
annotation[which(nchar(annotation)<2)]<-NA ### When no attotation
annotation <- strsplit(annotation,"\\;")
annotation <- sapply(annotation, "[[", 1)
Delta$log2FC <- log2(Delta$Fold.change)
direction <- ifelse(Delta$log2FC >0,"Hyper","Hypo")
foldchange <- Delta$log2FC
AUC <- Delta$AUC
relation <- Delta$RELATION_TO_UCSC_CPG_ISLAND
logP <- -log10(Delta$FDR.p.Val)
delta <- Delta$PercentageDifference/100

########################
set.seed(123)
ha = HeatmapAnnotation(df = data.frame(Samples = c(rep("Control", 10), rep("Autism", 14))), col = list(Samples = c("Control" = "blue4", "Autism"="red4")))
#column_tree = hclust(dist(t(meth)), method = "ward.D2")
ht_list = Heatmap(Beta, clustering_method_columns = "ward.D2", name = "Methylation", col = colorRamp2(c(-3, -1.5,0, 1, 2, 3, 5), c("green4","green2","white","red1","firebrick","red3","red4")),
                  cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 5, show_row_names = FALSE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(AUC, name = "AUC", col = colorRamp2(c(0.5, 0.75, 1.0), c("blue", "purple", "red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(foldchange, name = "log2(FC)", col = colorRamp2(c(3, 1, -1, -3), c("red4","red","purple", "blue")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(direction, name = "Direction", col = c("Hyper" = "red4", "Hypo" = "blue4"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(delta, name = "Delta Beta", col = colorRamp2(c(0.5, 0.4, 0.3,0.2), c("red4","red", "blue1","blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(logP, name = "-log10 (FDR)", col = colorRamp2(c(0, 20, 40, 60), c("blue4","red1","red3","red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(annotation, name = "Annotation", col = c("Body" = "aquamarine", "TSS200" = "red4", "TSS1500" = "violetred", "1stExon" = "black", "5'UTR" = "blue", "3'UTR" = "green"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(relation, name = "Probe Relationship", column_names_gp = gpar(fontsize = 8))
pdf("Autism New 230 August 2019.pdf", width = 10, height = 10)
#png("Autism New 230 August 2019.png", width = 10, height = 10, units="in", res=500)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()


###############################

###############################
###### Select top 20 CpG based on delta beta values ###########
draw_rownames <- function(.data) .data %>% do(mutate(.,rownames=rownames(.)))
tmp <- Delta %>% draw_rownames() %>% arrange(desc(PercentageDifference)) %>% slice(1:25)

Beta1 <- Meth[rownames(Meth)%in%tmp$rownames,]
Beta1 = t(scale(t(Beta1)))

TableInput1 <- TableInput[TableInput$TargetID%in%rownames(Beta1),]
Delta1 <- Delta[TableInput1$TargetID,]


annotation <- Delta1$UCSC_REFGENE_GROUP
annotation <- as.character(annotation)
annotation[which(nchar(annotation)<2)]<-NA ### When no attotation
annotation <- strsplit(annotation,"\\;")
annotation <- sapply(annotation, "[[", 1)
Delta1$log2FC <- log2(Delta1$Fold.change)
direction <- ifelse(Delta1$log2FC >0,"Hyper","Hypo")
foldchange <- Delta1$log2FC
AUC <- Delta1$AUC
relation <- Delta1$RELATION_TO_UCSC_CPG_ISLAND
logP <- -log10(Delta1$FDR.p.Val)
delta <- Delta1$PercentageDifference/100
##########################################
set.seed(121)
ha = HeatmapAnnotation(df = data.frame(Samples = c(rep("Control", 10), rep("Autism", 14))), col = list(Samples = c("Control" = "blue4", "Autism"="red4")))
#column_tree = hclust(dist(t(meth)), method = "ward.D2")
ht_list = Heatmap(Beta1, clustering_method_columns = "ward.D2", name = "Methylation", col = colorRamp2(c(-4, -2, 0, 1, 2, 3, 5), c("green4","green2","white","red","firebrick","red3","red4")),
                  cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 2, show_row_names = FALSE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(AUC, name = "AUC", col = colorRamp2(c(0.5, 0.75, 1.0), c("blue", "purple", "red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(foldchange, name = "log2(FC)", col = colorRamp2(c(3, 1, -1, -3), c("red4","red","purple", "blue")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(direction, name = "Direction", col = c("Hyper" = "red4", "Hypo" = "blue4"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(delta, name = "Delta Beta", col = colorRamp2(c(0.5, 0.4, 0.3,0.2), c("red4","red", "blue1","blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(logP, name = "-log10 (FDR)", col = colorRamp2(c(0, 20, 40, 60), c("blue4","red1","red3","red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(annotation, name = "Annotation", col = c("Body" = "aquamarine", "TSS200" = "red4", "TSS1500" = "violetred", "1stExon" = "black", "5'UTR" = "blue", "3'UTR" = "green"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(relation, name = "Probe Relationship", column_names_gp = gpar(fontsize = 8))
pdf("Autism New top 25.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()



