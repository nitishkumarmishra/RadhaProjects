library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

setwd("F:/OneDrive - University of Nebraska Medical Center/Radha/Myopia/")

TableInput <- read.csv("5-HM new Individual data-May 16 1541 Calculated for Heatmap to Nithish.txt", header = TRUE, sep = "\t", row.names = 2)
select <- grep("AVG_Beta", names(TableInput), value = TRUE)
Meth <- TableInput[, select]
colnames(Meth) <- gsub(".AVG_Beta", "", colnames(Meth))
colnames(Meth) <- gsub("UR.", "", colnames(Meth))

Delta <- read.csv("1-High myopia_Table 1_May 16 2018.txt", header = TRUE, sep = "\t", row.names = 1)
###### Select top 20 CpG based on delta beta values ###########
draw_rownames <- function(.data) .data %>% do(mutate(.,rownames=rownames(.)))
tmp <- Delta %>% draw_rownames() %>% arrange(desc(PercentageDifference)) %>% slice(1:20)

Beta <- Meth[tmp$rownames,]
Beta1 = t(scale(t(Beta)))


annotation <- TableInput[rownames(Beta), c("UCSC_REFGENE_GROUP")]
annotation <- as.character(annotation)
annotation[which(nchar(annotation)<2)]<-NA ### When no attotation
annotation <- strsplit(annotation,"\\;")
annotation <- sapply(annotation, "[[", 1)

direction <- ifelse(tmp$Fold.change >0,"Hyper","Hypo")
foldchange <- tmp$Fold.change
AUC <- tmp$AUC

logP <- -log10(tmp$FDR.p.Val)
delta <- tmp$PercentageDifference/100

########################
set.seed(121)
ha = HeatmapAnnotation(df = data.frame(Samples = c(rep("HM", 18), rep("Control", 18))), col = list(Samples = c("HM" = "red4", "Control"="blue4")))
#column_tree = hclust(dist(t(meth)), method = "ward.D2")
ht_list = Heatmap(Beta1, clustering_method_columns = "ward.D2", name = "Methylation", col = colorRamp2(c(-2.0, -1, 0, 2, 5), c("green4","green1","red1", "firebrick","red4")),
                  cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 2, show_row_names = TRUE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(AUC, name = "AUC", col = colorRamp2(c(0.75, 0.88, 1.0), c("blue4", "purple", "red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(foldchange, name = "log2(FC)", col = colorRamp2(c(5, 2, 1, 0), c("red4","red","purple", "blue")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(direction, name = "Direction", col = c("Hyper" = "red", "Hypo" = "blue"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(delta, name = "Delta Beta", col = colorRamp2(c(0.3, 0.15, 0), c("red4","red", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(logP, name = "-log10 (FDR)", col = colorRamp2(c(0, 25, 50), c("red4","red", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(annotation, name = "Annotation", col = c("Body" = "aquamarine", "TSS200" = "red4", "TSS1500" = "violetred", "1stExon" = "black", "5'UTR" = "blue", "3'UTR" = "green"), column_names_gp = gpar(fontsize = 8))
#Heatmap(relation, name = "Probe Relationship", column_names_gp = gpar(fontsize = 8))
pdf("Myopia ComplexHeatmap Top 20.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

save.image("Myopia_Heatmap.RData")
