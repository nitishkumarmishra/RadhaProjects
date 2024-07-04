library(ComplexHeatmap)
library(circlize)
setwd("F:/OneDrive - University of Nebraska Medical Center/Radha/HS-marker-data-Heat Map/HS Network/Latest/New Data/")
methylation <- read.csv("4-HS targets Individual 277 targets Feb 21 2018 to Nitish for Heatmap.txt", header = TRUE, sep = "\t")
results <- read.csv("4-HS targets Table-1 for MS 277 targets Group data Feb 21 2018 to Nitish.txt", header = TRUE, sep = "\t")
MergeData <- merge(methylation, results, by="TargetID")
rownames(MergeData) <- MergeData$TargetID
Meth <- MergeData[,grep("AVG_Beta", colnames(MergeData))]
colnames(Meth) <- gsub(".AVG", "_AVG", colnames(Meth))
tmp <- strsplit(colnames(Meth),"\\.")
colnames(Meth) <- sapply(tmp, "[[", 2)
colnames(Meth) <- gsub("_AVG_Beta", "", colnames(Meth))
##########################

case <- Meth[,grep("HS", colnames(Meth))]
Control <- Meth[,grep("CTR", colnames(Meth))]
Meth <- cbind(Control, case)
M=log2(Meth/(1-  Meth))
status <- c(rep("control", ncol(Control)), rep("case", ncol(case)))
design <- model.matrix(~0 + factor(c(rep(1, ncol(Control)), rep(2, ncol(case)))))
colnames(design) <- c("Normal", "HS")
cont.matrix <- makeContrasts("HS-Normal", levels = design)
fit <- lmFit(M, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
table <- topTable(fit2, adjust.method = "BH", number = nrow(Meth), sort.by='logFC')


case$CaseMean <- rowMeans(case)
Control$ControlMean <- rowMeans(Control)

Meth$ControlMean <- case$CaseMean; Meth$ControlMean <- Control$ControlMean
Meth$deltaBeta <- case$CaseMean - Control$ControlMean

anno = data.frame(Class =  c(rep("control", 24), rep("case",24)))
NMF::aheatmap(M, annCol = anno)
########################################################################
Meth1 <- Meth[which(Meth$deltaBeta > 0.075),]
MergeData1 <- MergeData[rownames(Meth1),]

annotation <- MergeData1$UCSC_REFGENE_GROUP
annotation <- as.character(annotation)
annotation[which(nchar(annotation)<2)]<-NA ### When no attotation
annotation <- strsplit(annotation,"\\;")
annotation <- sapply(annotation, "[[", 1)
direction <- ifelse(MergeData1$Fold.change >1,"Hyper","Hypo")
foldchange <- MergeData1$Fold.change
#relation <- MergeData1$
AUC <- MergeData1$AUC
logP <- -log10(MergeData1$FDR.p.Val)
delta <- Meth1$deltaBeta
#########################################################################
Meth1$CaseMean <- NULL; Meth1$ControlMean <- NULL; Meth1$deltaBeta <- NULL
Meth1 <- as.matrix(Meth1)

####################### Heatmap plot on beta value ######################
set.seed(121)
ha = HeatmapAnnotation(df = data.frame(type = c(rep("Control", 24), rep("HS", 24))), col = list(type = c("HS" = "red4", "Control"="blue4")))
#column_tree = hclust(dist(t(meth)), method = "ward.D2")
ht_list = Heatmap(Meth1, clustering_method_columns = "ward.D2", name = "Methylation", col = colorRamp2(c(0, 0.15, 0.22, 0.30, 0.45, 0.6, 0.75, 1), c("green4","green", "gray", "violet","violetred3","firebrick","red1","red4")),
                  cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 2, show_row_names = TRUE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(AUC, name = "AUC", col = colorRamp2(c(0.75, 0.88, 1.0), c("blue4", "purple", "red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(foldchange, name = "Fold Change", col = colorRamp2(c(3.5, 2.5, 2), c("red","purple", "blue")), column_names_gp = gpar(fontsize = 8))+
  #Heatmap(direction, name = "Methylation Direction", col = c("Hyper" = "red", "Hypo" = "blue"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(delta, name = "Delta beta", col = colorRamp2(c(0, 0.1, 0.2), c("red4","gray", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(logP, name = "- log10 (P-value)", col = colorRamp2(c(0, 20, 40), c("red4","gray", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(annotation, name = "Probe Annotation", col = c("Body" = "aquamarine", "TSS200" = "red4", "TSS1500" = "violetred", "1stExon" = "black", "5'UTR" = "blue", "3'UTR" = "green"), column_names_gp = gpar(fontsize = 8))
#Heatmap(relation, name = "Probe Relationship", column_names_gp = gpar(fontsize = 8))
pdf("HS Beta plot.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom")
dev.off()


######################### Heatmap plot on M value ########################
M1 <- log2(Meth1/(1-  Meth1))

set.seed(121)
ha = HeatmapAnnotation(df = data.frame(type = c(rep("Control", 24), rep("HS", 24))), col = list(type = c("HS" = "red4", "Control"="blue4")))
#column_tree = hclust(dist(t(meth)), method = "ward.D2")
#clustering_method_columns = "single"# we can use "average", "median", "centroid" anyone from hclust 
ht_list = Heatmap(M1, clustering_method_columns = "ward.D2", name = "Methylation", col = colorRamp2(c(-5, -2.5, 0, 1), c("green4", "gray", "firebrick","red4")),
                  cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 1, show_row_names = TRUE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(AUC, name = "AUC", col = colorRamp2(c(0.75, 0.88, 1.0), c("blue4", "purple", "red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(foldchange, name = "Fold Change", col = colorRamp2(c(3.5, 2.5, 2), c("red","purple", "blue")), column_names_gp = gpar(fontsize = 8))+
  #Heatmap(direction, name = "Methylation Direction", col = c("Hyper" = "red", "Hypo" = "blue"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(delta, name = "Delta beta", col = colorRamp2(c(0, 0.1, 0.2), c("red4","gray", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(logP, name = "- log10 (P-value)", col = colorRamp2(c(0, 20, 40), c("red4","gray", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(annotation, name = "Probe Annotation", col = c("Body" = "aquamarine", "TSS200" = "red4", "TSS1500" = "violetred", "1stExon" = "black", "5'UTR" = "blue", "3'UTR" = "green"), column_names_gp = gpar(fontsize = 8))
#Heatmap(relation, name = "Probe Relationship", column_names_gp = gpar(fontsize = 8))
pdf("HS M plot.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom")
dev.off()
################################
save.image("HS_Latest.RData")
