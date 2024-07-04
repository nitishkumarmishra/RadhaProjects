setwd("~/Desktop")
Radha <- read.csv("10-NSCLP Individual 638 targets march 26 2018 methyl to Nitish.txt", header = TRUE, row.names = 2, sep = "\t")
Radha$Index <- NULL
Beta <- Radha[,grep("Beta", colnames(Radha))]
colnames(Beta) <- gsub("AVG_", "", colnames(Beta))
colnames(Beta) <- gsub("NSCLP_n", "N", colnames(Beta))
colnames(Beta) <- gsub(".Beta", "", colnames(Beta))
#colnames(Beta) <- gsub("Affected\\d+.","", colnames(Beta))
colnames(Beta) <- gsub("NCLP|NSCP", "NSCLP", colnames(Beta))
############################
library(limma)
library(ComplexHeatmap)
library(circlize)

case <- Beta[,grep("Affected", colnames(Beta))]
Control <- Beta[,grep("Normal", colnames(Beta))]
Beta <- cbind(Control, case)
M=log2(Beta/(1-  Beta))
status <- c(rep("control", ncol(Control)), rep("case", ncol(case)))
design <- model.matrix(~0 + factor(c(rep(1, ncol(Control)), rep(2, ncol(case)))))
colnames(design) <- c("Normal", "NSLP")
cont.matrix <- makeContrasts("NSLP-Normal", levels = design)
fit <- lmFit(Beta, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
table <- topTable(fit2, adjust.method = "BH", number = nrow(Beta), sort.by='logFC')

anno = data.frame(Class =  c(rep("control", 21), rep("case",21)))
NMF::aheatmap(Beta[1:100,], annCol = anno)
heatmap(as.matrix(Beta), scale = "row")

CaseMean <- rowMeans(case, na.rm = TRUE); ControlMean <- rowMeans(Control, na.rm = TRUE)
Delta <- as.data.frame(cbind(CaseMean, ControlMean)) 
Delta$Delta <- Delta$CaseMean - Delta$ControlMean
Delta$FC <- Delta$CaseMean/Delta$ControlMean
Delta$log2FC <- log2(Delta$CaseMean/Delta$ControlMean)

library(dplyr)
draw_rownames <- function(.data) .data %>% do(mutate(.,rownames=rownames(.)))
tmp <- Delta %>% draw_rownames() %>% arrange(desc(Delta)) %>% slice(1:20)
#tmp <- Delta %>% draw_rownames() %>% arrange(desc(FC)) %>% slice(1:25)
#Delta.0.4 <- Delta %>% draw_rownames() %>% filter(Delta > 0.4)
#Delta.0.4 <- Delta[abs(Delta$Delta) >= 0.4,]
#Beta.1 <- as.matrix(Beta)
#heatmap(Beta.1[rownames(Delta.0.4),], scale = "row")

##########################################################
Beta1 <- Beta[tmp$rownames,]
Beta1[!is.finite(as.matrix(Beta1))] <- 0 ## Replace Inf/-inf with "0"
Beta1 = t(scale(t(Beta1)))
table1 <- table[rownames(Beta1),]

annotation <- Radha[rownames(Beta1), c("UCSC_REFGENE_GROUP")]
annotation <- as.character(annotation)
annotation[which(nchar(annotation)<2)]<-NA ### When no attotation
annotation <- strsplit(annotation,"\\;")
annotation <- sapply(annotation, "[[", 1)

direction <- ifelse(tmp$log2FC >0,"Hyper","Hypo")
foldchange <- tmp$log2FC
#relation <- MergeData1$

Radha1 <- read.csv("NSCLP-Table Genes 638 targets March 26 2018 to Nitish.txt", header = TRUE, row.names = 2, sep = "\t")
AUC <- Radha1[tmp$rownames, "AUC"]
#AUC <- MergeData1$AUC
logP <- -log10(table1$adj.P.Val)
delta <- tmp$Delta

####################### Heatmap plot on beta value ######################
set.seed(121)
ha = HeatmapAnnotation(df = data.frame(Samples = c(rep("Control", 21), rep("NSCLP", 21))), col = list(Samples = c("NSCLP" = "red4", "Control"="blue4")))
#column_tree = hclust(dist(t(meth)), method = "ward.D2")
ht_list = Heatmap(Beta1, clustering_method_columns = "ward.D2", name = "Methylation", col = colorRamp2(c(-2.0, -1.0, 0, 1, 2, 3), c("green4","green","red1","red2","firebrick","red4")),
                  cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 2, show_row_names = TRUE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(AUC, name = "AUC", col = colorRamp2(c(0.75, 0.88, 1.0), c("blue4", "purple", "red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(foldchange, name = "log2(FC)", col = colorRamp2(c(5, 2, 1, 0), c("red4","red","purple", "blue")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(direction, name = "Direction", col = c("Hyper" = "red", "Hypo" = "blue"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(delta, name = "Delta Beta", col = colorRamp2(c(0.5, 0.25, 0), c("red4","red", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(logP, name = "-log10 (FDR)", col = colorRamp2(c(0, 10, 20), c("red4","red", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(annotation, name = "Annotation", col = c("Body" = "aquamarine", "TSS200" = "red4", "TSS1500" = "violetred", "1stExon" = "black", "5'UTR" = "blue", "3'UTR" = "green"), column_names_gp = gpar(fontsize = 8))
#Heatmap(relation, name = "Probe Relationship", column_names_gp = gpar(fontsize = 8))
pdf("NSCLP plot 20 delta beta.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

save.image("NSLP.RData")
