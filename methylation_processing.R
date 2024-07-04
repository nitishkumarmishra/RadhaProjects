setwd("C:/Users/Siddesh.Southekal/OneDrive - University of Nebraska Medical Center/nitish_radhakrishna/NAS")

library(ChAMP)
library(doParallel)

# upload raw human methylation 450K IDAT files
# you should have sample sheet - only one .csv file in the folder
myLoad <- champ.load(directory = getwd(),
                     method="ChAMP",
                     methValue="B", # Beta values
                     autoimpute=TRUE,
                     filterDetP=TRUE,
                     ProbeCutoff=0.25,
                     SampleCutoff=0.25,
                     detPcut=0.05,
                     filterNoCG=TRUE,
                     arraytype="EPIC")


setwd("C:/Users/Siddesh.Southekal/OneDrive - University of Nebraska Medical Center/nitish_radhakrishna/NAS")


myNorm.BMIQ <- champ.norm(arraytype="EPIC",plotBMIQ=TRUE, resultsDir="./BMIQ_withoutXY/")  # normalization

setwd("C:/Users/Siddesh.Southekal/OneDrive - University of Nebraska Medical Center/nitish_radhakrishna/NAS")

myNorm.PBC <- champ.norm(arraytype="EPIC",method = "PBC", resultsDir="./PBC_withoutXY/")  # normalization

setwd("C:/Users/Siddesh.Southekal/OneDrive - University of Nebraska Medical Center/nitish_radhakrishna/NAS")

save.image("NAS_withoutXY.RData")


#####################################################################################################
# Quality check after normalization
# champ.QC(beta = myNorm,
#          pheno=myLoad$pd$Sample_Plate,
#          mdsPlot=FALSE,
#          densityPlot=TRUE,
#          dendrogram=TRUE)


# SVD (singular value decomposition) plot 
# champ.SVD(beta=myNorm,pd=myLoad$pd)


# Batch effect analysis
# myCombat <- champ.runCombat(beta=myNorm,
#                      pd=myLoad$pd,
#                      variablename="Sample_Group",
#                      batchname=c("Slide"),
#                      logitTrans=TRUE)


###################################################################################
# DMP and DMR analysis

# Affected <- myCombat[,grep("PanCancerAff", colnames(myCombat))]
# Control <- myCombat[,grep("PanCancerNor", colnames(myCombat))]
# Meth <- cbind(Affected, Control)


# Affected <- myNorm[,grep("-HLH", colnames(myNorm))]
# Control <- myNorm[,!grepl("-HLH", colnames(myNorm))]
# Meth <- cbind(Affected, Control)
# 
# 
# pd <- c(rep("Disease", ncol(Affected)), rep("Control", ncol(Control)))
# design <- model.matrix(~0+pd)
# colnames(design)=c("Control", "Disease")
# contrast.matrix <- makeContrasts(Disease-Control, levels = design)
# #
# #
# myannotation <- cpg.annotate("array", Meth, analysis.type = "differential", design = design, contrasts = TRUE, what="M", fdr = 0.01, cont.matrix = contrast.matrix, coef = "Disease - Control")
# dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2, p.adjust.method = "BH", pcutoff = "fdr", betacutoff=0)
# # #### Annotate overlapping promoter regions (+/- 2000 bp from TSS)
#  results.ranges <- extractRanges(dmrcoutput, genome = "hg38")
# groups <- c(Affected="magenta", Control="forestgreen")
# type <- as.factor(c(rep("Disease", ncol(Affected)), rep("Control", ncol(Control))))
# cols <- groups[as.character(type)]
# samps <- c(1:23, 23+(1:21))

#############################################
#DMR.plot(ranges=results.ranges, dmr=2, CpGs=Meth, phen.col=cols, genome="hg38", arraytype = "EPIC", what="Beta", samps=samps, showSampleNames = FALSE, col.axis="Black", col.title="Black",cex.sampleNames = 0.6, col.sampleNames="Black", separator = 1.75, pch=18, toscale=TRUE, plotmedians=TRUE,legend = TRUE, cex = 0.6, fontcolor="Black", fontcolor.group="Black", fontsize.group=8, fontcolor.feature = 1, cex.feature = 0.6, min.width = 3, min.distance = 5, cex.title=0.5, rotation.title=360)

## Plot is not working. To many samples and to many CpGs
#par(mfrow=c(1,1))


# df1 <- data.frame(seqnames=seqnames(results.ranges),
#                   starts=start(results.ranges),
#                   ends=end(results.ranges),
#                   names=c(rep(".", length(results.ranges))),
#                   strand=strand(results.ranges))
# df <- mcols(results.ranges)
# df2 <- cbind(df1, df)
# write.csv(df2, file = "DMRcateDMR3.txt")
# 


# myDMR.Combat <- champ.DMR(beta = myNorm,pheno = myLoad$pd$Sample_Group, adjPvalDmr = 0.01,method = "DMRcate")
# myDMP.Combat <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Sample_Group, adjPVal=0.01)
# myDMP.Combat$Control_to_Disease$log2FC = log2(myDMP.Combat$Control_to_Disease$Disease_AVG/myDMP.Combat$Control_to_Disease$Control_AVG)
# 
# 
# write.csv(myDMP.Combat$Disease_to_Control, file = "DiffMeth.txt")





##############################################################################################################
# PCA PBC 


library(data.table)
library(ggplot2)
#Meth <- as.data.frame(Meth); Meth <- na.omit(Meth); rownames(Meth) <- Meth$TargetID; Meth$TargetID <- NULL
pca <- prcomp(t(myNorm.PBC), scale. = TRUE)
pd <- ifelse(grepl("noOPIOD", colnames(myNorm.PBC)), "noOPIOD", ifelse(grepl("OPIOD_No_NAS", colnames(myNorm.PBC)), "OPIOD_No_NAS", "OPIOD_NAS"))
#color <- as.factor(c(rep("noOPIOD", 32),rep("OPIOD_NAS", 32),rep("OPIOD_No_NAS", 32)))
color <- as.factor(pd)
PCi<-data.frame(pca$x,Sample=color)


ggplot(PCi,aes(x=PC1,y=PC2, col=Sample))+
  geom_point(size=4,alpha=0.9)+ #Size and alpha just for fun
  scale_color_manual(values = c("blue4","red4","green4"))+ #your colors here
  theme_classic()

ggsave("2D_PCA_withoutXY.pdf", dpi = 600, width = 5, height = 5)

################### 3D PCA plot ######################
groups <- levels(color)
#group.color <- c("noOPIOD"="blue4","OPIOD_NAS"="red4","OPIOD_No_NAS"="green4")
#colors <- c(rep("red4", 23), rep("blue4", 24))
#pca$pcolor <- colors


pca$pcolor <- pd

pca$pcolor <- ifelse(grepl("noOPIOD", pca$pcolor), "blue4", ifelse(grepl("OPIOD_No_NAS", pca$pcolor), "green4", "red4"))

s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=1.2, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), angle = 45)
s3d.coords <- s3d$xyz.convert(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100)
legend("top", inset = 0.05, bty = "n", cex = 1.2, groups, col = c("blue4", "red4","green4"), pch = 19)
dev.print(pdf, '3D_PCA_withoutXY.pdf', width = 10, height = 10)



########### 3D-PCA plot with sample name #############
s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=0.01, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), angle = 45)

s3d.coords <- s3d$xyz.convert(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80)
legend("top", inset = 0.05, bty = "n", cex = 1.2, groups, col = c("blue4", "red4","green4"), pch = 19)
text(s3d$xyz.convert(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80),labels = rownames(pca$x), cex= 0.6, col = pca$pcolor)
dev.print(pdf, '3D_PCA_withoutXY_withName.pdf', width = 12, height = 12)

############################################################################################################
# Heatmap plot
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))

myDMP1 <- myDMP$Disease_to_Control[abs(myDMP$Disease_to_Control$deltaBeta) >=0.2,]
annotation <- myDMP1$feature
foldchange <- myDMP1$foldChange
relation <- myDMP1$cgi
logP <- -log10(myDMP1$adj.P.Val)
direction <- ifelse(myDMP1$foldChange >1,"Hyper","Hypo")
#delta <- abs(rowMeans(BetaCORTA.affected) - rowMeans(BetaCORTA.normal))
delta <- abs(myDMP1$deltaBeta)

Meth1 <- myNorm[rownames(myDMP1),]


set.seed(1110)
ha = HeatmapAnnotation(df = data.frame(type = c(rep("HLHS", 23), rep("Control", 24))), col = list(type = c("HLHS" = "red4", "Control"="blue4")))
#column_tree = hclust(dist(t(meth)), method = "ward.D2")
ht_list = Heatmap(Meth1, clustering_method_columns = "ward.D2", name = "Methylation", col = colorRamp2(c(0, 0.2, 0.4, 0.5, 0.8, 1), c("green4","lightseagreen", "violet","violetred3","firebrick","red4")),
                  cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 3, show_row_names = TRUE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  #Heatmap(AUC, name = "AUC", col = colorRamp2(c(0.75, 0.88, 1.0), c("blue4", "purple", "red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(foldchange, name = "Fold Change", col = colorRamp2(c(2.5, 1, 0), c("red","purple", "blue")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(direction, name = "Methylation Direction", col = c("Hyper" = "red", "Hypo" = "blue"), column_names_gp = gpar(fontsize = 8))+
  Heatmap(delta, name = "Delta beta", col = colorRamp2(c(0, 0.1, 0.3), c("gray", "red", "red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(logP, name = "log P-value", col = colorRamp2(c(2, 4, 6), c("red4","gray", "blue4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(annotation, name = "Probe Annotation", column_names_gp = gpar(fontsize = 8))
#Heatmap(relation, name = "Probe Relationship", column_names_gp = gpar(fontsize = 8))
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom")

pdf("HLHS Delta 0.2 .pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential methylation analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom")
dev.off()

