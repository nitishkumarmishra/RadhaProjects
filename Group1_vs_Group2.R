############### Group2 vs Group1 DESeq2 analysis ########
library(dplyr)
library(DESeq2)
##############################################
##############################################
setwd("C:/Users/nitis/Desktop/HTseqReadCounts")
library(data.table)
Sample_list <-  list.files(path = ".",pattern=".tab")

for (file in Sample_list){
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.table(file, header=TRUE, sep="\t",skip = 3, stringsAsFactors = FALSE)
    colnames(dataset) <- c("gene_id","Counts")
    dataset <- subset(dataset,select = c("gene_id","Counts"))
    setnames(dataset,old = c("Counts"),new =(file))}
  # if the merged dataset does exist, append to it
  else { 
    temp_dataset <- read.table(file, header=TRUE, sep="\t",skip = 3, stringsAsFactors = FALSE)
    colnames(temp_dataset) <- c("gene_id","Counts")
    temp_dataset <- subset(temp_dataset,select = c("gene_id","Counts"))
    setnames(temp_dataset,old = c("Counts"),new =(file))
    dataset<-merge(dataset, temp_dataset,by = 'gene_id')
    rm(temp_dataset)}}

colnames(dataset) <- gsub("_R_ReadsPerGene.out.tab", "", colnames(dataset))
colnames(dataset) <- gsub("_-_|_", "-", colnames(dataset))
rownames(dataset) <- dataset$`gene-id`
dataset$`gene-id` <- NULL
colnames(dataset) <- substr(colnames(dataset), 1, 8)

##############################################
dataset1 <- dataset[, -grep("OBBM", colnames(dataset))]
dataset1 <- dataset1[apply(dataset1,1,function(x) sum(x==0))<ncol(dataset1)*0.75,]##Remove genes which have over 25% zero
dataset1 <- dataset1[apply(dataset1,2,function(x) sum(x==0))<nrow(dataset1)*0.75,]##Remove samples which have over 25% zero
keep <- rowSums(dataset1) >= 10 ## Remove very low readcount genes.
dataset1 <- dataset1[keep,]

PhenoData <- read.csv("PhenoData_NAS_Vs_Normal.txt", header = TRUE, sep = ",")
##############################################
## Grup1 vs Group2
PhenoData <- PhenoData %>%
  filter(Group=="Group 1"| Group=="Group 2") 

condition <- factor(ifelse(PhenoData$Group=="Group 1", "exp", "control"), levels = c("control", "exp"))
dataset1 <- dataset1[, as.character(PhenoData$Pheno)]
countdata <- dataset1
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
#dds
dds <- DESeq(dds)
res <- results(dds)
table(res$padj<0.05)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=FALSE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
#head(resdata)
resdata <- na.omit(resdata) ## Remove genes with Pvalue = NA

## Write results
write.csv(resdata, file="Group1 Vs Group2 DEG.csv")
resdata_0.05 <- resdata[(abs(resdata$log2FoldChange)> 1 & resdata$padj <= 0.05) ,]
write.csv(resdata_0.05, file="Group1 Vs Group2 DEG log2FC 1.0 FDR 0.05.csv")
resdata_0.01 <- resdata[(abs(resdata$log2FoldChange)> 1.5 & resdata$padj <= 0.01) ,]
write.csv(resdata_0.01, file="Group1 Vs Group2 DEG log2FC 1.5 FDR 0.01.csv")


# Plot dispersions
png("Group1 Vs Group2 qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
#rld <- rlogTransformation(dds) ## Its very slow, so I used VST transformation
rld <- varianceStabilizingTransformation(dds)

# Colors for plots below
## Ugly:
## (mycols <- 1:length(unique(condition)))
## Use RColorBrewer, better
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
png("Group1 Vs Group2 qc-heatmap.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "red", "blue"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")

legend("topright",
       legend = unique(condition),
       col = unique(mycols[condition]), 
       lty= 1,             
       lwd = 5,           
       cex=.7)
dev.off()

#######################
## MA plot
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
## I like mine better:
library(calibrate)
maplot <- function (res, thresh=0.01, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("Group1 Vs Group2 diffexpr-maplot.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()
##################################################
## Volcano plot with "significant" genes labeled

volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, xlab="log2(Fold-Change)", ylab="log10(P-value)", main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "Both"), pch=20, col=c("red","orange","green"))
}
png("Group1 Vs Group2 diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-5, 5))
dev.off()

##############################################
rm(file, Sample_list)
save.image("Group1_vs_Group2.RData")
