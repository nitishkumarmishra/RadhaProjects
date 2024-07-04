### R code to read HTSeq-count files
setwd("C:/Users/nitis/Desktop/HTseqReadCounts")
library(data.table)
Sample_list <-  list.files(path = ".",pattern=".tab")

for (file in Sample_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.table(file, header=TRUE, sep="\t",skip = 3, stringsAsFactors = FALSE)
    colnames(dataset) <- c("gene_id","Counts")
    dataset <- subset(dataset,select = c("gene_id","Counts"))
    setnames(dataset,old = c("Counts"),new =(file))
  }
  
  # if the merged dataset does exist, append to it
  else {
    temp_dataset <- read.table(file, header=TRUE, sep="\t",skip = 3, stringsAsFactors = FALSE)
    colnames(temp_dataset) <- c("gene_id","Counts")
    temp_dataset <- subset(temp_dataset,select = c("gene_id","Counts"))
    setnames(temp_dataset,old = c("Counts"),new =(file))
    dataset<-merge(dataset, temp_dataset,by = 'gene_id')
    rm(temp_dataset)
  }
  
}

colnames(dataset) <- gsub("_R_ReadsPerGene.out.tab", "", colnames(dataset))
colnames(dataset) <- gsub("_-_|_", "-", colnames(dataset))

rownames(dataset) <- dataset$`gene-id`
dataset$`gene-id` <- NULL
colnames(dataset) <- substr(colnames(dataset), 1, 8)

##############################################
dataset1 <- dataset[, -grep("OBBM", colnames(dataset))]
dataset1 <- dataset1[apply(dataset1,1,function(x) sum(x==0))<ncol(dataset1)*0.75,]
dataset1 <- dataset1[apply(dataset1,2,function(x) sum(x==0))<nrow(dataset1)*0.75,]

PhenoData <- read.csv("PhenoData_NAS_Vs_Normal.txt", header = TRUE, sep = ",")

library(dplyr)
Ctl_Vs_NorCtl <- 
   PhenoData %>% filter(Status %in% c("Control","NormalControl"))

############### DESeq2 analysis ########
condition <- factor(Ctl_Vs_NorCtl$Condition)
dataset1 <- dataset1[, as.character(Ctl_Vs_NorCtl$Pheno)]

library(DESeq2)
countdata <- dataset1
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds

# Run the DESeq pipeline
dds <- DESeq(dds)

# Get differential expression results
res <- results(dds)
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=FALSE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
resdata <- na.omit(resdata) ## Remove genes with Pvalue = NA
## Write results
write.csv(resdata, file="NormalCtl_Vs_Ctr_diffexpr-results.csv")

resdata_0.05 <- resdata[(abs(resdata$log2FoldChange)> 1 & resdata$padj <= 0.05) ,]
write.csv(resdata_0.05, file="NormalCtl_Vs_Ctr_diffexpr_log2FC_1.0_FDR_0.05.csv")

resdata_0.01 <- resdata[(abs(resdata$log2FoldChange)> 1.5 & resdata$padj <= 0.01) ,]
write.csv(resdata_0.01, file="NormalCtl_Vs_Ctr_diffexpr_log2FC_1.5_FDR_0.01.csv")

## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")


# Plot dispersions
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
#rld <- rlogTransformation(dds)
rld <- varianceStabilizingTransformation(dds)
head(assay(rld))
hist(assay(rld))

# Colors for plots below
## Ugly:
## (mycols <- 1:length(unique(condition)))
## Use RColorBrewer, better
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
png("qc-heatmap-sample_NormalCtl_Vs_Ctr_diffexpr.png", w=1000, h=1000, pointsize=20)
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
png("diffexpr-maplot_NormalCtl_Vs_Ctr_diffexpr.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()
##################################################
## Volcano plot with "significant" genes labeled

volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot_NormalCtl_Vs_Ctr_diffexpr.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()

##############################################
rm(file, Sample_list)
#save(dataset, file = "HTseqRedCount.rda")
save.image("NormalCtl_Vs_Ctr_diffexpr.RData")

