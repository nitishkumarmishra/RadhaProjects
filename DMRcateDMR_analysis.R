## DMRcate analysis worked on R 4.3.3 (not on R 5.1)
library(DMRcate)
library(limma)
setwd("F:/OneDrive - University of Nebraska Medical Center/Radha/CP DMRcate")
#load("PAAD_Methylation_Oct-10.RData")
#rm(list=setdiff(ls(), c("PAAD.Clinical", "BMIQ.Meth")))
Radha <- read.csv("1-CP_Individual_analysis March 21 2019 Nitish DMR.txt", header = TRUE, row.names = 2, sep = "\t")

Beta <- Radha[,grep("Beta", colnames(Radha))]
colnames(Beta) <- gsub(".AVG_Beta", "", colnames(Beta))
colnames(Beta) <- gsub("Autism_", "", colnames(Beta))
colnames(Beta) <- gsub("CP_Control", "Control", colnames(Beta))
colnames(Beta) <- gsub("NSCLP_n", "N", colnames(Beta))


case <- Beta[,grep("Case", colnames(Beta))]
Control <- Beta[,grep("Control", colnames(Beta))]
Beta <- cbind(Control, case)
meth <- as.matrix(Beta)

myMs <- logit2(Beta)
myMs.noSNPs <- rmSNPandCH(as.matrix(myMs), dist=2, mafcut=0.05)


status <- c(rep("control", ncol(Control)), rep("case", ncol(case)))
design <- model.matrix(~0 + factor(c(rep(1, ncol(Control)), rep(2, ncol(case)))))
colnames(design) <- c("Normal", "NSLP")
cont.matrix <- makeContrasts("NSLP-Normal", levels = design)

cancerID <- grep("Case", colnames(myMs.noSNPs))
normalID <- grep("Control", colnames(myMs.noSNPs))

c1 <- length(cancerID)
c2 <- length(normalID)

Tumor.BMIQ <- myMs.noSNPs[,cancerID]
Normal.BMIQ <- myMs.noSNPs[, normalID]

groups <- as.factor(c(rep("Case",c1),rep("Control",c2)))


design<-model.matrix(~0+groups)
colnames(design)=levels(groups)
contrast.matrix <- makeContrasts(Case-Control, levels = design)


myannotation <- cpg.annotate("array",myMs.noSNPs, analysis.type = "differential", design = design, contrasts = TRUE, what="M", fdr = 0.01, cont.matrix = contrast.matrix, coef = "Case - Control")
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2, p.adjust.method = "BH", pcutoff = "fdr", betacutoff=0)
#### Annotate overlapping promoter regions (+/- 2000 bp from TSS)
results.ranges <- extractRanges(dmrcoutput, genome = "hg38")
groups <- c(Case="magenta", Control="forestgreen")
type <- as.factor(c(rep("Case", c1), rep("Control", c2)))
cols <- groups[as.character(type)]
samps <- c(1:23, 23+(1:21))


DMR.plot(ranges=results.ranges, dmr=2, CpGs=meth, phen.col=cols, genome="hg38", arraytype = "450K", what="Beta", samps=samps, showSampleNames = FALSE, col.axis="Black", col.title="Black",cex.sampleNames = 0.6, col.sampleNames="Black", 
         separator = 1.75, pch=18, toscale=TRUE, plotmedians=TRUE,legend = TRUE, cex = 0.6, fontcolor="Black", fontcolor.group="Black", fontsize.group=8, fontcolor.feature = 1, cex.feature = 0.6, min.width = 3, min.distance = 5, cex.title=0.5, rotation.title=360)

## Plot is not working. To many samples and to many CpGs
par(mfrow=c(1,1))


#DMR.plot(ranges=results.ranges, dmr=797, CpGs=BMIQ.Meth, phen.col=cols, genome="hg38", arraytype = "450K", what="Beta", samps=samps, showSampleNames = TRUE, cex.sampleNames = 0.8, separator = 1, pch=18, toscale=TRUE, plotmedians=TRUE,legend = TRUE)

rm(c1, c2, numNAs, myNorm, Meth, tx.hg19, tx.hg38, tx.mm10, probe.features, probeInfoALL.lv)


df1 <- data.frame(seqnames=seqnames(results.ranges),
                  starts=start(results.ranges),
                  ends=end(results.ranges),
                  names=c(rep(".", length(results.ranges))),
                  strand=strand(results.ranges)
)
df <- mcols(results.ranges)
df2 <- cbind(df1, df)
write.csv(df2, file = "DMRcateDMR.txt")

save.image("DMRcate_analysis")

###################################


if (!dir.exists("DMRcate_Plot")) dir.create("DMRcate_Plot")
for (i in 1:100) {
  
  DMR.plot(ranges=results.ranges, dmr=i, CpGs=meth, phen.col=cols, genome="hg38", arraytype = "450K", what="Beta", samps=samps, showSampleNames = FALSE, col.axis="Black", col.title="Black",cex.sampleNames = 0.6, col.sampleNames="Black", 
           separator = 1.75, pch=18, toscale=TRUE, plotmedians=TRUE,legend = TRUE, cex = 0.6, fontcolor="Black", fontcolor.group="Black", fontsize.group=8, fontcolor.feature = 1, cex.feature = 0.6, min.width = 3, min.distance = 5, cex.title=0.5, rotation.title=360)
  
  file_name = paste0("DMRcate_Plot", "/DMR_", i, ".pdf")
  dev.print(pdf, file_name, width = 9, height = 9)
}
