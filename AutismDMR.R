setwd("F:/OneDrive - University of Nebraska Medical Center/Radha/Autism/NewAutism/DMR analysis")
suppressPackageStartupMessages(library(ChAMP))
library("doParallel")
library(DMRcate)

myCombat <- read.csv("0-Autism Blood Individual data Aug 9 2018 reanalysis.txt", header = TRUE, sep = "\t", row.names = 2)
myCombat <- myCombat[,grep(".AVG_Beta", colnames(myCombat))]
colnames(myCombat) <- gsub("X", "", colnames(myCombat))
colnames(myCombat) <- gsub(".AVG_Beta", "", colnames(myCombat))

Affected <- myCombat[,11:24]
Control <- myCombat[,1:10]
Meth <- cbind(Affected, Control)
# 
pd <- c(rep("Disease", ncol(Affected)), rep("Control", ncol(Control)))
design <- model.matrix(~0+pd)
colnames(design)=c("Control", "Disease")
contrast.matrix <- makeContrasts(Disease-Control, levels = design)

Meth <- as.matrix(Meth)
Meth1 <- impute.knn(Meth, k=15, rng.seed = 123)
Meth1 <- Meth1$data

myannotation <- cpg.annotate("array", Meth1, analysis.type = "differential", design = design, contrasts = TRUE, what="M", fdr = 0.01, cont.matrix = contrast.matrix, coef = "Disease - Control")
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2, p.adjust.method = "BH", pcutoff = "fdr", betacutoff=0)
# #### Annotate overlapping promoter regions (+/- 2000 bp from TSS)
results.ranges <- extractRanges(dmrcoutput, genome = "hg38")
groups <- c(Affected="magenta", Control="forestgreen")
type <- as.factor(c(rep("Disease", ncol(Affected)), rep("Control", ncol(Control))))
cols <- groups[as.character(type)]
samps <- c(1:14, 14+(1:10))

#############################################
DMR.plot(ranges=results.ranges, dmr=1, CpGs=Meth1, phen.col=cols, genome="hg38", arraytype = "EPIC", what="Beta", samps=samps, showSampleNames = FALSE, col.axis="Black", col.title="Black",cex.sampleNames = 0.6, col.sampleNames="Black", separator = 1.75, pch=18, toscale=TRUE, plotmedians=TRUE,legend = TRUE, cex = 0.6, fontcolor="Black", fontcolor.group="Black", fontsize.group=8, fontcolor.feature = 1, cex.feature = 0.6, min.width = 3, min.distance = 5, cex.title=0.5, rotation.title=360)

## Plot is not working. To many samples and to many CpGs
par(mfrow=c(1,1))


df1 <- data.frame(seqnames=seqnames(results.ranges),
                 starts=start(results.ranges),
                 ends=end(results.ranges),
                 names=c(rep(".", length(results.ranges))),
                 strand=strand(results.ranges)
)
df <- mcols(results.ranges)
df2 <- cbind(df1, df)
write.csv(df2, file = "DMRcateDMR.txt")


DMR <- champ.DMR(Meth1, pheno = pd, arraytype = "450k", method = "DMRcate", minProbes=3, adjPvalDmr=0.1)
DMR.Bumphunter <- champ.DMR(Meth1, pheno = pd, arraytype = "450k", method = "Bumphunter", minProbes=3, adjPvalDmr=0.05)

save.image("AutismDMR.RData")
