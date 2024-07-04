setwd("F:/OneDrive - University of Nebraska Medical Center/Radha/Alzheimer/")
Radha <- read.csv("0-Alheimers_48_Individual_data Feb_15 2018 Original.txt", header = TRUE, row.names = 2, sep = "\t")
Radha$Index <- NULL
Beta <- Radha[,grep("Beta", colnames(Radha))]
colnames(Beta) <- gsub("Alzheimers.", "", colnames(Beta))
#colnames(Beta) <- gsub("NSCLP_n", "N", colnames(Beta))
colnames(Beta) <- gsub(".AVG_Beta", "", colnames(Beta))
#colnames(Beta) <- gsub("Affected\\d+.","", colnames(Beta))
#colnames(Beta) <- gsub("NCLP|NSCP", "NSCLP", colnames(Beta))
############################
library(limma)
library(ComplexHeatmap)
library(circlize)

case <- Beta[,grep("Aff", colnames(Beta))]
Control <- Beta[,grep("Contr", colnames(Beta))]
Beta <- cbind(case, Control)
#M=log2(Beta/(1-  Beta))
status <- c(rep("Alzheimers", ncol(case)), rep("Normal", ncol(Control)))
design <- model.matrix(~0 + factor(c(rep(1, ncol(case)), rep(2, ncol(Control)))))
colnames(design) <- c("Alzheimers", "Normal")
cont.matrix <- makeContrasts("Alzheimers-Normal", levels = design)
fit <- lmFit(Beta, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
table <- topTable(fit2, adjust.method = "BH", number = nrow(Beta), sort.by='logFC')

save.image("Alzheimers_Pvalue.RData")
