library(samr)
#Beta.val <- read.table(file="1-AS-dataBeta1-51.txt", header=TRUE, sep=" ", row.names=1)
#Beta.val <- as.matrix(read.table(file="1-AS-dataBeta1-51.txt", header=TRUE, sep=" ", row.names=1, as.is=TRUE))
Beta.val <- read.table("1-AS-selected-1.csv", header = TRUE, sep = ",")
#sampleIDs <- colnames(Beta.val)
#sampleIDs <- gsub("\\.", "-", sampleIDs)
#colnames(Beta.val) <- sampleIDs
Beta.val$ProbeID_A <- NULL
Beta.val$ProbeID_B <- NULL
sampleIDs <- colnames(Beta.val)
Beta.val <- data.matrix(Beta.val)
######
#sampleIDs <- colnames(Beta.val)
AS.Samples <- grep("AS", colnames(Beta.val), value = TRUE)
Non.AS.Samples <- setdiff(sampleIDs, AS.Samples)
MethMat <- Beta.val[, c(AS.Samples, Non.AS.Samples)]
y <- c(rep(2,length(AS.Samples)),rep(1,length(Non.AS.Samples)))
data <- list(x=MethMat,y=y,logged2=FALSE,genenames=paste(row.names(MethMat)), geneid=paste(row.names(MethMat)))
samr.obj <- samr(data, resp.type = "Two class unpaired", nperms = 1000, random.seed = 123456)
delta.table <- samr.compute.delta.table(samr.obj, nvals=50)
del=0
