
## This notebook for batcheffect analysis
setwd("C:/Users/nitish.mishra/Desktop/HS-Radha")
library("doParallel")
detectCores()
myLoad <- champ.load(directory = getwd())
myNorm <- champ.norm()

load("Radha_HS_idat.Rdata")

t.BMIQ <- t(myNorm)
data <- as.matrix(t.BMIQ)
batch <- as.numeric(gsub("Group", "", myLoad$pd$Sample_Group))

data[1:5,1:5]

library(gPCA)
set.seed(101)
out<-gPCA.batchdetect(x=data,batch=batch,center=FALSE,nperm=250)
out$delta ; out$p.val

gDist(out)
CumulativeVarPlot(out,ug="unguided",col="blue")

PCplot(out,ug="unguided",type="1v2")
PCplot(out,ug="guided",type="comp",npcs=2)

#### Now running for comBat normalized data
# By using SVA tool
library(sva)
batch <- myLoad$pd$Sample_Group
modcombat = model.matrix(~1, data=myLoad$pd)
combat_edata = ComBat(dat=myNorm, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
data1 <- t(combat_edata)
set.seed(101)
out<-gPCA.batchdetect(x=data1,batch=batch,center=FALSE,nperm=250)
out$delta ; out$p.val

combat_edata[1:5,1:5]

gDist(out)
CumulativeVarPlot(out,ug="unguided",col="blue")

PCplot(out,ug="unguided",type="1v2")
PCplot(out,ug="guided",type="comp",npcs=2)
