library(ChAMP)
setwd("C:/Users/nitish.mishra/Desktop/HS-Radha")
library("doParallel")
detectCores()
myLoad <- champ.load(directory = getwd())
myNorm <- champ.norm()
set.seed(101)
champ.SVD() ## to check batch effect

champ.SVD(beta=myNorm,pd=myLoad$pd)


myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd, batchname="Slide", variablename="Sample_Group")



#### comBat batcheffect
library(sva)
batch <- myLoad$pd$Sample_Group
modcombat = model.matrix(~1, data=myLoad$pd)
combat_edata = ComBat(dat=myNorm, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)


mod = model.matrix(~as.factor(Sample_Group), data=myLoad$pd)
