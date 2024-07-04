library(gPCA)
data(caseDat)
batch<-caseDat$batch
data<-caseDat$data
out<-gPCA.batchdetect(x=data,batch=batch,center=FALSE,nperm=250)
out$delta ; out$p.val
## Plots:
gDist(out)
CumulativeVarPlot(out,ug="unguided",col="blue")
PCplot(out,ug="unguided",type="1v2")
PCplot(out,ug="unguided",type="comp",npcs=4)


out<-gPCA.batchdetect(x=data,batch=batch,center=FALSE, scaleY=FALSE,filt=NULL,nperm=1000,seed=NULL)
gDist(out)
PCplot(out,ug="guided",type="1v2")
PCplot(out,ug="guided",type="comp",npcs=3)
CumulativeVarPlot(out,ug="unguided",col="blue")

data(caseDat)
batch<-caseDat$batch
data<-caseDat$data
out<-gPCA.batchdetect(x=data,batch=batch,center=TRUE)
out$delta ; out$p.val