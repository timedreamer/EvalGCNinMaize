#


# read data
setwd("C:\\WORK\\EvalGCNinMaize\\data\\AUC_rank\\maize\\indi")

sampleSize <- read.table("maize_indi_mean.txt",sep="\t",header=T)
rownames(sampleSize) <- sampleSize[,1]

sampleSize <- sampleSize[,-1]
sampleSize <- sampleSize[1:15,]

# save as svg fit a logorithm model
par(mfrow=c(4,2),mar=c(4,6,4,1)+0.1)
for (i in c(2:9)){
  plot(log(sampleSize[,1]),sampleSize[,i],xlab="",ylab="average AUROC")
  mtext("log(sample size)",side=1,line=2.5)
  abline(lm(sampleSize[,i]~ log(sampleSize[,1])),col='red')
  t <- summary(lm(sampleSize[,i]~log(sampleSize[,1])))
  pv <- paste("p-value:",round(t$sigma,2));rs <- paste("rsquare:",round(t$r.squared,2));
  legend("bottomright",legend=c(rs,pv))
}

