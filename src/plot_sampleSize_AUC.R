#


# read data
setwd("C:\\WORK\\EvalGCNinMaize\\data\\AUC_rank\\maize\\indi")

sampleSize <- read.table("maize_indi_mean.txt",sep="\t",header=T)
rownames(sampleSize) <- sampleSize[,1]

sampleSize <- sampleSize[,-1]
sampleSize <- sampleSize[1:15,]

par(mfrow=c(3,3),mar=c(7,6,4,1)+0.1)
par(mfrow=c(4,2),mar=c(4,6,4,1)+0.1)
for (i in c(2:9)){
  plot(sampleSize[,1],sampleSize[,i],xlab="",ylab="average AUROC")
  mtext("sample size",side=1,line=2.5)
  abline(lm(sampleSize[,i]~ sampleSize[,1]),col='red')
  
}

t <- cor(sampleSize[,c(3:6)])
t
t1 <- cor(sampleSize[,c(7:10)])
t1
