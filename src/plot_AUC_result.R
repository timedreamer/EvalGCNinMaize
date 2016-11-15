# This scipt is to plot the result for AUC. The results have been organized in seperate files
# in the result folder.

setwd("C:\\WORK\\EvalGCNinMaize\\results")

result_go <- read.delim("result_GO_AUC.txt",sep="\t",header=T)
colnames(result_go)
pcc <- result_go$pcc
mrnet <- result_go$mrnet
clr <- result_go$clr

par(mar=c(6,6,3,1))
plot(pcc,type="b",col="black",xaxt="n",ylim=c(0.5,0.75), main="",
     ylab="mean AUC",xlab="network",lwd=2.5,cex.axis =1.5,cex.lab=1.5,cex.main=2)
lines(mrnet,type="b",col="red",lwd=2.5)
lines(clr,type="b",col="blue",lwd=2.5)
axis(1, at=1:9, labels=c(12,36,65,108,270,404,1266,"agg","protein"),
     cex.axis=1.5,cex.lab=2)
legend("topleft",c("pcc","mrnet","clr"),lty=c(1,1,1),
       col=c("black","red","blue"),lwd=2,cex=1.2) 
# export pdf as 10*9


result_pp <- read.delim("result_PP_PTY_AUC.txt",sep="\t",header=T)
pcc <- result_pp$pcc
mrnet <- result_pp$mrnet
clr <- result_pp$clr

par(mar=c(6,6,3,1))
plot(pcc,type="b",col="black",xaxt="n",ylim=c(0.5,0.6), main="",
     ylab="mean AUC",xlab="network",lwd=2.5,cex.axis =1.5,cex.lab=1.5,cex.main=2)
lines(mrnet,type="b",col="red",lwd=2.5)
lines(clr,type="b",col="blue",lwd=2.5)
axis(1, at=1:9, labels=c(12,36,65,108,270,404,1266,"agg","protein"),
     cex.axis=1.5,cex.lab=2)
legend("topleft",c("pcc","mrnet","clr"),lty=c(1,1,1),
       col=c("black","red","blue"),lwd=2,cex=1.2)
# export pdf as 10*9

