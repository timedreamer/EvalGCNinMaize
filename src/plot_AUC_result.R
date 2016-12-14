# This scipt is to plot the result for AUC. The results have been organized in seperate files
# in the result folder.

setwd("C:\\WORK\\EvalGCNinMaize\\results")

###############################################################################################
##                           1266 Libs Total(Barplot)
###############################################################################################
#PP_PTY
result_pp <- read.delim("fouMethod_PPPTY_AUC.txt",sep="\t",header=T)
rownames(result_pp) <- result_pp[,1]
result_pp[,1] <- NULL

layout(rbind(1,2), heights=c(7,1))  # put legend on bottom 1/8th of the chart

par(mar=c(4,5,3,1))
barplot(as.matrix(result_pp),beside = T,col = c("darkred", "orange1", "yellow1","green2"),
        ylim=c(0,0.6),lwd=2.5,cex.axis =1.5,cex.lab=1.5,cex.names = 1.8)
abline(h=0.53,lty=3,col="red",lwd=2)
par(mar=c(0, 0, 0, 0))
plot.new()
legend("center",c("vst","cpm","rc","rpkm"),bty="o",ncol=4,
       fill=c("darkred", "orange1", "yellow1","green2"),cex=1.2) # export as pdf 10*9




# Gene Ontology
result_go <- read.delim("fouMethod_GO_AUC.txt",sep="\t",header=T)
rownames(result_go) <- result_go[,1]
result_go[,1] <- NULL

layout(rbind(1,2), heights=c(7,1))  # put legend on bottom 1/8th of the chart

par(mar=c(4,5,3,1))
barplot(as.matrix(result_go),beside = T,col = c("darkred", "orange1", "yellow1","green2"),
        ylim=c(0,0.8),lwd=2.5,cex.axis =1.5,cex.lab=1.5,cex.names = 1.8)
abline(h=0.58,lty=3,col="red",lwd=2)
par(mar=c(0, 0, 0, 0))
plot.new()
legend("center",c("vst","cpm","rc","rpkm"),bty="o",ncol=4,
       fill=c("darkred", "orange1", "yellow1","green2"),cex=1.2) # export as pdf 10*9





# HDA101
result_hda <- read.delim("fouMethod_HDA101_AUC.txt",sep="\t",header=T)
rownames(result_hda) <- result_hda[,1]
result_hda[,1] <- NULL

layout(rbind(1,2), heights=c(7,1))  # put legend on bottom 1/8th of the chart

par(mar=c(4,5,3,1))
barplot(as.matrix(result_hda),beside = T,col = c("darkred", "orange1", "yellow1","green2"),
        ylim=c(0,0.68),lwd=2.5,cex.axis =1.5,cex.lab=1.5,cex.names = 1.8)
abline(h=0.56,lty=3,col="red",lwd=2)
par(mar=c(0, 0, 0, 0))
plot.new()
legend("center",c("vst","cpm","rc","rpkm"),bty="o",ncol=4,
       fill=c("darkred", "orange1", "yellow1","green2"),cex=1.2) # export as pdf 10*9




###############################################################################################
##                          Plot Aggregate Network(Line chart)
###############################################################################################

# Gene Ontology
result_go <- read.delim("result_GO_AUC.txt",sep="\t",header=T)
colnames(result_go)
pcc <- result_go$pcc
mrnet <- result_go$mrnet
clr <- result_go$clr

par(mar=c(6,6,3,1))
plot(pcc,type="b",col="black",xaxt="n",ylim=c(0.5,0.76), main="",
     ylab="mean AUC",xlab="network",lwd=2.5,cex.axis =1.5,cex.lab=1.5)
lines(mrnet,type="b",col="red",lwd=2.5)
lines(clr,type="b",col="blue",lwd=2.5)
axis(1, at=1:10, labels=c(12,36,65,108,270,404,1266,"six","all","protein"),
     cex.axis=1.5,cex.lab=2)
legend("topleft",c("pcc","mrnet","clr"),lty=c(1,1,1),
       col=c("black","red","blue"),lwd=2,cex=1.5) 
# export pdf as 10*9

# PP_PTY
result_pp <- read.delim("result_PP_PTY_AUC.txt",sep="\t",header=T)
pcc <- result_pp$pcc
mrnet <- result_pp$mrnet
clr <- result_pp$clr

par(mar=c(6,6,3,1))
plot(pcc,type="b",col="black",xaxt="n",ylim=c(0.5,0.63), main="",
     ylab="mean AUC",xlab="network",lwd=2.5,cex.axis =1.5,cex.lab=1.5)
lines(mrnet,type="b",col="red",lwd=2.5)
lines(clr,type="b",col="blue",lwd=2.5)
axis(1, at=1:10, labels=c(12,36,65,108,270,404,1266,"six","all","protein"),
     cex.axis=1.5,cex.lab=2)
legend("topleft",c("pcc","mrnet","clr"),lty=c(1,1,1),
       col=c("black","red","blue"),lwd=2,cex=1.5)
# export pdf as 10*9



