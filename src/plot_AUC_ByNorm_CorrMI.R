# This script is draw AUC for different normalization method, correlation vs MI method
# individually. 

# First three figures are to show difference between correlation and MI, to show MI are better.
# Second two are to show normalization method did not change too much.
# take out RC

setwd("C:\\WORK\\EvalGCNinMaize\\results")

result_pp <- read.delim("fouMethod_PPPTY_AUC.txt",sep="\t",header=T)
rownames(result_pp) <- result_pp[,1]
result_pp[,1] <- NULL
result_pp <- result_pp[-3,] # delete RC result

result_go <- read.delim("fouMethod_GO_AUC.txt",sep="\t",header=T)
rownames(result_go) <- result_go[,1]
result_go[,1] <- NULL
result_go <- result_go[-3,]

# HDA101
result_hda <- read.delim("fouMethod_HDA101_AUC.txt",sep="\t",header=T)
rownames(result_hda) <- result_hda[,1]
result_hda[,1] <- NULL
###PP_PTY

# Plot individual AUC result by Normalization method, show diff between
# Correlation and MI
# export as pdf 10*

par(mfrow=c(3,3),mar=c(7,6,4,1)+0.1)
barplot(as.matrix(result_pp[1,]),main="VST",xpd=F,las=2,
        ylim=c(0.4,0.62),lwd=2.5,cex.axis = 2,cex.names = 2,col=c(rep("grey",9)))
abline(h=min(result_pp[1,]),lty=2);abline(h=max(result_pp[1,]),lty=3)

barplot(as.matrix(result_pp[2,]),main="CPM",xpd=F,las=2,
        ylim=c(0.4,0.62),lwd=2.5,cex.axis = 2,cex.names = 2,col=c(rep("grey",9)))
abline(h=min(result_pp[2,]),lty=2);abline(h=max(result_pp[2,]),lty=3)

barplot(as.matrix(result_pp[3,]),main="RPKM",xpd=F,las=2,
        ylim=c(0.4,0.62),lwd=2.5,cex.axis = 2,cex.names = 2,col=c(rep("grey",9)))
abline(h=min(result_pp[3,]),lty=2);abline(h=max(result_pp[3,]),lty=3)

###GO

# Plot individual AUC result by Normalization method, show diff between
# Correlation and MI
# export as svg 800*800

####Plot individually, Show diff in Correlation and MI
#par(mfrow=c(1,3),mar=c(7,6,4,1))
barplot(as.matrix(result_go[1,]),main="VST",xpd=F,las=2,
        ylim=c(0.4,0.8),lwd=2.5,cex.axis = 2,cex.names = 2,col=c(rep("grey",9)))
abline(h=min(result_go[1,]),lty=2);abline(h=max(result_go[1,]),lty=3)

barplot(as.matrix(result_go[2,]),main="CPM",xpd=F,las=2,
        ylim=c(0.4,0.8),lwd=2.5,cex.axis = 2,cex.names = 2,col=c(rep("grey",9)))
abline(h=min(result_go[2,]),lty=2);abline(h=max(result_go[2,]),lty=3)

barplot(as.matrix(result_go[3,]),main="RPKM",xpd=F,las=2,
        ylim=c(0.4,0.8),lwd=2.5,cex.axis = 2,cex.names = 2,col=c(rep("grey",9)))
abline(h=min(result_go[3,]),lty=2);abline(h=max(result_go[3,]),lty=3)

###HDA101


barplot(as.matrix(result_hda[1,]),main="VST",xpd=F,las=2,
        ylim=c(0.4,0.7),lwd=2.5,cex.axis = 2,cex.names = 2,col=c(rep("grey",9)))
abline(h=min(result_hda[1,]),lty=2);abline(h=max(result_hda[1,]),lty=3)

barplot(as.matrix(result_hda[2,]),main="CPM",xpd=F,las=2,
        ylim=c(0.4,0.7),lwd=2.5,cex.axis = 2,cex.names = 2,col=c(rep("grey",9)))
abline(h=min(result_hda[2,]),lty=2);abline(h=max(result_hda[2,]),lty=3)

barplot(as.matrix(result_hda[3,]),main="RPKM",xpd=F,las=2,
        ylim=c(0.4,0.7),lwd=2.5,cex.axis = 2,cex.names = 2,col=c(rep("grey",9)))
abline(h=min(result_hda[3,]),lty=2);abline(h=max(result_hda[3,]),lty=3)



################################################################################################
###PP_PTY

# Plot individual AUC result by Correlation and MI
# show similar between normalization method
#export svg 800*800


par(mfrow=c(3,3),mar=c(5,4,4,1)+0.1)
barplot(as.matrix(t(result_pp[,1])),main="PCC",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.65),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

barplot(as.matrix(t(result_pp[,2])),main="SCC",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.65),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

barplot(as.matrix(t(result_pp[,3])),main="GCC",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.65),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

barplot(as.matrix(t(result_pp[,4])),main="KCC",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.65),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

barplot(as.matrix(t(result_pp[,5])),main="Bicor",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.65),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

barplot(as.matrix(t(result_pp[,6])),main="AA",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.65),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

barplot(as.matrix(t(result_pp[,7])),main="MA",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.65),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

barplot(as.matrix(t(result_pp[,8])),main="MRNET",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.65),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

barplot(as.matrix(t(result_pp[,9])),main="CLR",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.65),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

###GO  

# Plot individual AUC result by Correlation and MI
# show similar between normalization method
#export svg 800*800

par(mfrow=c(3,3),mar=c(5,4,4,1)+0.1)
barplot(as.matrix(t(result_go[,1])),main="PCC",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.8),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

barplot(as.matrix(t(result_go[,2])),main="SCC",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.8),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

barplot(as.matrix(t(result_go[,3])),main="GCC",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.8),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

barplot(as.matrix(t(result_go[,4])),main="KCC",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.8),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

barplot(as.matrix(t(result_go[,5])),main="Bicor",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.8),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

barplot(as.matrix(t(result_go[,6])),main="AA",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.8),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

barplot(as.matrix(t(result_go[,7])),main="MA",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F, ylim=c(0.4,0.8),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

barplot(as.matrix(t(result_go[,8])),main="MRNET",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.8),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

barplot(as.matrix(t(result_go[,9])),main="CLR",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.8),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)


##HDA101
par(mfrow=c(3,3),mar=c(5,4,4,1)+0.1)
barplot(as.matrix(t(result_hda[,1])),main="PCC",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.7),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

barplot(as.matrix(t(result_hda[,2])),main="SCC",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.7),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

barplot(as.matrix(t(result_hda[,3])),main="GCC",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.7),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

barplot(as.matrix(t(result_hda[,4])),main="KCC",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.7),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

barplot(as.matrix(t(result_hda[,5])),main="Bicor",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.7),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

barplot(as.matrix(t(result_hda[,6])),main="AA",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.7),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

barplot(as.matrix(t(result_hda[,7])),main="MA",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F, ylim=c(0.4,0.7),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

barplot(as.matrix(t(result_hda[,8])),main="MRNET",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.7),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

barplot(as.matrix(t(result_hda[,9])),main="CLR",cex.main=2,names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.7),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey",las=2)

