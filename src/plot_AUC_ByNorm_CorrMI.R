# This script is draw AUC for different normalization method, correlation vs MI method
# individually. 

# First two figures are to show difference between correlation and MI, to show MI are better.
# Second two are to show normalization method did not change too much.
# But yes, RC is wield, if not find a good explaination, then take out.

setwd("C:\\WORK\\EvalGCNinMaize\\results")

result_pp <- read.delim("fouMethod_PPPTY_AUC.txt",sep="\t",header=T)
rownames(result_pp) <- result_pp[,1]
result_pp[,1] <- NULL
result_pp <- result_pp[-3,] # delete RC result

result_go <- read.delim("fouMethod_GO_AUC.txt",sep="\t",header=T)
rownames(result_go) <- result_go[,1]
result_go[,1] <- NULL
result_go <- result_go[-3,]

###PP_PTY

# Plot individual AUC result by Normalization method, show diff between
# Correlation and MI
# export as pdf 10*7

par(mfrow=c(1,3),mar=c(7,6,4,1))
barplot(as.matrix(result_pp[1,]),main="VST",xpd=F,las=2,
        ylim=c(0.4,0.6),lwd=2.5,cex.axis = 2,cex.names = 2,col=c(rep("grey",9)))

barplot(as.matrix(result_pp[2,]),main="CPM",xpd=F,las=2,
        ylim=c(0.4,0.6),lwd=2.5,cex.axis = 2,cex.names = 2,col=c(rep("grey",9)))

# barplot(as.matrix(result_pp[3,]),main="RC",xpd=F,
#         ylim=c(0.4,0.6),lwd=2.5,cex.axis = 1.5,cex.names = 1.2,col=c(rep("grey",9)))

barplot(as.matrix(result_pp[3,]),main="RPKM",xpd=F,las=2,
        ylim=c(0.4,0.6),lwd=2.5,cex.axis = 2,cex.names = 2,col=c(rep("grey",9)))


###GO

# Plot individual AUC result by Normalization method, show diff between
# Correlation and MI
# export as pdf 10*7

####Plot individually, Show diff in Correlation and MI
par(mfrow=c(1,3),mar=c(7,6,4,1))
barplot(as.matrix(result_go[1,]),main="VST",xpd=F,las=2,
        ylim=c(0.4,0.74),lwd=2.5,cex.axis = 2,cex.names = 2,col=c(rep("grey",9)))

barplot(as.matrix(result_go[2,]),main="CPM",xpd=F,las=2,
        ylim=c(0.4,0.74),lwd=2.5,cex.axis = 2,cex.names = 2,col=c(rep("grey",9)))

# barplot(as.matrix(result_go[3,]),main="RC",xpd=F,
#         ylim=c(0.4,0.8),lwd=2.5,cex.axis = 1.5,cex.names = 1.2,col=c(rep("grey",9)))

barplot(as.matrix(result_go[3,]),main="RPKM",xpd=F,las=2,
        ylim=c(0.4,0.74),lwd=2.5,cex.axis = 2,cex.names = 2,col=c(rep("grey",9)))





################################################################################################
###PP_PTY

# Plot individual AUC result by Correlation and MI
# show similar between normalization method
#export pdf 8*8


par(mfrow=c(3,3))
barplot(as.matrix(t(result_pp[,1])),main="PCC",names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.6),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey")

barplot(as.matrix(t(result_pp[,2])),main="SCC",names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.6),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey")

barplot(as.matrix(t(result_pp[,3])),main="GCC",names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.6),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey")

barplot(as.matrix(t(result_pp[,4])),main="KCC",names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.6),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey")

barplot(as.matrix(t(result_pp[,5])),main="Bicor",names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.6),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey")

barplot(as.matrix(t(result_pp[,6])),main="AA",names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.6),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey")

barplot(as.matrix(t(result_pp[,7])),main="MA",names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.6),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey")

barplot(as.matrix(t(result_pp[,8])),main="MRNET",names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.6),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey")

barplot(as.matrix(t(result_pp[,9])),main="CLR",names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.6),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey")

###GO  

# Plot individual AUC result by Correlation and MI
# show similar between normalization method
#export pdf 8*8

par(mfrow=c(3,3))
barplot(as.matrix(t(result_go[,1])),main="PCC",names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.71),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey")

barplot(as.matrix(t(result_go[,2])),main="SCC",names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.71),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey")

barplot(as.matrix(t(result_go[,3])),main="GCC",names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.71),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey")

barplot(as.matrix(t(result_go[,4])),main="KCC",names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.71),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey")

barplot(as.matrix(t(result_go[,5])),main="Bicor",names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.71),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey")

barplot(as.matrix(t(result_go[,6])),main="AA",names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.63),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey")

barplot(as.matrix(t(result_go[,7])),main="MA",names.arg = c("VST","CPM","RPKM"),
        xpd=F, ylim=c(0.4,0.71),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey")

barplot(as.matrix(t(result_go[,8])),main="MRNET",names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.77),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey")

barplot(as.matrix(t(result_go[,9])),main="CLR",names.arg = c("VST","CPM","RPKM"),
        xpd=F,ylim=c(0.4,0.79),lwd=2.5,cex.axis = 1.5,cex.names = 1.5,col="grey")


