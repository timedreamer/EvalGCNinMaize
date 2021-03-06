
setwd("D:\\Users\\jhuang\\Documents\\Co-expression\\Normalization result")
setwd("C:\\WORK\\EvalGCNinMaize\\data")

# plot all gene expression into density plot
raw_data  <- read.delim("ALL_FC_noDuplicateLib_biggerThan5Million_70allignmentRate_1266.txt",
                        stringsAsFactors = FALSE,check.names = FALSE)

vst <- read.delim("VST_result.txt",sep="\t",header=T)
cpm <- read.delim("CPM_log2.txt",sep="\t",header=T)
rc <- read.delim("rawCount_log2_result.txt",sep="\t",header=T)
rpkm <- read.delim("RPKM_log2.txt",sep="\t",header=T)
dim(vst)
vst[1:5,1:5]


par(mfrow=c(2,2),mar=c(6,6,3,1))
plot(density(as.matrix(vst)),main="",
     lwd=2.5,cex.axis =1.5,cex.lab=1.5,cex.main=1,ylab="Density",xlab="Expr");
text(18,0.18,"VST",cex=2)
curve(dnorm(x, mean=mean(as.matrix(vst)), sd=sd(as.matrix(vst))), add=TRUE,col=2,lwd=2)

plot(density(as.matrix(cpm[1:1266])),main="",
     lwd=2.5,cex.axis =1.5,cex.lab=1.5,cex.main=1,ylab="Density",xlab="Expr");
text(17,0.22,"CPM",cex=2)
curve(dnorm(x, mean=mean(as.matrix(cpm[1:1266])), sd=sd(as.matrix(cpm[1:1266]))), add=TRUE,col=2,lwd=2)

plot(density(as.matrix(rpkm[1:1266])),main="",
     lwd=2.5,cex.axis =1.5,cex.lab=1.5,cex.main=1,ylab="Density",xlab="Expr");
text(15,0.22,"RPKM",cex=2)
curve(dnorm(x, mean=mean(as.matrix(rpkm[1:1266])), sd=sd(as.matrix(rpkm[1:1266]))), add=TRUE,col=2,lwd=2)

plot(density(as.matrix(rc)),main="",
     lwd=2.5,cex.axis =1.5,cex.lab=1.5,cex.main=1,ylab="Density",xlab="Expr");
text(18,0.18,"RC",cex=2)
curve(dnorm(x, mean=mean(as.matrix(rc)), sd=sd(as.matrix(rc))), add=TRUE,col=2,lwd=2)

# plot gene length distribution
par(mar=c(6,6,3,1))
plot(density(as.numeric(cpm[,1267])),main="geneLength",lwd=2.5,
     cex.axis=1.2,cex.lab=1.5,cex.main=1,ylab="Density",xlab="geneLength")

plot(density(as.numeric(raw_data[,6])),main="geneLength",lwd=2.5,
     cex.axis=1.2,cex.lab=1.5,cex.main=1,ylab="Density",xlab="geneLength")
#export as 9*9 pdf


###############################################################################################
# gene average expresion
par(mfrow=c(2,2),mar=c(6,6,3,1))

avg_expr_vst <- apply(X = unname(vst[1:1266]),MARGIN = 1,FUN = mean)
plot(density(avg_expr_vst),main="",
     lwd=2.5,cex.axis =1.2,cex.lab=1.5,cex.main=1,ylim=c(0,0.24),ylab="Density",xlab="Expression")
text(15,0.2,"VST",cex=2)
curve(dnorm(x, mean=mean(avg_expr_vst), sd=sd(avg_expr_vst)), add=TRUE,col=2,lwd=2)

avg_expr_cpm <- apply(X = unname(cpm[1:1266]),MARGIN = 1,FUN = mean)
plot(density(avg_expr_cpm),main="",
     lwd=2.5,cex.axis =1.2,cex.lab=1.5,cex.main=1,ylim=c(0,0.29),ylab="Density",xlab="Expression")
text(11,0.25,"CPM",cex=2)
curve(dnorm(x, mean=mean(avg_expr_cpm), sd=sd(avg_expr_cpm)), add=TRUE,col=2,lwd=2)

avg_expr_rpkm <- apply(X = unname(rpkm[1:1266]),MARGIN = 1,FUN = mean)
plot(density(avg_expr_rpkm),main="",
     lwd=2.5,cex.axis =1.2,cex.lab=1.5,cex.main=1,ylab="Density",xlab="Expression")
text(10,0.28,"RPKM",cex=2)
curve(dnorm(x, mean=mean(avg_expr_rpkm), sd=sd(avg_expr_rpkm)), add=TRUE,col=2,lwd=2)

avg_expr_rc <- apply(X = unname(rc[1:1266]),MARGIN = 1,FUN = mean)
plot(density(avg_expr_rc),main="",
     lwd=2.5,cex.axis =1.2,cex.lab=1.5,cex.main=1,ylim=c(0,0.29),ylab="Density",xlab="Expression")
text(15,0.25,"RC",cex=2)
curve(dnorm(x, mean=mean(avg_expr_rc), sd=sd(avg_expr_rc)), add=TRUE,col=2,lwd=2)


# Plot gene length vs gene avg expresion
par(mfrow=c(2,2))

plot(cpm$geneLength,avg_expr_vst,
     lwd=2.5,cex.axis =1.5,cex.lab=1.5,cex.main=1,ylab="AvgExpr",xlab="Length")
text(12500,14,"VST",cex=2.5)
plot(cpm$geneLength,avg_expr_cpm,
     lwd=2.5,cex.axis =1.5,cex.lab=1.5,cex.main=1,ylab="AvgExpr",xlab="Length")
text(12500,10.5,"CPM",cex=2.5)
plot(cpm$geneLength,avg_expr_rpkm,
     lwd=2.5,cex.axis =1.5,cex.lab=1.5,cex.main=1,ylab="AvgExpr",xlab="Length")
text(11000,9.5,"RPKM",cex=2.5)
plot(cpm$geneLength,avg_expr_rc,
     lwd=2,cex.axis =1.5,cex.lab=1.5,cex.main=1,ylab="AvgExpr",xlab="Length")
text(12000,14,"RC",cex=2.5)
