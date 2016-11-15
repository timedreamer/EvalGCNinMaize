
setwd("D:\\Users\\jhuang\\Documents\\Co-expression\\Normalization result")


# plot all gene expression into density plot
vst <- read.delim("VST_result.txt",sep="\t",header=T)
cpm <- read.delim("CPM_log2.txt",sep="\t",header=T)
rc <- read.delim("rawCount_log2_result.txt",sep="\t",header=T)
rpkm <- read.delim("RPKM_log2.txt",sep="\t",header=T)
dim(vst)
vst[1:5,1:5]
vst_expr <- assay(vst)

par(mfrow=c(2,2))
plot(density(as.matrix(vst)),main="VST_allGeneExprDistribution",
     lwd=2.5,cex.axis =1.5,cex.lab=1.5,cex.main=1,ylab="Density",xlab="")
plot(density(as.matrix(cpm[1:1266])),main="CPM_allGeneExprDistribution",
     lwd=2.5,cex.axis =1.5,cex.lab=1.5,cex.main=1,ylab="Density",xlab="")
plot(density(as.matrix(rpkm[1:1266])),main="RPKM_allGeneExprDistribution",
     lwd=2.5,cex.axis =1.5,cex.lab=1.5,cex.main=1,ylab="Density",xlab="")
plot(density(as.matrix(rc)),main="RC_allGeneExprDistribution",
     lwd=2.5,cex.axis =1.5,cex.lab=1.5,cex.main=1,ylab="Density",xlab="")

# plot gene length distribution
par(mar=c(6,6,3,1))
plot(density(as.numeric(cpm[,1267])),main="geneLength",lwd=2.5,
     cex.axis=1.2,cex.lab=1.5,cex.main=1,ylab="Density",xlab="geneLength")

# gene average expresion
par(mfrow=c(2,2))

avg_expr_vst <- apply(X = unname(vst[1:1266]),MARGIN = 1,FUN = mean)
plot(density(avg_expr_vst),main="VST_avgExpr",
     lwd=2.5,cex.axis =1.2,cex.lab=1.5,cex.main=1,ylab="Density",xlab="Expression")

avg_expr_cpm <- apply(X = unname(cpm[1:1266]),MARGIN = 1,FUN = mean)
plot(density(avg_expr_cpm),main="CPM_avgExpr",
     lwd=2.5,cex.axis =1.2,cex.lab=1.5,cex.main=1,ylab="Density",xlab="Expression")

avg_expr_rpkm <- apply(X = unname(rpkm[1:1266]),MARGIN = 1,FUN = mean)
plot(density(avg_expr_rpkm),main="RPKM_avgExpr",
     lwd=2.5,cex.axis =1.2,cex.lab=1.5,cex.main=1,ylab="Density",xlab="Expression")

avg_expr_rc <- apply(X = unname(rc[1:1266]),MARGIN = 1,FUN = mean)
plot(density(avg_expr_rc),main="RC_avgExpr",
     lwd=2.5,cex.axis =1.2,cex.lab=1.5,cex.main=1,ylab="Density",xlab="Expression")


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
