# PCA analysis for RNA-Seq samples in GCN building.
# follow this http://www-huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html

setwd("C:\\WORK\\EvalGCNinMaize\\data")

cpm_norm <- read.delim("CPM_log2.txt",sep="\t",header=T) # geneLength last column
cpm_norm[1267] <- NULL # need to delete last column
cpm_norm[1:5,1:5]

# Draw by default R plot
fit <- prcomp(t(cpm_norm),scale. = F,retx =T)
# this is how calculate the percentage of each principle.
percentVar <- round(100*fit$sdev^2/sum(fit$sdev^2),1)
plot(fit$x[,1],fit$x[,2],xlab=paste0("PC1, VarExp:", round(percentVar[1],4)),
     ylab = paste0("PC2, VarExp:", round(percentVar[2],4)))
plot(fit) # show variance explained.
# an naive way to show what those samples are.
pc1_higher100 <- names(fit$x[,1][which(fit$x[,1] > 100)])


# Better to draw using ggplot2
library(ggfortify)
setwd("D:\\Dropbox\\Project_related\\Library_SRA_info")
# the order of libs are the same with cpm_norm
libs <- read.table("ALL_FILE_INCLUDE_1266_lib_info.txt",header=T,sep="\t")

# Color-blind friendly palette. https://goo.gl/6rX2j8
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")

# Color different tissue
p <- autoplot(prcomp(t(cpm_norm),scale. = F,retx =T), data = libs, colour = 'Tissue_short')
p+scale_color_manual(values=cbPalette) + labs(x = "PC1(28%)",y="PC2(14.6%)") + 
  theme(axis.title=element_text(size="22"))

# Color different genotypes
p_genotype <- autoplot(prcomp(t(cpm_norm),scale. = F,retx =T), data = libs, colour = 'Genotype_short')
p_genotype + labs(x = "PC1(28%)",y="PC2(14.6%)") + theme(axis.title=element_text(size="22"))

# Plot no color.
p_noColor <- autoplot(prcomp(t(cpm_norm),scale. = F,retx = T))
p_noColor + labs(x = "PC1(28%)",y="PC2(14.6%)") + theme(axis.title=element_text(size="22"))


######################################################################

#PCA for PPPTY and GO evaluation METHODS
setwd("C:\\WORK\\EvalGCNinMaize\\data\\AUC_rank\\maize\\methods")
vst <- read.table("maize_vst_1266_allGOMatrix.txt",sep="\t",header=T)
cpm <- read.table("maize_cpm_1266_allGOMatrix.txt",sep="\t",header=T)
rpkm <- read.table("maize_rpkm_1266_allGOMatrix.txt",sep="\t",header=T)

vst_p <- read.table("maize_vst_1266_allPPPTYMatrix__nullDelete.txt",sep="\t",header=T)
cpm_p <- read.table("maize_cpm_1266_allPPPTYMatrix_nullDelete.txt",sep="\t",header=T)
rpkm_p <- read.table("maize_rpkm_1266_allPPPTYMatrix_nullDelete.txt",sep="\t",header=T)


pcPlotGO <- function(name) {
  row.names(name) <- name[,1]
  name <- name[,-1]
  fit <- prcomp(t(name),scale. = T,retx =T)
  percentVar <- round(100*fit$sdev^2/sum(fit$sdev^2),1)
  plot(fit$x[,1],fit$x[,2],xlab=paste0("PC1, VarExp:", round(percentVar[1],4)),
     ylab = paste0("PC2, VarExp:", round(percentVar[2],4)),pch=19,
     main=paste("Methods GO"),ylim=c(-20,20),xlim=c(-20,30))
  text(fit$x[,1], fit$x[,2], colnames(vst[,-1]), pos= 1)
}

pcPlotPPPTY <- function(name) {
  row.names(name) <- name[,1]
  name <- name[,-1]
  fit <- prcomp(t(name),scale. = T,retx =T)
  percentVar <- round(100*fit$sdev^2/sum(fit$sdev^2),1)
  plot(fit$x[,1],fit$x[,2],xlab=paste0("PC1, VarExp:", round(percentVar[1],4)),
       ylab = paste0("PC2, VarExp:", round(percentVar[2],4)),pch=19,
       main=paste("Methods PPPTY"),ylim=c(-60,40),xlim=c(-40,60))
  text(fit$x[,1], fit$x[,2], colnames(vst[,-1]), pos= 1)
}

# export as svg 700*500
par(mfrow=c(2,3))
pcPlotGO(vst);pcPlotGO(cpm);pcPlotGO(rpkm)
pcPlotPPPTY(vst_p);pcPlotPPPTY(cpm_p);pcPlotPPPTY(rpkm_p)

##########################################################################
##PCA FOR INDI samples. using CPM and four methods.
setwd("C:\\WORK\\EvalGCNinMaize\\data\\AUC_rank\\maize\\indi")

allgo <- read.table("all_GO_indi.txt",sep="\t",header = T)
method <- as.factor(all[,1])
m <- as.factor(c("pcc","scc","mrnet","clr"))

fit <- prcomp(allgo[,-c(1,2)],scale. = T,retx =T)
percentVar <- round(100*fit$sdev^2/sum(fit$sdev^2),1)
plot(fit$x[,1],fit$x[,2],xlab=paste0("PC1, VarExp:", round(percentVar[1],4)),
     ylab = paste0("PC2, VarExp:", round(percentVar[2],4)),pch=19,
     main=paste("PCA_GO"),col=method)
legend("topleft", legend =m ,col=m, pch=19)
text(fit$x[,1], fit$x[,2], allgo[,2], pos= 3,cex=0.8)

pcc <- subset(allgo,allgo$method == "pcc")
pcc$average <- apply(pcc[3:279],1,mean)
pcc$average
scc <- subset(allgo,allgo$method == "scc");scc$average <- apply(scc[3:279],1,mean)
mrnet <- subset(allgo,allgo$method == "mrnet");mrnet$average <- apply(mrnet[3:279],1,mean)
clr <- subset(allgo,allgo$method == "clr");clr$average <- apply(clr[3:279],1,mean)
plot(pcc)
lines(mrnet$average,type="b",col="red",lwd=2.5)
lines(scc$average,type="b",col="green",lwd=2.5)
lines(clr$average,type="b",col="blue",lwd=2.5)
