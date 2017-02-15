library(pheatmap);library(ggplot2)
library(corrplot);library(reshape2)

setwd("C:\\WORK\\EvalGCNinMaize\\data\\AUC_rank\\maize\\methods")
vst <- read.table("maize_vst_1266_allGOMatrix.txt",sep="\t",header=T)
row.names(vst) <- vst[,1]
vst <- vst[,-1]

# save as svg, 700*400
pheatmap(vst,scale = "row",show_rownames = F,cutree_rows = 3,
         cellwidth = 20,treeheight_row =0,treeheight_col = 20,main="VST_GO")

t2 <- read.table("maize_vst_1266_allPPPTYMatrix__nullDelete.txt",sep="\t",header=T)
row.names(t2) <- t2[,1]
t2 <- t2[,-1]
pheatmap(t2,scale = "row",cutree_rows = 3,show_rownames = F,fontsize_col =20,
         cellwidth = 20,treeheight_row =0,treeheight_col = 20,main="RPKM_PPPTY")


# pairwise wilcoxon test
t2 <- read.table("maize_rpkm_1266_allGOMatrix.txt",sep="\t",header=T)
row.names(t2) <- t2[,1]
t2 <- t2[,-1]
t <- stack(t2)
pairwise.wilcox.test(t[,1],t[,2],p.adjust.method = "b",alternative="g")


###########################################
#plot indi heatmap
setwd("C:\\WORK\\EvalGCNinMaize\\data\\AUC_rank\\maize\\indi")
pcc <- read.table("maize_indiPCC_GO.txt",sep="\t",header=T)
row.names(pcc) <- pcc[,1]
pcc <- pcc[,-1]

pheatmap(pcc,scale='row',show_rownames = F,cutree_rows = 3,
         cellwidth = 20,treeheight_row =0,treeheight_col = 20,main="PCC_GO")

pcc <- read.table("maize_indiPCC_PPPTY.txt",sep="\t",header=T)
row.names(pcc) <- pcc[,1]
pcc <- pcc[,-1]
pheatmap(pcc,scale='row',show_rownames = F,cutree_rows = 3,
         cellwidth = 20,treeheight_row =0,treeheight_col = 20,main="PCC_GO")


pcc <- read.table("maize_indiPCC_GO_delete7.txt",sep="\t",header=T)
row.names(pcc) <- pcc[,1]
pcc <- pcc[,-1]

boxplot(pcc,main="GO_PCC",las=2,
        color=c("white"),ylim=c(0.3,1),
        outcol="grey",whiskcol="grey",staplecol="grey",cex.axis=2,
        names=c("S12_1","S12_2","S32","S36","S15","S94_2","S84","S21","S18","S24",
                "S65","S108","S270","S404","S1266","agg_15"));box(lwd=2)
points(mean_pcc,col="black",pch=8)

t <- pcc[,1:3]
boxplot(t)
colnames(t) <- factor(colnames(t),level=c("S32","S84","S94_2"))
boxplot(t)
t <- pcc[,c("S32","S84","S94_2")]

# order colname 400*600
cname <- c("S12_1","S12_2","S15","S18","S21","S24","S32","S36","S65","S84","S94_2","S108",
           "S270","S404","S1266","agg_15")
pcc <- pcc[,cname]
boxplot(pcc,las=2,color=c("white"),ylim=c(0.4,1),
        outcol="grey",whiskcol="grey",staplecol="grey",cex.axis=1.5);box(lwd=2)

pcc <- read.table("maize_indiCLR_GO_delete7.txt",sep="\t",header=T)
row.names(pcc) <- pcc[,1]
pcc <- pcc[,-1]
