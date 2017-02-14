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

