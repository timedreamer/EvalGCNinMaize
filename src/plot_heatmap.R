library(pheatmap)
library("corrplot")

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

t2 <- read.table("maize_rpkm_1266_allGOMatrix.txt",sep="\t",header=T)
row.names(t2) <- t2[,1]
t2 <- t2[,-1]
t <- stack(t2)
pairwise.wilcox.test(t[,1],t[,2],p.adjust.method = "b",alternative="g")


fd <- "C:\\Users\\jhuang\\OneDrive\\Network\\FROM3077Windows\\Figures\\final_AUC\\arabidopsis"

setwd(fd)

t2 <- read.table("arab_GOMatrix.txt",sep="\t",header = T)
wilcox.test(t2$X1363pcc,t2$X1363clr,paired = T,alternative = 't',conf.int = T,exact = T)
apply(X = t2[,-1],2, mean)
colnames(t2)
t2[1:5,1:5]




get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}


ggplot(corData, aes(x=Var1, y=Var2, fill=value)) +geom_tile() + xlab("") + ylab("") + 
  scale_fill_gradient2(low = "navyblue", high = "red", mid = "white",midpoint = 0.5, limit = c(0,1), space = "Lab") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) # export as pd
