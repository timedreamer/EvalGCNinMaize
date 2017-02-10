library(pheatmap)

setwd("D:\\Users\\jhuang\\Desktop")

t1 <- read.csv("t1.csv")
t1 <- t1[,-1]

t2 <- read.table("maize_cpm_1266_allGOMatrix.txt",sep="\t",header=T)
row.names(t2) <- t2[,1]
t2 <- t2[,-1]

# save as svg, 700*400
pheatmap(t2,scale = "row",cutree_rows = 3,show_rownames = F,
         cellwidth = 20,treeheight_row =0,treeheight_col = 20,main="cpm_GO")

t2 <- read.table("maize_cpm_1266_allPPPTYMatrix_nullDelete.txt",sep="\t",header=T)
row.names(t2) <- t2[,1]
t2 <- t2[,-1]
pheatmap(t2,scale = "row",cutree_rows = 3,show_rownames = F,fontsize_col =20,
         cellwidth = 20,treeheight_row =0,treeheight_col = 20,main="cpm_PPPTY")


cogene_result["GRMZM2G016890"]
