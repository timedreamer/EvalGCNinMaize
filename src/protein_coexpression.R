# This script is to calcualte the ROC value for the protein expression dataset,
# just published on Science journal. The protein expression was downloaded as table from
# Science website in the supplememntal materials. They also provide gene expression
# table, I calculate PCC, but did not go on. I don't think it is necessary to calculate
# gene level co-expression from their data.

# in order to run, need some functions in GO_eval.R.

#First try PCC because it's fast.
setwd("./Co-expression/")
pr_norm <- read.delim("./Protein/TableS2_Protein.txt",sep="\t",header=T) # 148 varaibles
pr_norm[1:5,1:5]
rownames(pr_norm) <- pr_norm[,1]
pr_norm[,1] <- NULL

# need to change nrows to the right dim before run. Here the table contains
# 17862 genes.
calc_auc <- function(corMatrix,x){
  t1 <- abs(corMatrix[,x])
  m1 <- data.frame(matrix(nrow = 17862,ncol=3)) # change nrows to the right dim
  m1[,1] <- as.character();m1[,2] <- as.double();m1[,3] <- 0
  m1[,1] <- names(t1);m1[,2] <- t1;m1[which(m1[,1] %in% cogene_result[[x]]),3] <- 1
  if (length(table(m1[,3]))==2){
    pred <- prediction(m1[,2],m1[,3])
    auc <- performance(pred,"auc");auc <- unlist(slot(auc, "y.values"))
  } else {
    auc <- NULL
  }
  return(auc)
}

# calculate AUC for PP_PTY
pr_norm_ft <- pr_pcc[,which(colnames(pr_pcc) %in% filter_node_name)]                       #
auc_pr_norm_pcc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=pr_norm_ft)
mean(unlist(auc_pr_norm_pcc)) # 0.565453

pr_norm_ft <- pr_mrnet[,which(colnames(pr_mrnet) %in% filter_node_name)]                       #
auc_pr_norm_mrnet <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=pr_norm_ft)
mean(unlist(auc_pr_norm_mrnet)) # 0.5366992

pr_norm_ft <- pr_clr[,which(colnames(pr_clr) %in% filter_node_name)]                       #
auc_pr_norm_clr <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=pr_norm_ft)
mean(unlist(auc_pr_norm_clr)) #0.5415833

save(auc_pr_norm_pcc,auc_pr_norm_mrnet,auc_pr_norm_clr,file="protein_auc_PPPTY.RData")

# Calculate AUC for Gene ontology
GO_pr_pcc <- run_GBA(pr_pcc,annotations_sub) # 0.599105
GO_pr_pcc[[3]]

GO_pr_mrnet <- run_GBA(pr_mrnet,annotations_sub) # 0.6364826
GO_pr_mrnet[[3]]

GO_pr_clr <- run_GBA(pr_clr,annotations_sub) # 0.6215432
GO_pr_clr[[3]]

GO_pr_mrnet <- GO_tr_mrnet;GO_pr_clr <- GO_tr_clr;
rm(GO_tr_mrnet,GO_tr_clr)
save(GO_pr_clr,GO_pr_mrnet,GO_pr_pcc,file="GO_AUROC_protein.RData")




#######################################################################################
# if want to use the gene expression they provided, the following scipt can be used.

setwd("./Co-expression/")
tr_norm <- read.delim("./Protein/TableS1_Gene.txt",sep="\t",header = T)
tr_norm[1:5,1:5]
rownames(tr_norm) <- tr_norm[,1]
tr_norm[,1] <- NULL
pr_gene <- rownames(pr_norm)

tr_norm_filter <- tr_norm[which(rownames(tr_norm) %in% pr_gene),]


tr_norm_ft <- tr_pcc[,which(colnames(tr_pcc) %in% filter_node_name)]                       #
auc_tr_norm_pcc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=tr_norm_ft)
mean(unlist(auc_tr_norm_pcc))

GO_tr_pcc <- run_GBA(tr_pcc,annotations_sub) # 0.6209497
GO_tr_pcc[[3]]
