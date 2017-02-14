# agg by ranking maize. KIN3077 Windows

rm(list=ls())

setwd("D:\\Users\\jhuang\\Documents\\Co-expression\\aggregate_ntwk")

###PCC

# load data
load("./12/12_cpm_pcc.RData");load("./36/36_cpm_pcc.RData");load("./65/65_cpm_pcc.RData")
load("./108/108_cpm_pcc.RData");load("270/270_cpm_pcc.RData");load("404/404_cpm_pcc.RData")
load("Exp7/exp7_cpm_pcc.RData");load("Exp8/exp8_cpm_pcc.RData");load("Exp8/exp8_cpm_pcc.RData")
load("Exp9/exp9_cpm_pcc.RData");load("Exp10/exp10_cpm_pcc.RData");load("Exp11/exp11_cpm_pcc.RData")
load("Exp12/exp12_cpm_pcc.RData");load("Exp13/exp13_cpm_pcc.RData");load("Exp13/exp13_cpm_pcc.RData")
load("Exp15/exp15_cpm_pcc.RData");
load("Exp14/exp14_cpm_pcc.RData"); # somehow the exp14_pcc has NA in the matrix.
cpmexp14_pcc[is.na(cpmexp14_pcc)] <- 0 # need to change NA to 0.

#load MRNET
load("./12/12_cpm_mrnet.RData");load("./36/36_cpm_mrnet.RData");load("./65/65_cpm_mrnet.RData")
load("./108/108_cpm_mrnet.RData");load("270/270_cpm_mrnet.RData");load("404/404_cpm_mrnet.RData")
load("Exp7/exp7_cpm_mrnet.RData");load("Exp8/exp8_cpm_mrnet.RData");load("Exp8/exp8_cpm_mrnet.RData")
load("Exp9/exp9_cpm_mrnet.RData");load("Exp10/exp10_cpm_mrnet.RData");load("Exp11/exp11_cpm_mrnet.RData")
load("Exp12/exp12_cpm_mrnet.RData");load("Exp13/exp13_cpm_mrnet.RData");load("Exp14/exp14_cpm_mrnet.RData")
load("Exp15/exp15_cpm_mrnet.RData")


# load
load("./12/12_cpm_clr.RData");load("./36/36_cpm_clr.RData");load("./65/65_cpm_clr.RData")
load("./108/108_cpm_clr.RData");load("270/270_cpm_clr.RData");load("404/404_cpm_clr.RData")
load("Exp7/exp7_cpm_clr.RData");load("Exp8/exp8_cpm_clr.RData");load("Exp8/exp8_cpm_clr.RData")
load("Exp9/exp9_cpm_clr.RData");load("Exp10/exp10_cpm_clr.RData");load("Exp11/exp11_cpm_clr.RData")
load("Exp12/exp12_cpm_clr.RData");load("Exp13/exp13_cpm_clr.RData");load("Exp14/exp14_cpm_clr.RData")
load("Exp15/exp15_cpm_clr.RData")


# rank standalized. code similar with EGAD
ntwk_ranking <- function(gene.corr,output){
  n <- nrow(gene.corr)
  net <- matrix(rank(gene.corr, na.last = "keep", ties.method = "average"), nrow = n, ncol = n)
  rownames(net) <- rownames(gene.corr)
  colnames(net) <- colnames(gene.corr)
  net <- net/max(net, na.rm = TRUE)
  diag(net) <- 1
  assign(paste0(output,"_rank"),net)
  save(list=paste0(output,"_rank"),file=paste0(output,"rank.RData"))

}

#  PCC
system.time(ntwk_ranking(cpm12_pcc,"cpm_12pcc"))
system.time(ntwk_ranking(cpm36_pcc,"cpm_36pcc"))
system.time(ntwk_ranking(cpm65_pcc,"cpm_65pcc"))
system.time(ntwk_ranking(cpm108_pcc,"cpm_108pcc"))
system.time(ntwk_ranking(cpm270_pcc,"cpm_270pcc"))
system.time(ntwk_ranking(cpm404_pcc,"cpm_404pcc"))
system.time(ntwk_ranking(cpmexp7_pcc,"cpmexp7_pcc"))
system.time(ntwk_ranking(cpmexp8_pcc,"cpmexp8_pcc"))
system.time(ntwk_ranking(cpmexp9_pcc,"cpmexp9_pcc"))
system.time(ntwk_ranking(cpmexp10_pcc,"cpmexp10_pcc"))
system.time(ntwk_ranking(cpmexp11_pcc,"cpmexp11_pcc"))
system.time(ntwk_ranking(cpmexp12_pcc,"cpmexp12_pcc"))
system.time(ntwk_ranking(cpmexp13_pcc,"cpmexp13_pcc"))
system.time(ntwk_ranking(cpmexp14_pcc,"cpmexp14_pcc"))
system.time(ntwk_ranking(cpmexp15_pcc,"cpmexp15_pcc"))

#MRNET
system.time(ntwk_ranking(cpm12_mrnet,"cpm_12mrnet"))
system.time(ntwk_ranking(cpm36_mrnet,"cpm_36mrnet"))
system.time(ntwk_ranking(cpm65_mrnet,"cpm_65mrnet"))
system.time(ntwk_ranking(cpm108_mrnet,"cpm_108mrnet"))
system.time(ntwk_ranking(cpm270_mrnet,"cpm_270mrnet"))
system.time(ntwk_ranking(cpm404_mrnet,"cpm_404mrnet"))
system.time(ntwk_ranking(cpmexp7_mrnet,"cpmexp7_mrnet"))
system.time(ntwk_ranking(cpmexp8_mrnet,"cpmexp8_mrnet"))
system.time(ntwk_ranking(cpmexp9_mrnet,"cpmexp9_mrnet"))
system.time(ntwk_ranking(cpmexp10_mrnet,"cpmexp10_mrnet"))
system.time(ntwk_ranking(cpmexp11_mrnet,"cpmexp11_mrnet"))
system.time(ntwk_ranking(cpmexp12_mrnet,"cpmexp12_mrnet"))
system.time(ntwk_ranking(cpmexp13_mrnet,"cpmexp13_mrnet"))
system.time(ntwk_ranking(cpmexp14_mrnet,"cpmexp14_mrnet"))
system.time(ntwk_ranking(cpmexp15_mrnet,"cpmexp15_mrnet"))

#CLR

system.time(ntwk_ranking(cpm12_clr,"cpm_12clr"))
system.time(ntwk_ranking(cpm36_clr,"cpm_36clr"))
system.time(ntwk_ranking(cpm65_clr,"cpm_65clr"))
system.time(ntwk_ranking(cpm108_clr,"cpm_108clr"))
system.time(ntwk_ranking(cpm270_clr,"cpm_270clr"))
system.time(ntwk_ranking(cpm404_clr,"cpm_404clr"))
system.time(ntwk_ranking(cpmexp7_clr,"cpmexp7_clr"))
system.time(ntwk_ranking(cpmexp8_clr,"cpmexp8_clr"))
system.time(ntwk_ranking(cpmexp9_clr,"cpmexp9_clr"))
system.time(ntwk_ranking(cpmexp10_clr,"cpmexp10_clr"))
system.time(ntwk_ranking(cpmexp11_clr,"cpmexp11_clr"))
system.time(ntwk_ranking(cpmexp12_clr,"cpmexp12_clr"))
system.time(ntwk_ranking(cpmexp13_clr,"cpmexp13_clr"))
system.time(ntwk_ranking(cpmexp14_clr,"cpmexp14_clr"))
system.time(ntwk_ranking(cpmexp15_clr,"cpmexp15_clr"))



# GO eval on Agg_rank
library(EGAD)
setwd("D:\\Users\\jhuang\\Documents\\Co-expression\\aggregate_ntwk\\MaizeRank")
load("D:/Users/jhuang/Documents/Co-expression/GO_eval_annotationSub_maize.RData")

file_names=as.list(dir(pattern="cpm*"))
file_names
lapply(file_names,load,.GlobalEnv) # load all files

agg_RankntwkPcc <- cpm_12pcc_rank+cpm_36pcc_rank+cpm_65pcc_rank+cpm_108pcc_rank+cpm_270pcc_rank +
  cpm_404pcc_rank + cpmexp8_pcc_rank + cpmexp9_pcc_rank + cpmexp10_pcc_rank + cpmexp11_pcc_rank +
  cpmexp12_pcc_rank + cpmexp13_pcc_rank + cpmexp14_pcc_rank + cpmexp15_pcc_rank + cpmexp7_pcc_rank

save(agg_RankntwkPcc,file="maizeAgg_Rankpcc.RData")

GO_agg_rankPccMaize <- run_GBA(agg_RankntwkPcc,annotations_sub)
save(GO_agg_rankPccMaize,file="GO_agg_rankPcc_Maize.RData")
rm(list=ls(pattern = "cpm*"))
GO_agg_rankPccMaize[[3]] # 0.75811
GO_agg_pcc_noRank[[3]] # 0.7201086


setwd("D:\\Users\\jhuang\\Documents\\Co-expression\\aggregate_ntwk\\MaizeRank\\mrnet/")
file_names=as.list(dir(pattern="cpm*"))
file_names
lapply(file_names,load,.GlobalEnv)

agg_Rankntwkmrnet <- cpm_12mrnet_rank+cpm_36mrnet_rank+cpm_65mrnet_rank+cpm_108mrnet_rank+cpm_270mrnet_rank +
  cpm_404mrnet_rank + cpmexp8_mrnet_rank + cpmexp9_mrnet_rank + cpmexp10_mrnet_rank + cpmexp11_mrnet_rank +
  cpmexp12_mrnet_rank + cpmexp13_mrnet_rank + cpmexp14_mrnet_rank + cpmexp15_mrnet_rank + cpmexp7_mrnet_rank

save(agg_Rankntwkmrnet,file="maizeAgg_Rankmrnet.RData")

GO_agg_rankMrnetMaize <- run_GBA(agg_Rankntwkmrnet,annotations_sub)
save(GO_agg_rankMrnetMaize,file="GO_agg_rankMrnet_Maize.RData")
rm(list=ls(pattern = "cpm*"))

setwd("D:\\Users\\jhuang\\Documents\\Co-expression\\aggregate_ntwk\\MaizeRank\\clr//")
file_names=as.list(dir(pattern="cpm*"))
file_names
lapply(file_names,load,.GlobalEnv)

agg_Rankntwkclr <- cpm_12clr_rank+cpm_36clr_rank+cpm_65clr_rank+cpm_108clr_rank+cpm_270clr_rank +
  cpm_404clr_rank + cpmexp8_clr_rank + cpmexp9_clr_rank + cpmexp10_clr_rank + cpmexp11_clr_rank +
  cpmexp12_clr_rank + cpmexp13_clr_rank + cpmexp14_clr_rank + cpmexp15_clr_rank + cpmexp7_clr_rank

save(agg_Rankntwkclr,file="maizeAgg_Rankclr.RData")

GO_agg_rankClrMaize <- run_GBA(agg_Rankntwkclr,annotations_sub)
save(GO_agg_rankClrMaize,file="GO_agg_rankClr_Maize.RData")
rm(list=ls(pattern = "cpm*"))


##PPPTY evaluation
library(ROCR)
setwd("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP")
load("4704_nodeName.RData") # from higher than 5 genes connected
load("4704_PP_PTY_combine_list.RData")


calc_auc <- function(corMatrix,x){
  t1 <- abs(corMatrix[,x])
  m1 <- data.frame(matrix(nrow = 15116,ncol=3))
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

load("~/Co-expression/aggregate_ntwk/MaizeRank/clr/maizeAgg_Rankclr.RData")
load("~/Co-expression/aggregate_ntwk/MaizeRank/mrnet/maizeAgg_Rankmrnet.RData")
load("~/Co-expression/aggregate_ntwk/MaizeRank/pcc/maizeAgg_Rankpcc.RData")

# may change the data here to get filter_node_name
filter_node_name <- node_name[which(node_name %in% colnames(agg_Rankntwkclr))]
cogene_result <- cogene_result[names(cogene_result) %in% filter_node_name]

agg_rank_pcc_maize <- agg_RankntwkPcc[,which(colnames(agg_RankntwkPcc) %in% filter_node_name)] #                      #
auc_aggRank_pcc_maize <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=agg_rank_pcc_maize)                   #
mean(unlist(auc_aggRank_pcc_maize)) #

agg_rank_mrnet_maize <- agg_Rankntwkmrnet[,which(colnames(agg_Rankntwkmrnet) %in% filter_node_name)] #                      #
auc_aggRank_mrnet_maize <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=agg_rank_mrnet_maize)                   #
mean(unlist(auc_aggRank_mrnet_maize)) #

agg_rank_clr_maize <- agg_Rankntwkclr[,which(colnames(agg_Rankntwkclr) %in% filter_node_name)] #                      #
auc_aggRank_clr_maize <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=agg_rank_clr_maize)                   #
mean(unlist(auc_aggRank_clr_maize))

save(auc_aggRank_pcc_maize,auc_aggRank_mrnet_maize,auc_aggRank_clr_maize,
     file="agg_rank_maize_PPPTYResult.RData")
