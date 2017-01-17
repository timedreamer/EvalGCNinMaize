# this file is to calculat arabidopsis aggrgated networks.
library(EGAD)

# pcc
setwd("D:/Users/jhuang/Documents/Co-expression/arabidopsis/Text_arab/arab_all_pcc")

load("exp01_cpm_pcc.RData");load("exp02_cpm_pcc.RData");load("exp03_cpm_pcc.RData");load("exp04_cpm_pcc.RData");load("exp05_cpm_pcc.RData");
load("exp06_cpm_pcc.RData");load("exp07_cpm_pcc.RData");load("exp08_cpm_pcc.RData");load("exp09_cpm_pcc.RData");load("exp10_cpm_pcc.RData");
load("exp11_cpm_pcc.RData");load("exp12_cpm_pcc.RData");load("exp13_cpm_pcc.RData");load("exp14_cpm_pcc.RData");load("exp15_cpm_pcc.RData");
load("exp16_cpm_pcc.RData");load("exp17_cpm_pcc.RData");load("exp18_cpm_pcc.RData");load("exp19_cpm_pcc.RData");load("exp20_cpm_pcc.RData");
load("exp21_cpm_pcc.RData");load("exp22_cpm_pcc.RData");load("exp23_cpm_pcc.RData");load("exp24_cpm_pcc.RData");load("exp25_cpm_pcc.RData");
load("exp26_cpm_pcc.RData");load("exp27_cpm_pcc.RData");load("exp28_cpm_pcc.RData");load("exp29_cpm_pcc.RData");load("exp30_cpm_pcc.RData");
load("exp31_cpm_pcc.RData");load("exp32_cpm_pcc.RData");load("exp33_cpm_pcc.RData");load("exp34_cpm_pcc.RData");load("exp35_cpm_pcc.RData");
load("exp36_cpm_pcc.RData");load("exp37_cpm_pcc.RData");load("exp38_cpm_pcc.RData");load("exp39_cpm_pcc.RData");load("exp40_cpm_pcc.RData");
load("exp41_cpm_pcc.RData");load("exp42_cpm_pcc.RData");load("exp43_cpm_pcc.RData")

cpmexp01_pcc[is.na(cpmexp01_pcc)] <- 0;cpmexp02_pcc[is.na(cpmexp02_pcc)] <- 0;cpmexp03_pcc[is.na(cpmexp03_pcc)] <- 0
cpmexp04_pcc[is.na(cpmexp04_pcc)] <- 0;cpmexp05_pcc[is.na(cpmexp05_pcc)] <- 0;cpmexp06_pcc[is.na(cpmexp06_pcc)] <- 0
cpmexp07_pcc[is.na(cpmexp07_pcc)] <- 0;cpmexp08_pcc[is.na(cpmexp08_pcc)] <- 0;cpmexp09_pcc[is.na(cpmexp09_pcc)] <- 0;cpmexp10_pcc[is.na(cpmexp10_pcc)] <- 0
cpmexp11_pcc[is.na(cpmexp11_pcc)] <- 0;cpmexp12_pcc[is.na(cpmexp12_pcc)] <- 0;cpmexp13_pcc[is.na(cpmexp13_pcc)] <- 0;cpmexp14_pcc[is.na(cpmexp14_pcc)] <- 0
cpmexp15_pcc[is.na(cpmexp15_pcc)] <- 0;cpmexp16_pcc[is.na(cpmexp16_pcc)] <- 0;cpmexp17_pcc[is.na(cpmexp17_pcc)] <- 0;cpmexp18_pcc[is.na(cpmexp18_pcc)] <- 0
cpmexp19_pcc[is.na(cpmexp19_pcc)] <- 0;cpmexp20_pcc[is.na(cpmexp20_pcc)] <- 0;cpmexp21_pcc[is.na(cpmexp21_pcc)] <- 0;cpmexp22_pcc[is.na(cpmexp22_pcc)] <- 0
cpmexp23_pcc[is.na(cpmexp23_pcc)] <- 0;cpmexp24_pcc[is.na(cpmexp24_pcc)] <- 0;cpmexp25_pcc[is.na(cpmexp25_pcc)] <- 0;cpmexp26_pcc[is.na(cpmexp26_pcc)] <- 0
cpmexp27_pcc[is.na(cpmexp27_pcc)] <- 0;cpmexp28_pcc[is.na(cpmexp28_pcc)] <- 0;cpmexp29_pcc[is.na(cpmexp29_pcc)] <- 0;cpmexp30_pcc[is.na(cpmexp30_pcc)] <- 0
cpmexp31_pcc[is.na(cpmexp31_pcc)] <- 0;cpmexp32_pcc[is.na(cpmexp32_pcc)] <- 0;cpmexp33_pcc[is.na(cpmexp33_pcc)] <- 0;cpmexp34_pcc[is.na(cpmexp34_pcc)] <- 0
cpmexp35_pcc[is.na(cpmexp35_pcc)] <- 0;cpmexp36_pcc[is.na(cpmexp36_pcc)] <- 0;cpmexp37_pcc[is.na(cpmexp37_pcc)] <- 0;cpmexp38_pcc[is.na(cpmexp38_pcc)] <- 0
cpmexp39_pcc[is.na(cpmexp39_pcc)] <- 0;cpmexp40_pcc[is.na(cpmexp40_pcc)] <- 0;cpmexp41_pcc[is.na(cpmexp41_pcc)] <- 0;cpmexp42_pcc[is.na(cpmexp42_pcc)] <- 0
cpmexp43_pcc[is.na(cpmexp43_pcc)] <- 0

agg_ntwk_pcc <- abs(cpmexp01_pcc) + abs(cpmexp02_pcc) +abs(cpmexp03_pcc) +abs(cpmexp04_pcc) +abs(cpmexp05_pcc) +abs(cpmexp06_pcc) +abs(cpmexp07_pcc) +
  abs(cpmexp08_pcc) +abs(cpmexp09_pcc) +abs(cpmexp10_pcc) +abs(cpmexp11_pcc) +abs(cpmexp12_pcc) +abs(cpmexp13_pcc) +abs(cpmexp14_pcc) +
  abs(cpmexp15_pcc) +abs(cpmexp16_pcc) +abs(cpmexp17_pcc) +abs(cpmexp18_pcc) +abs(cpmexp19_pcc) +abs(cpmexp20_pcc) +abs(cpmexp21_pcc) +abs(cpmexp22_pcc) +
  abs(cpmexp23_pcc) +abs(cpmexp24_pcc) +abs(cpmexp25_pcc) +abs(cpmexp26_pcc) +abs(cpmexp27_pcc) +abs(cpmexp28_pcc) +abs(cpmexp29_pcc) +
  abs(cpmexp30_pcc) +abs(cpmexp31_pcc) +abs(cpmexp32_pcc) +abs(cpmexp33_pcc) +abs(cpmexp34_pcc) +abs(cpmexp35_pcc) +abs(cpmexp36_pcc) +
  abs(cpmexp37_pcc) +abs(cpmexp38_pcc) +abs(cpmexp39_pcc) +abs(cpmexp40_pcc) +abs(cpmexp41_pcc) +abs(cpmexp42_pcc) +abs(cpmexp43_pcc)


save(agg_ntwk_pcc,file="arab_agg_ntwk_pcc.Rdata") # save the aggregated network pcc
rm(list = ls(pattern = "*_pcc"))

load("~/Co-expression/arabidopsis/Text_arab/arab_all_pcc/arab_agg_ntwk_pcc.Rdata")
GO_arab_agg_pcc <- run_GBA(agg_ntwk_pcc,annotations_sub)
GO_arab_agg_pcc[[3]]
save(GO_arab_agg_pcc,file="GO_arab_agg_pcc.RData")

#indi 
load("~/Co-expression/arabidopsis/arabidopsis_GO_pcc.RData")
GO_arab_pcc[[3]]

pcc_indi <- GO_arab_pcc[[1]][,1]
pcc_agg <- GO_arab_agg_pcc[[1]][,1]
wilcox.test(pcc_indi,pcc_agg)
t1 <- pcc_indi - pcc_agg
table(t1>0)


# mrnet
setwd("D:/Users/jhuang/Documents/Co-expression/arabidopsis/Text_arab/arab_all_mrnet")

load("exp01_cpm_mrnet.RData");load("exp02_cpm_mrnet.RData");load("exp03_cpm_mrnet.RData");load("exp04_cpm_mrnet.RData");load("exp05_cpm_mrnet.RData");
load("exp06_cpm_mrnet.RData");load("exp07_cpm_mrnet.RData");load("exp08_cpm_mrnet.RData");load("exp09_cpm_mrnet.RData");load("exp10_cpm_mrnet.RData");
load("exp11_cpm_mrnet.RData");load("exp12_cpm_mrnet.RData");load("exp13_cpm_mrnet.RData");load("exp14_cpm_mrnet.RData");load("exp15_cpm_mrnet.RData");
load("exp16_cpm_mrnet.RData");load("exp17_cpm_mrnet.RData");load("exp18_cpm_mrnet.RData");load("exp19_cpm_mrnet.RData");load("exp20_cpm_mrnet.RData");
load("exp21_cpm_mrnet.RData");load("exp22_cpm_mrnet.RData");load("exp23_cpm_mrnet.RData");load("exp24_cpm_mrnet.RData");load("exp25_cpm_mrnet.RData");
load("exp26_cpm_mrnet.RData");load("exp27_cpm_mrnet.RData");load("exp28_cpm_mrnet.RData");load("exp29_cpm_mrnet.RData");load("exp30_cpm_mrnet.RData");
load("exp31_cpm_mrnet.RData");load("exp32_cpm_mrnet.RData");load("exp33_cpm_mrnet.RData");load("exp34_cpm_mrnet.RData");load("exp35_cpm_mrnet.RData");
load("exp36_cpm_mrnet.RData");load("exp37_cpm_mrnet.RData");load("exp38_cpm_mrnet.RData");load("exp39_cpm_mrnet.RData");load("exp40_cpm_mrnet.RData");
load("exp41_cpm_mrnet.RData");load("exp42_cpm_mrnet.RData");load("exp43_cpm_mrnet.RData")


agg_ntwk_mrnet <- abs(cpmexp01_mrnet) + abs(cpmexp02_mrnet) +abs(cpmexp03_mrnet) +abs(cpmexp04_mrnet) +abs(cpmexp05_mrnet) +abs(cpmexp06_mrnet) +abs(cpmexp07_mrnet) +
  abs(cpmexp08_mrnet) +abs(cpmexp09_mrnet) +abs(cpmexp10_mrnet) +abs(cpmexp11_mrnet) +abs(cpmexp12_mrnet) +abs(cpmexp13_mrnet) +abs(cpmexp14_mrnet) +
  abs(cpmexp15_mrnet) +abs(cpmexp16_mrnet) +abs(cpmexp17_mrnet) +abs(cpmexp18_mrnet) +abs(cpmexp19_mrnet) +abs(cpmexp20_mrnet) +abs(cpmexp21_mrnet) +abs(cpmexp22_mrnet) +
  abs(cpmexp23_mrnet) +abs(cpmexp24_mrnet) +abs(cpmexp25_mrnet) +abs(cpmexp26_mrnet) +abs(cpmexp27_mrnet) +abs(cpmexp28_mrnet) +abs(cpmexp29_mrnet) +
  abs(cpmexp30_mrnet) +abs(cpmexp31_mrnet) +abs(cpmexp32_mrnet) +abs(cpmexp33_mrnet) +abs(cpmexp34_mrnet) +abs(cpmexp35_mrnet) +abs(cpmexp36_mrnet) +
  abs(cpmexp37_mrnet) +abs(cpmexp38_mrnet) +abs(cpmexp39_mrnet) +abs(cpmexp40_mrnet) +abs(cpmexp41_mrnet) +abs(cpmexp42_mrnet) +abs(cpmexp43_mrnet)



save(agg_ntwk_mrnet,file="arab_agg_ntwk_mrnet.Rdata") # save the aggregated network mrnet

rm(list = ls(pattern = "*_mrnet"))

load("~/Co-expression/arabidopsis/Text_arab/arab_all_mrnet/arab_agg_ntwk_mrnet.Rdata")
GO_arab_agg_mrnet <- run_GBA(agg_ntwk_mrnet,annotations_sub)

save(GO_arab_agg_mrnet,file="GO_arab_agg_mrnet.RData")

GO_arab_agg_pcc[[3]]
GO_arab_agg_mrnet[[3]]
GO_arab_agg_clr[[3]]

GO_arab_pcc[[3]]
GO_arab_mrnet[[3]]
GO_arab_clr[[3]]

wilcox.test(GO_arab_agg_pcc[[1]][,1], GO_arab_agg_clr[[1]][,1])

# CLR
setwd("D:/Users/jhuang/Documents/Co-expression/arabidopsis/Text_arab/arab_all_clr")

load("exp01_cpm_clr.RData");load("exp02_cpm_clr.RData");load("exp03_cpm_clr.RData");load("exp04_cpm_clr.RData");load("exp05_cpm_clr.RData");
load("exp06_cpm_clr.RData");load("exp07_cpm_clr.RData");load("exp08_cpm_clr.RData");load("exp09_cpm_clr.RData");load("exp10_cpm_clr.RData");
load("exp11_cpm_clr.RData");load("exp12_cpm_clr.RData");load("exp13_cpm_clr.RData");load("exp14_cpm_clr.RData");load("exp15_cpm_clr.RData");
load("exp16_cpm_clr.RData");load("exp17_cpm_clr.RData");load("exp18_cpm_clr.RData");load("exp19_cpm_clr.RData");load("exp20_cpm_clr.RData");
load("exp21_cpm_clr.RData");load("exp22_cpm_clr.RData");load("exp23_cpm_clr.RData");load("exp24_cpm_clr.RData");load("exp25_cpm_clr.RData");
load("exp26_cpm_clr.RData");load("exp27_cpm_clr.RData");load("exp28_cpm_clr.RData");load("exp29_cpm_clr.RData");load("exp30_cpm_clr.RData");
load("exp31_cpm_clr.RData");load("exp32_cpm_clr.RData");load("exp33_cpm_clr.RData");load("exp34_cpm_clr.RData");load("exp35_cpm_clr.RData");
load("exp36_cpm_clr.RData");load("exp37_cpm_clr.RData");load("exp38_cpm_clr.RData");load("exp39_cpm_clr.RData");load("exp40_cpm_clr.RData");
load("exp41_cpm_clr.RData");load("exp42_cpm_clr.RData");load("exp43_cpm_clr.RData")


agg_ntwk_clr <- abs(cpmexp01_clr) + abs(cpmexp02_clr) +abs(cpmexp03_clr) +abs(cpmexp04_clr) +abs(cpmexp05_clr) +abs(cpmexp06_clr) +abs(cpmexp07_clr) +
  abs(cpmexp08_clr) +abs(cpmexp09_clr) +abs(cpmexp10_clr) +abs(cpmexp11_clr) +abs(cpmexp12_clr) +abs(cpmexp13_clr) +abs(cpmexp14_clr) +
  abs(cpmexp15_clr) +abs(cpmexp16_clr) +abs(cpmexp17_clr) +abs(cpmexp18_clr) +abs(cpmexp19_clr) +abs(cpmexp20_clr) +abs(cpmexp21_clr) +abs(cpmexp22_clr) +
  abs(cpmexp23_clr) +abs(cpmexp24_clr) +abs(cpmexp25_clr) +abs(cpmexp26_clr) +abs(cpmexp27_clr) +abs(cpmexp28_clr) +abs(cpmexp29_clr) +
  abs(cpmexp30_clr) +abs(cpmexp31_clr) +abs(cpmexp32_clr) +abs(cpmexp33_clr) +abs(cpmexp34_clr) +abs(cpmexp35_clr) +abs(cpmexp36_clr) +
  abs(cpmexp37_clr) +abs(cpmexp38_clr) +abs(cpmexp39_clr) +abs(cpmexp40_clr) +abs(cpmexp41_clr) +abs(cpmexp42_clr) +abs(cpmexp43_clr)



save(agg_ntwk_clr,file="arab_agg_ntwk_clr.Rdata") # save the aggregated network Clr


rm(list = ls(pattern = "*_clr"))
GO_arab_agg_clr <- run_GBA(agg_ntwk_clr,annotations_sub)

save(GO_arab_agg_clr,file="GO_arab_agg_clr.RData")




###Protein-protein 

library(ROCR)
calc_auc <- function(corMatrix,x){
  t1 <- abs(corMatrix[,x])
  m1 <- data.frame(matrix(nrow = 13519,ncol=3))
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

filter_node_name <- node_name[which(node_name %in% colnames(agg_ntwk_pcc))]
cogene_result <- cogene_result[names(cogene_result) %in% filter_node_name]

ara_pcc_new <- (agg_ntwk_pcc - mean(agg_ntwk_pcc))/sd(agg_ntwk_pcc) 

# the pcc correlaition matrix was right shifted a little.
# this cause AUROC higher, like what happend for RC.So I did a normalization.
ara_pcc_ft <- ara_pcc_new[,which(colnames(ara_pcc_new) %in% filter_node_name)] # dim vst_gcc_ft 15116*753                        #
auc_ara_pcc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=ara_pcc_new);
mean(unlist(auc_ara_pcc));

ara_mrnet_ft <- agg_ntwk_mrnet[,which(colnames(agg_ntwk_mrnet) %in% filter_node_name)] # dim vst_gcc_ft 15116*753                        #
auc_ara_mrnet <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=agg_ntwk_mrnet);
mean(unlist(auc_ara_mrnet));

ara_clr_ft <- agg_ntwk_clr[,which(colnames(agg_ntwk_clr) %in% filter_node_name)] # dim vst_gcc_ft 15116*753                        #
auc_ara_clr <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=agg_ntwk_clr);
mean(unlist(auc_ara_clr));

save(auc_ara_pcc,auc_ara_mrnet,auc_ara_clr,file="auc_AggPP_3Methods.RData")

ap <- unlist(auc_ara_pcc);am <- unlist(auc_ara_mrnet);ac <- unlist(auc_ara_clr)

agg_ap <- unlist(auc_ara_pcc);agg_am <- unlist(auc_ara_mrnet);agg_ac <- unlist(auc_ara_clr)
mean(ac);mean(agg_ac)
wilcox.test(ac,agg_ac,paired = F,alternative = "g")

agg <- cbind(agg_ap,agg_am,agg_ac);agg_name <- c(rep("pcc",3138),rep("mrnet",3138),rep("clr",3138))
pairwise.wilcox.test(agg,agg_name,p.adjust.method = 'b')

total_1363 <- cbind(ap,am,ac);total_1363Name <- c(rep("pcc",3138),rep("mrnet",3138),rep("clr",3138))
pairwise.wilcox.test(total_1363,total_1363Name,p.adjust.method = 'b')
