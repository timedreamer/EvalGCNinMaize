# This script is to calculate the AUC, both PPPTY and GO for the aggregation of all 15 experiments
# within the 1266 libs. Testeds.
# PS: exp14_cpm_pcc need to change NA to 0 before aggregation.


######################################################################################################
#### WITHOUT RANKING STEP
#####################################################################################################

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

agg_ntwk <- abs(cpm12_pcc) + abs(cpm36_pcc) + abs(cpm65_pcc) + abs(cpm108_pcc) +
  abs(cpm270_pcc) + abs(cpm404_pcc) + abs(cpmexp7_pcc) +abs(cpmexp8_pcc) +
  abs(cpmexp9_pcc) +abs(cpmexp10_pcc) +abs(cpmexp11_pcc) +abs(cpmexp12_pcc) +
  abs(cpmexp13_pcc) +abs(cpmexp15_pcc) + abs(cpmexp14_pcc)

rm(cpm12_pcc,cpm36_pcc,cpm65_pcc,cpm108_pcc,cpm270_pcc,cpm404_pcc,cpmexp7_pcc,cpmexp8_pcc,
   cpmexp9_pcc,cpmexp10_pcc,cpmexp11_pcc,cpmexp12_pcc,cpmexp13_pcc,cpmexp14_pcc,cpmexp15_pcc)

agg15_ntwkPcc <- agg_ntwk
save(agg15_ntwkPcc,file="agg15_ntwkPCC.Rdata") # save the aggregated network PCC

cpmagg_test_ft <- agg_ntwk[,which(colnames(agg_ntwk) %in% filter_node_name)]                      #
auc_agg_test <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpmagg_test_ft)                   #
mean(unlist(auc_agg_test)) #
auc_agg_pcc_noRank <- auc_agg_test
save(auc_agg_pcc_noRank,file="agg_cpm_pccAll_PPPTY_noRank.RData")


system.time(GO_pcc_agg_ntwk <- run_GBA(agg_ntwk,annotations_sub))
GO_pcc_agg_ntwk[[3]]
save(GO_pcc_agg_ntwk,file="agg_cpmPccAll_GO_noRank.RData")

rm(GO_pcc_agg_ntwk,agg_ntwk,agg_ntwk)


###MRNET

# load data
load("./12/12_cpm_mrnet.RData");load("./36/36_cpm_mrnet.RData");load("./65/65_cpm_mrnet.RData")
load("./108/108_cpm_mrnet.RData");load("270/270_cpm_mrnet.RData");load("404/404_cpm_mrnet.RData")
load("Exp7/exp7_cpm_mrnet.RData");load("Exp8/exp8_cpm_mrnet.RData");load("Exp8/exp8_cpm_mrnet.RData")
load("Exp9/exp9_cpm_mrnet.RData");load("Exp10/exp10_cpm_mrnet.RData");load("Exp11/exp11_cpm_mrnet.RData")
load("Exp12/exp12_cpm_mrnet.RData");load("Exp13/exp13_cpm_mrnet.RData");load("Exp13/exp13_cpm_mrnet.RData")
load("Exp15/exp15_cpm_mrnet.RData")

agg_ntwk <- abs(cpm12_mrnet) + abs(cpm36_mrnet) + abs(cpm65_mrnet) + abs(cpm108_mrnet) +
  abs(cpm270_mrnet) + abs(cpm404_mrnet) + abs(cpmexp7_mrnet) +abs(cpmexp8_mrnet) +
  abs(cpmexp9_mrnet) +abs(cpmexp10_mrnet) +abs(cpmexp11_mrnet) +abs(cpmexp12_mrnet) +
  abs(cpmexp13_mrnet) +abs(cpmexp15_mrnet) + abs(cpmexp14_mrnet)


rm(cpm12_mrnet,cpm36_mrnet,cpm65_mrnet,cpm108_mrnet,cpm270_mrnet,cpm404_mrnet,cpmexp7_mrnet,cpmexp8_mrnet,
   cpmexp9_mrnet,cpmexp10_mrnet,cpmexp11_mrnet,cpmexp12_mrnet,cpmexp13_mrnet,cpmexp14_mrnet,cpmexp15_mrnet)

agg15_ntwkMrnet <- agg_ntwk
save(agg15_ntwkMrnet,file="agg15_ntwkMrnet.Rdata") # save the aggregated network Mrnet

cpmagg_test_ft <- agg_ntwk[,which(colnames(agg_ntwk) %in% filter_node_name)] # dim vst_gcc_ft 15116*753                        #
auc_agg_test <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpmagg_test_ft)                   #
mean(unlist(auc_agg_test)) #
auc_agg_mrnet_noRank <- auc_agg_test
save(auc_agg_mrnet_noRank,file="agg_cpm_mrnetAll_PPPTY_noRank.RData")

GO_mrnet_agg_ntwk <- run_GBA(agg_ntwk,annotations_sub)
GO_mrnet_agg_ntwk[[3]]
GO_agg_mrnet_noRank <- GO_mrnet_agg_ntwk
save(GO_agg_mrnet_noRank,file="agg_cpmMrnetAll_GO_noRank.RData")


###CLR

# load data
load("./12/12_cpm_clr.RData");load("./36/36_cpm_clr.RData");load("./65/65_cpm_clr.RData")
load("./108/108_cpm_clr.RData");load("270/270_cpm_clr.RData");load("404/404_cpm_clr.RData")
load("Exp7/exp7_cpm_clr.RData");load("Exp8/exp8_cpm_clr.RData");load("Exp8/exp8_cpm_clr.RData")
load("Exp9/exp9_cpm_clr.RData");load("Exp10/exp10_cpm_clr.RData");load("Exp11/exp11_cpm_clr.RData")
load("Exp12/exp12_cpm_clr.RData");load("Exp13/exp13_cpm_clr.RData");load("Exp13/exp13_cpm_clr.RData")
load("Exp15/exp15_cpm_clr.RData")


agg_ntwk <- abs(cpm12_clr) + abs(cpm36_clr) + abs(cpm65_clr) + abs(cpm108_clr) +
  abs(cpm270_clr) + abs(cpm404_clr) + abs(cpmexp7_clr) +abs(cpmexp8_clr) +
  abs(cpmexp9_clr) +abs(cpmexp10_clr) +abs(cpmexp11_clr) +abs(cpmexp12_clr) +
  abs(cpmexp13_clr) +abs(cpmexp15_clr) + abs(cpmexp14_clr)


rm(cpm12_clr,cpm36_clr,cpm65_clr,cpm108_clr,cpm270_clr,cpm404_clr,cpmexp7_clr,cpmexp8_clr,
   cpmexp9_clr,cpmexp10_clr,cpmexp11_clr,cpmexp12_clr,cpmexp13_clr,cpmexp15_clr,cpmexp14_clr)

agg15_ntwkClr <- agg_ntwk
save(agg15_ntwkClr,file="agg15_ntwkClr.Rdata") # save the aggregated network Clr

cpmagg_test_ft <- agg_ntwk[,which(colnames(agg_ntwk) %in% filter_node_name)] # dim vst_gcc_ft 15116*753                        #
auc_agg_test <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpmagg_test_ft)                   #
mean(unlist(auc_agg_test)) #
auc_agg_clr_noRank <- auc_agg_test
save(auc_agg_clr_noRank,file="agg_cpm_clrAll_PPPTY_noRank.RData")

GO_clr_agg_ntwk <- run_GBA(agg_ntwk,annotations_sub)
GO_clr_agg_ntwk[[3]]
GO_agg_clr_noRank <- GO_clr_agg_ntwk
save(GO_agg_clr_noRank,file="agg_cpmClrAll_GO_noRank.RData")
