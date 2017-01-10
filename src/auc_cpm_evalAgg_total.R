# This scipt is to combine all six individual network together
# then calculate AUC for the "aggregated network" from Protein-protein
# interaction dataset as well as Gene Ontology.
# run GO_eval.R first

############################################################################################################
####PART I NORMALIZED BY RANKING
########################################################################################################
# use absolute value for PCC, then convert PCC to ranking(same method used in EGAD)
cpm12_pcc <- abs(cpm12_pcc);cpm12_pcc <- cpm12_pcc/max(cpm12_pcc,na.rm=T);diag(cpm12_pcc) <- 1
cpm36_pcc <- abs(cpm36_pcc);cpm36_pcc <- cpm36_pcc/max(cpm36_pcc,na.rm=T);diag(cpm36_pcc) <- 1
cpm65_pcc <- abs(cpm65_pcc);cpm65_pcc <- cpm65_pcc/max(cpm65_pcc,na.rm=T);diag(cpm65_pcc) <- 1
cpm108_pcc <- abs(cpm108_pcc);cpm108_pcc <- cpm108_pcc/max(cpm108_pcc,na.rm=T);diag(cpm108_pcc) <- 1
cpm270_pcc <- abs(cpm270_pcc);cpm270_pcc <- cpm270_pcc/max(cpm270_pcc,na.rm=T);diag(cpm270_pcc) <- 1
cpm404_pcc <- abs(cpm404_pcc);cpm404_pcc <- cpm404_pcc/max(cpm404_pcc,na.rm=T);diag(cpm404_pcc) <- 1

# aggregated network is the sum of six individual ones
agg_ntwk <- cpm12_pcc + cpm36_pcc + cpm65_pcc + cpm108_pcc + cpm270_pcc + cpm404_pcc

# calculate AUC for PP_PTY.
# need to run calc_auc() and get filter_node_name.
cpmagg_test_ft <- agg_ntwk[,which(colnames(agg_ntwk) %in% filter_node_name)]               #
auc_agg_test <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpmagg_test_ft)                   #
mean(unlist(auc_agg_test)) # 0.5567192
save(auc_agg_test,file="agg_cpm_pcc_PPPTY.RData")

# calculate AUC for GO.
GO_pcc_agg_ntwk <- run_GBA(agg_ntwk,annotations_sub)
GO_pcc_agg_ntwk[[3]]
save(GO_pcc_agg_ntwk,file="agg_cpm_pcc_GO.RData")


##########################################################################################################
# For MRNET method. no need to calc absolute value.
cpm12_mrnet <- cpm12_mrnet/max(cpm12_mrnet,na.rm=T);diag(cpm12_mrnet) <- 1
cpm36_mrnet <- cpm36_mrnet/max(cpm36_mrnet,na.rm=T);diag(cpm36_mrnet) <- 1
cpm65_mrnet <- cpm65_mrnet/max(cpm65_mrnet,na.rm=T);diag(cpm65_mrnet) <- 1
cpm108_mrnet <- cpm108_mrnet/max(cpm108_mrnet,na.rm=T);diag(cpm108_mrnet) <- 1
cpm270_mrnet <- cpm270_mrnet/max(cpm270_mrnet,na.rm=T);diag(cpm270_mrnet) <- 1
cpm404_mrnet <- cpm404_mrnet/max(cpm404_mrnet,na.rm=T);diag(cpm404_mrnet) <- 1

agg_ntwk <- cpm12_mrnet + cpm36_mrnet + cpm65_mrnet + cpm108_mrnet + cpm270_mrnet + cpm404_mrnet

# calculate AUC PP_PTY for MRNET
cpmagg_test_ft <- agg_ntwk[,which(colnames(agg_ntwk) %in% filter_node_name)] #                         #
auc_agg_test <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpmagg_test_ft)                   #
mean(unlist(auc_agg_test)) #
save(auc_agg_test,file="agg_cpm_mrnet_PPPTY.RData")

# calculate AUC GO for MRNET
GO_mrnet_agg_ntwk <- run_GBA(agg_ntwk,annotations_sub)
GO_mrnet_agg_ntwk[[3]]
save(GO_mrnet_agg_ntwk,file="agg_cpm_mrnet_GO.RData")


############################################################################################################
# For CLR method. no need to calc absolute value.
cpm12_clr <- cpm12_clr/max(cpm12_clr,na.rm=T);diag(cpm12_clr) <- 1
cpm36_clr <- cpm36_clr/max(cpm36_clr,na.rm=T);diag(cpm36_clr) <- 1
cpm65_clr <- cpm65_clr/max(cpm65_clr,na.rm=T);diag(cpm65_clr) <- 1
cpm108_clr <- cpm108_clr/max(cpm108_clr,na.rm=T);diag(cpm108_clr) <- 1
cpm270_clr <- cpm270_clr/max(cpm270_clr,na.rm=T);diag(cpm270_clr) <- 1
cpm404_clr <- cpm404_clr/max(cpm404_clr,na.rm=T);diag(cpm404_clr) <- 1

agg_ntwk <- cpm12_clr + cpm36_clr + cpm65_clr + cpm108_clr + cpm270_clr + cpm404_clr

cpmagg_test_ft <- agg_ntwk[,which(colnames(agg_ntwk) %in% filter_node_name)] #                      #
auc_agg_test <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpmagg_test_ft)                   #
mean(unlist(auc_agg_test)) #
save(auc_agg_test,file="agg_cpm_clr_PPPTY.RData")


GO_clr_agg_ntwk <- run_GBA(agg_ntwk,annotations_sub)
GO_clr_agg_ntwk[[3]]
save(GO_clr_agg_ntwk,file="agg_cpm_clr_GO.RData")



######################################################################################################
#### WITHOUT RANKING
#####################################################################################################

#PCC
agg_ntwk <- abs(cpm12_pcc) + abs(cpm36_pcc) + abs(cpm65_pcc) + abs(cpm108_pcc) + abs(cpm270_pcc) + abs(cpm404_pcc)

cpmagg_test_ft <- agg_ntwk[,which(colnames(agg_ntwk) %in% filter_node_name)] # dim vst_gcc_ft 15116*753                        #
auc_agg_test <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpmagg_test_ft)                   #
mean(unlist(auc_agg_test)) #
auc_agg_pcc_noRank <- auc_agg_test
save(auc_agg_pcc_noRank,file="agg_cpm_pcc_PPPTY_noRank.RData")

GO_pcc_agg_ntwk <- run_GBA(agg_ntwk,annotations_sub)
GO_pcc_agg_ntwk[[3]]
GO_agg_pcc_noRank <- GO_pcc_agg_ntwk
save(GO_agg_pcc_noRank,file="agg_cpm_pcc_GO_noRank.RData")


#MRNET
agg_ntwk <- cpm12_mrnet + cpm36_mrnet + cpm65_mrnet + cpm108_mrnet + cpm270_mrnet + cpm404_mrnet

cpmagg_test_ft <- agg_ntwk[,which(colnames(agg_ntwk) %in% filter_node_name)] # dim vst_gcc_ft 15116*753                        #
auc_agg_test <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpmagg_test_ft)                   #
mean(unlist(auc_agg_test)) #
auc_agg_mrnet_noRank <- auc_agg_test
save(auc_agg_mrnet_noRank,file="agg_cpm_mrnet_PPPTY_noRank.RData")

GO_mrnet_agg_ntwk <- run_GBA(agg_ntwk,annotations_sub)
GO_mrnet_agg_ntwk[[3]]
GO_agg_mrnet_noRank <- GO_mrnet_agg_ntwk
save(GO_agg_mrnet_noRank,file="agg_cpm_mrnet_GO_noRank.RData")

#CLR
agg_ntwk <- cpm12_clr + cpm36_clr + cpm65_clr + cpm108_clr + cpm270_clr + cpm404_clr

cpmagg_test_ft <- agg_ntwk[,which(colnames(agg_ntwk) %in% filter_node_name)] # dim vst_gcc_ft 15116*753                        #
auc_agg_test <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpmagg_test_ft)                   #
mean(unlist(auc_agg_test)) #
auc_agg_clr_noRank <- auc_agg_test
save(auc_agg_clr_noRank,file="agg_cpm_clr_PPPTY_noRank.RData")

GO_clr_agg_ntwk <- run_GBA(agg_ntwk,annotations_sub)
GO_clr_agg_ntwk[[3]]
GO_agg_clr_noRank <- GO_clr_agg_ntwk
save(GO_agg_clr_noRank,file="agg_cpm_clr_GO_noRank.RData")
