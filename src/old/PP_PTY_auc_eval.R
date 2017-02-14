# This script is to calculate AUC for Protein-protein and Pathway database.
# 4704 node_name means I used all genes with more than 5 interactions in the following analysis.
# I also tries 1777 node, which used genes with more than 71 intereactions.
# Final results are very similar, but 4704 gave a little bit higher mean AUC numbers.(0.1,0.2 higher)

library(ROCR)

load("4704_nodeName.RData")
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


length(which(node_name %in% colnames(rc_pcc)))
filter_node_name <- node_name[which(node_name %in% colnames(rc_pcc))]
#OR
filter_node_name <- colnames(rc_scc_ft)

cogene_result <- cogene_result[names(cogene_result) %in% filter_node_name]


#vst
##########################################################################################################################
vst_gcc_ft <- vst_gcc[,which(colnames(vst_gcc) %in% filter_node_name)] # dim vst_gcc_ft 15116*753                        #
auc_vst_gcc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=vst_gcc_ft);rm(vst_gcc,vst_gcc_ft)                   #
                                                                                                                         #
vst_pcc_ft <- vst_pcc[,which(colnames(vst_pcc) %in% filter_node_name)]                                                   #
auc_vst_pcc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=vst_pcc_ft);rm(vst_pcc,vst_pcc_ft)                   #
                                                                                                                         #
vst_scc_ft <- vst_scc[,which(colnames(vst_scc) %in% filter_node_name)]                                                   #
auc_vst_scc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=vst_scc_ft);rm(vst_scc,vst_scc_ft)                   #
                                                                                                                         #
vst_kcc_ft <- vst_kcc[,which(colnames(vst_kcc) %in% filter_node_name)]                                                   #
auc_vst_kcc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=vst_kcc_ft);rm(vst_kcc,vst_kcc_ft)                   #
                                                                                                                         #
vst_bicor_ft <- vst_bicor[,which(colnames(vst_bicor) %in% filter_node_name)]                                             #
auc_vst_bicor <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=vst_bicor_ft);rm(vst_bicor,vst_bicor_ft)           #
                                                                                                                         #
vst_aa_ft <- vst_aa[,which(colnames(vst_aa) %in% filter_node_name)]                                       #
auc_vst_aa <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=vst_aa_ft);rm(vst_aa,vst_aa_ft)   #
                                                                                                                         #
vst_ma_ft <- vst_ma[,which(colnames(vst_ma) %in% filter_node_name)]                                       #
auc_vst_ma <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=vst_ma_ft);rm(vst_ma,vst_ma_ft)   #
                                                                                                                         #
vst_mrnet_ft <- vst_mrnet[,which(colnames(vst_mrnet) %in% filter_node_name)]                                             #
auc_vst_mrnet <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=vst_mrnet_ft);rm(vst_mrnet,vst_mrnet_ft)           #
                                                                                                                         #
vst_clr_ft <- vst_clr[,which(colnames(vst_clr) %in% filter_node_name)]                                                   #
auc_vst_clr <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=vst_clr_ft);rm(vst_clr,vst_clr_ft)                   #
                                                                                                                         #
method <- c("gcc","pcc","scc","kcc","bicor","aa","marance","mrnet","clr")                                           #
save(paste0("auc_vst_",method),file="AUC_vst.RData")                                                                     #
##########################################################################################################################

#cpm
##########################################################################################################################
cpm_gcc_ft <- cpm_gcc[,which(colnames(cpm_gcc) %in% filter_node_name)] # dim cpm_gcc_ft 15116*753                        #
auc_cpm_gcc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm_gcc_ft);rm(cpm_gcc,cpm_gcc_ft)                   #
                                                                                                                         #
cpm_pcc_ft <- cpm_pcc[,which(colnames(cpm_pcc) %in% filter_node_name)]                                                   #
auc_cpm_pcc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm_pcc_ft);rm(cpm_pcc,cpm_pcc_ft)                   #
                                                                                                                         #
cpm_scc_ft <- cpm_scc[,which(colnames(cpm_scc) %in% filter_node_name)]                                                   #
auc_cpm_scc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm_scc_ft);rm(cpm_scc,cpm_scc_ft)                   #
                                                                                                                         #
cpm_kcc_ft <- cpm_kcc[,which(colnames(cpm_kcc) %in% filter_node_name)]                                                   #
auc_cpm_kcc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm_kcc_ft);rm(cpm_kcc,cpm_kcc_ft)                   #
                                                                                                                         #
cpm_bicor_ft <- cpm_bicor[,which(colnames(cpm_bicor) %in% filter_node_name)]                                             #
auc_cpm_bicor <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm_bicor_ft);rm(cpm_bicor,cpm_bicor_ft)           #
                                                                                                                         #
cpm_aa_ft <- cpm_aa[,which(colnames(cpm_aa) %in% filter_node_name)]                                       #
auc_cpm_aa <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm_aa_ft);rm(cpm_aa,cpm_aa_ft)   #
                                                                                                                         #
cpm_ma_ft <- cpm_ma[,which(colnames(cpm_ma) %in% filter_node_name)]                                       #
auc_cpm_ma <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm_ma_ft);rm(cpm_ma,cpm_ma_ft)   #
                                                                                                                         #
cpm_mrnet_ft <- cpm_mrnet[,which(colnames(cpm_mrnet) %in% filter_node_name)]                                             #
auc_cpm_mrnet <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm_mrnet_ft);rm(cpm_mrnet,cpm_mrnet_ft)           #
                                                                                                                         #
cpm_clr_ft <- cpm_clr[,which(colnames(cpm_clr) %in% filter_node_name)]                                                   #
auc_cpm_clr <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm_clr_ft);rm(cpm_clr,cpm_clr_ft)                   #
                                                                                                                         #
method <- c("gcc","pcc","scc","kcc","bicor","aa","marance","mrnet","clr")                                           #
save(list=c(paste0("auc_cpm_",method)),file="AUC_cpm.RData")                                                                     #
##########################################################################################################################


#rc
##########################################################################################################################
rc_gcc_ft <- rc_gcc[,which(colnames(rc_gcc) %in% filter_node_name)] # dim rc_gcc_ft 15116*753                            #
auc_rc_gcc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=rc_gcc_ft);rm(rc_gcc,rc_gcc_ft)                       #
                                                                                                                         #
rc_pcc_ft <- rc_pcc[,which(colnames(rc_pcc) %in% filter_node_name)]                                                      #
auc_rc_pcc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=rc_pcc_ft);rm(rc_pcc,rc_pcc_ft)                       #
                                                                                                                         #
rc_scc_ft <- rc_scc[,which(colnames(rc_scc) %in% filter_node_name)]                                                      #
auc_rc_scc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=rc_scc_ft);rm(rc_scc,rc_scc_ft)                       #
                                                                                                                         #
rc_kcc_ft <- rc_kcc[,which(colnames(rc_kcc) %in% filter_node_name)]                                                      #
auc_rc_kcc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=rc_kcc_ft);rm(rc_kcc,rc_kcc_ft)                       #
                                                                                                                         #
rc_bicor_ft <- rc_bicor[,which(colnames(rc_bicor) %in% filter_node_name)]                                                #
auc_rc_bicor <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=rc_bicor_ft);rm(rc_bicor,rc_bicor_ft)               #
                                                                                                                         #
rc_aa_ft <- rc_aa[,which(colnames(rc_aa) %in% filter_node_name)]                                          #
auc_rc_aa <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=rc_aa_ft);rm(rc_aa,rc_aa_ft)       #
                                                                                                                         #
rc_ma_ft <- rc_ma[,which(colnames(rc_ma) %in% filter_node_name)]                                          #
auc_rc_ma <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=rc_ma_ft);rm(rc_ma,rc_ma_ft)       #
                                                                                                                         #
rc_mrnet_ft <- rc_mrnet[,which(colnames(rc_mrnet) %in% filter_node_name)]                                                #
auc_rc_mrnet <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=rc_mrnet_ft);rm(rc_mrnet,rc_mrnet_ft)               #
                                                                                                                         #
rc_clr_ft <- rc_clr[,which(colnames(rc_clr) %in% filter_node_name)]                                                      #
auc_rc_clr <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=rc_clr_ft);rm(rc_clr,rc_clr_ft)                       #
                                                                                                                         #
method <- c("gcc","pcc","scc","kcc","bicor","aa","marance","mrnet","clr")                                           #
save(list=c(paste0("auc_rc_",method)),file="AUC_rc.RData")                                                                       #
##########################################################################################################################


#rpkm
##############################################################################################################################
rpkm_gcc_ft <- rpkm_gcc[,which(colnames(rpkm_gcc) %in% filter_node_name)] # dim rpkm_gcc_ft 15116*753                        #
auc_rpkm_gcc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=rpkm_gcc_ft);rm(rpkm_gcc,rpkm_gcc_ft)                   #
                                                                                                                             #
rpkm_pcc_ft <- rpkm_pcc[,which(colnames(rpkm_pcc) %in% filter_node_name)]                                                    #
auc_rpkm_pcc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=rpkm_pcc_ft);rm(rpkm_pcc,rpkm_pcc_ft)                   #
                                                                                                                             #
rpkm_scc_ft <- rpkm_scc[,which(colnames(rpkm_scc) %in% filter_node_name)]                                                    #
auc_rpkm_scc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=rpkm_scc_ft);rm(rpkm_scc,rpkm_scc_ft)                   #
                                                                                                                             #
rpkm_kcc_ft <- rpkm_kcc[,which(colnames(rpkm_kcc) %in% filter_node_name)]                                                    #
auc_rpkm_kcc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=rpkm_kcc_ft);rm(rpkm_kcc,rpkm_kcc_ft)                   #
                                                                                                                             #
rpkm_bicor_ft <- rpkm_bicor[,which(colnames(rpkm_bicor) %in% filter_node_name)]                                              #
auc_rpkm_bicor <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=rpkm_bicor_ft);rm(rpkm_bicor,rpkm_bicor_ft)           #
                                                                                                                             #
rpkm_aa_ft <- rpkm_aa[,which(colnames(rpkm_aa) %in% filter_node_name)]                                        #
auc_rpkm_aa <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=rpkm_aa_ft);rm(rpkm_aa,rpkm_aa_ft)   #
                                                                                                                             #
rpkm_ma_ft <- rpkm_ma[,which(colnames(rpkm_ma) %in% filter_node_name)]                                        #
auc_rpkm_ma <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=rpkm_ma_ft);rm(rpkm_ma,rpkm_ma_ft)   #
                                                                                                                             #
rpkm_mrnet_ft <- rpkm_mrnet[,which(colnames(rpkm_mrnet) %in% filter_node_name)]                                              #
auc_rpkm_mrnet <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=rpkm_mrnet_ft);rm(rpkm_mrnet,rpkm_mrnet_ft)           #
                                                                                                                             #
rpkm_clr_ft <- rpkm_clr[,which(colnames(rpkm_clr) %in% filter_node_name)]                                                    #
auc_rpkm_clr <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=rpkm_clr_ft);rm(rpkm_clr,rpkm_clr_ft)                   #
                                                                                                                             #
method <- c("gcc","pcc","scc","kcc","bicor","aa","marance","mrnet","clr")                                               #
save(list=c(paste0("auc_rpkm_",method)),file="AUC_rpkm.RData")                                                                       #
##############################################################################################################################


# organize result. combine all results together. MEAN.
# the previous results are lists, need to unlist() first then calcualate mean value.
auc_vst_all_mean <- c(mean(unlist(auc_vst_pcc)),mean(unlist(auc_vst_scc)),mean(unlist(auc_vst_gcc)),mean(unlist(auc_vst_kcc)),mean(unlist(auc_vst_bicor)),
                      mean(unlist(auc_vst_aa)),mean(unlist(auc_vst_ma)),mean(unlist(auc_vst_mrnet)),mean(unlist(auc_vst_clr)))

auc_cpm_all_mean <- c(mean(unlist(auc_cpm_pcc)),mean(unlist(auc_cpm_scc)),mean(unlist(auc_cpm_gcc)),mean(unlist(auc_cpm_kcc)),mean(unlist(auc_cpm_bicor)),
                      mean(unlist(auc_cpm_aa)),mean(unlist(auc_cpm_ma)),mean(unlist(auc_cpm_mrnet)),mean(unlist(auc_cpm_clr)))


auc_rc_all_mean <- c(mean(unlist(auc_rc_pcc)),mean(unlist(auc_rc_scc)),mean(unlist(auc_rc_gcc)),mean(unlist(auc_rc_kcc)),mean(unlist(auc_rc_bicor)),
                      mean(unlist(auc_rc_aa)),mean(unlist(auc_rc_ma)),mean(unlist(auc_rc_mrnet)),mean(unlist(auc_rc_clr)))

auc_rpkm_all_mean <- c(mean(unlist(auc_rpkm_pcc)),mean(unlist(auc_rpkm_scc)),mean(unlist(auc_rpkm_gcc)),mean(unlist(auc_rpkm_kcc)),mean(unlist(auc_rpkm_bicor)),
                      mean(unlist(auc_rpkm_aa)),mean(unlist(auc_rpkm_ma)),mean(unlist(auc_rpkm_mrnet)),mean(unlist(auc_rpkm_clr)))

save(auc_vst_all_mean,auc_cpm_all_mean,auc_rc_all_mean,auc_rpkm_all_mean,file="auc_mean_all_4077.RData")


# sort cogene_result by length of the elements in the list. just curious, not necessary.
len <- sapply(cogene_result, length)
cogene_result <- cogene_result[order(len,decreasing = T)]
