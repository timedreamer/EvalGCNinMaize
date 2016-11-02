# this script is to calcualte AUC from PP_PTY for individual networks.
# six different sizes.


library(ROCR)

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


# may change the data here to get filter_node_name
filter_node_name <- node_name[which(node_name %in% colnames(cpm12_pcc))]
cogene_result <- cogene_result[names(cogene_result) %in% filter_node_name]


#12
##########################################################################################################################
cpm12_pcc_ft <- cpm12_pcc[,which(colnames(cpm12_pcc) %in% filter_node_name)]                       #
auc_12_pcc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm12_pcc_ft);rm(cpm12_pcc,cpm12_pcc_ft)                   #
#
cpm12_mrnet_ft <- cpm12_mrnet[,which(colnames(cpm12_mrnet) %in% filter_node_name)]                   #
auc_12_mrnet <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm12_mrnet_ft);rm(cpm12_mrnet,cpm12_mrnet_ft)                   #
#
cpm12_clr_ft <- cpm12_clr[,which(colnames(cpm12_clr) %in% filter_node_name)]                       #
auc_12_clr <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm12_clr_ft);rm(cpm12_clr,cpm12_clr_ft)                   #
#

method <- c("pcc","mrnet","clr")                                           #
save(list=c(paste0("auc_12_",method)),file="AUC_cpm12.RData")                                                                      #
##########################################################################################################################

#36
##########################################################################################################################
cpm36_pcc_ft <- cpm36_pcc[,which(colnames(cpm36_pcc) %in% filter_node_name)]                 #
auc_36_pcc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm36_pcc_ft);rm(cpm36_pcc,cpm36_pcc_ft)                   #
#
cpm36_mrnet_ft <- cpm36_mrnet[,which(colnames(cpm36_mrnet) %in% filter_node_name)]                     #
auc_36_mrnet <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm36_mrnet_ft);rm(cpm36_mrnet,cpm36_mrnet_ft)                   #
#
cpm36_clr_ft <- cpm36_clr[,which(colnames(cpm36_clr) %in% filter_node_name)]                       #
auc_36_clr <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm36_clr_ft);rm(cpm36_clr,cpm36_clr_ft)                   #
#

method <- c("pcc","mrnet","clr")                                           #
save(list=c(paste0("auc_36_",method)),file="AUC_cpm36.RData")                                                                        #
##########################################################################################################################

#65
##########################################################################################################################
cpm65_pcc_ft <- cpm65_pcc[,which(colnames(cpm65_pcc) %in% filter_node_name)] # dim vst_gcc_ft 15116*753                        #
auc_65_pcc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm65_pcc_ft);rm(cpm65_pcc,cpm65_pcc_ft)                   #
#
cpm65_mrnet_ft <- cpm65_mrnet[,which(colnames(cpm65_mrnet) %in% filter_node_name)] # dim vst_gcc_ft 15116*753                        #
auc_65_mrnet <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm65_mrnet_ft);rm(cpm65_mrnet,cpm65_mrnet_ft)                   #
#
cpm65_clr_ft <- cpm65_clr[,which(colnames(cpm65_clr) %in% filter_node_name)] # dim vst_gcc_ft 15116*753                        #
auc_65_clr <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm65_clr_ft);rm(cpm65_clr,cpm65_clr_ft)                   #
#

method <- c("pcc","mrnet","clr")                                           #
save(list=c(paste0("auc_65_",method)),file="AUC_cpm65.RData")                                                                        #
##########################################################################################################################

#108
##########################################################################################################################
cpm108_pcc_ft <- cpm108_pcc[,which(colnames(cpm108_pcc) %in% filter_node_name)] # dim vst_gcc_ft 15116*753                        #
auc_108_pcc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm108_pcc_ft);rm(cpm108_pcc,cpm108_pcc_ft)                   #
#
cpm108_mrnet_ft <- cpm108_mrnet[,which(colnames(cpm108_mrnet) %in% filter_node_name)] # dim vst_gcc_ft 15116*753                        #
auc_108_mrnet <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm108_mrnet_ft);rm(cpm108_mrnet,cpm108_mrnet_ft)                   #
#
cpm108_clr_ft <- cpm108_clr[,which(colnames(cpm108_clr) %in% filter_node_name)] # dim vst_gcc_ft 15116*753                        #
auc_108_clr <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm108_clr_ft);rm(cpm108_clr,cpm108_clr_ft)                   #
#

method <- c("pcc","mrnet","clr")                                           #
save(list=c(paste0("auc_108_",method)),file="AUC_cpm108.RData")                                                                     #
##########################################################################################################################

#270
##########################################################################################################################
cpm270_pcc_ft <- cpm270_pcc[,which(colnames(cpm270_pcc) %in% filter_node_name)] # dim vst_gcc_ft 15116*753                        #
auc_270_pcc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm270_pcc_ft);rm(cpm270_pcc,cpm270_pcc_ft)                   #
#
cpm270_mrnet_ft <- cpm270_mrnet[,which(colnames(cpm270_mrnet) %in% filter_node_name)] # dim vst_gcc_ft 15116*753                        #
auc_270_mrnet <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm270_mrnet_ft);rm(cpm270_mrnet,cpm270_mrnet_ft)                   #
#
cpm270_clr_ft <- cpm270_clr[,which(colnames(cpm270_clr) %in% filter_node_name)] # dim vst_gcc_ft 15116*753                        #
auc_270_clr <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm270_clr_ft);rm(cpm270_clr,cpm270_clr_ft)                   #
#

method <- c("pcc","mrnet","clr")                                           #
save(list=c(paste0("auc_270_",method)),file="AUC_cpm270.RData")                                                                    #
##########################################################################################################################


#404
##########################################################################################################################
cpm404_pcc_ft <- cpm404_pcc[,which(colnames(cpm404_pcc) %in% filter_node_name)] # dim vst_gcc_ft 15116*753                        #
auc_404_pcc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm404_pcc_ft);rm(cpm404_pcc,cpm404_pcc_ft)                   #
#
cpm404_mrnet_ft <- cpm404_mrnet[,which(colnames(cpm404_mrnet) %in% filter_node_name)] # dim vst_gcc_ft 15116*753                        #
auc_404_mrnet <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm404_mrnet_ft);rm(cpm404_mrnet,cpm404_mrnet_ft)                   #
#
cpm404_clr_ft <- cpm404_clr[,which(colnames(cpm404_clr) %in% filter_node_name)] # dim vst_gcc_ft 15116*753                        #
auc_404_clr <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=cpm404_clr_ft);rm(cpm404_clr,cpm404_clr_ft)                   #
#

method <- c("pcc","mrnet","clr")                                           #
save(list=c(paste0("auc_404_",method)),file="AUC_cpm404.RData")                                                                     #
##########################################################################################################################


#calculate mean AUC of each network and save in a file.
auc_12_mean <- c(mean(unlist(auc_12_pcc)),mean(unlist(auc_12_mrnet)),mean(unlist(auc_12_clr)))
auc_36_mean <- c(mean(unlist(auc_36_pcc)),mean(unlist(auc_36_mrnet)),mean(unlist(auc_36_clr)))
auc_65_mean <- c(mean(unlist(auc_65_pcc)),mean(unlist(auc_65_mrnet)),mean(unlist(auc_65_clr)))
auc_108_mean <- c(mean(unlist(auc_108_pcc)),mean(unlist(auc_108_mrnet)),mean(unlist(auc_108_clr)))
auc_270_mean <- c(mean(unlist(auc_270_pcc)),mean(unlist(auc_270_mrnet)),mean(unlist(auc_270_clr)))
auc_404_mean <- c(mean(unlist(auc_404_pcc)),mean(unlist(auc_404_mrnet)),mean(unlist(auc_404_clr)))


save(auc_12_mean,auc_36_mean,auc_65_mean,auc_108_mean,auc_270_mean,auc_404_mean,
  file="cpm_meanAUC_aggregate_ntwk.RData")
