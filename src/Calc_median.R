

#######################################################################
# save in D:\Users\jhuang\Documents\Co-expression\Normalization result\FromNewScript_Oct2016
# PPPTY
## VST
vst_pcc <- unlist(unname(auc_vst_pcc));vst_scc <- unlist(unname(auc_vst_scc))
vst_kcc <- unlist(unname(auc_vst_kcc));vst_gcc <- unlist(unname(auc_vst_gcc))
vst_bicor <- unlist(unname(auc_vst_bicor))
vst_aa <- unlist(unname(auc_vst_aa));vst_ma <- unlist(unname(auc_vst_ma))
vst_mrnet <- unlist(unname(auc_vst_mrnet));vst_clr <- unlist(unname(auc_vst_clr))

vst_df <- cbind(vst_pcc,vst_scc,vst_kcc,vst_gcc,vst_bicor,
                vst_aa,vst_ma,vst_mrnet,vst_clr)
vst_median_method <- apply(vst_df,2,median)

## CPM
cpm_pcc <- unlist(unname(auc_cpm_pcc));cpm_scc <- unlist(unname(auc_cpm_scc))
cpm_kcc <- unlist(unname(auc_cpm_kcc));cpm_gcc <- unlist(unname(auc_cpm_gcc))
cpm_bicor <- unlist(unname(auc_cpm_bicor))
cpm_aa <- unlist(unname(auc_cpm_aa));cpm_ma <- unlist(unname(auc_cpm_ma))
cpm_mrnet <- unlist(unname(auc_cpm_mrnet));cpm_clr <- unlist(unname(auc_cpm_clr))

cpm_df <- cbind(cpm_pcc,cpm_scc,cpm_kcc,cpm_gcc,cpm_bicor,
                cpm_aa,cpm_ma,cpm_mrnet,cpm_clr)
cpm_median_method <- apply(cpm_df,2,median)

##RPKM
rpkm_pcc <- unlist(unname(auc_rpkm_pcc));rpkm_scc <- unlist(unname(auc_rpkm_scc))
rpkm_kcc <- unlist(unname(auc_rpkm_kcc));rpkm_gcc <- unlist(unname(auc_rpkm_gcc))
rpkm_bicor <- unlist(unname(auc_rpkm_bicor))
rpkm_aa <- unlist(unname(auc_rpkm_aa));rpkm_ma <- unlist(unname(auc_rpkm_ma))
rpkm_mrnet <- unlist(unname(auc_rpkm_mrnet));rpkm_clr <- unlist(unname(auc_rpkm_clr))

rpkm_df <- cbind(rpkm_pcc,rpkm_scc,rpkm_kcc,rpkm_gcc,rpkm_bicor,
                rpkm_aa,rpkm_ma,rpkm_mrnet,rpkm_clr)
rpkm_median_method <- apply(rpkm_df,2,median)

median_threeMethod <- rbind(vst_median_method,cpm_median_method,rpkm_median_method)
colnames(median_threeMethod) <- c("pcc","scc","kcc","gcc","bicor",
                                     "aa","ma","mrnet","clr")
write.table(median_threeMethod,file="AUCmedian_threePPPTYNormAllConst.txt",sep="\t",quote=F)

rm(list=ls(pattern="^auc_cpm_"));rm(list=ls(pattern="^auc_vst_"));rm(list=ls(pattern="^auc_rpkm_"))
rm(list=ls(pattern="^auc_rc_"));rm(list=ls(pattern="^cpm_"));rm(list=ls(pattern="^rpkm_"))
rm(list=ls(pattern="^vst_"))


# GO
## VST
vst_pcc <- GO_vst_pcc[[1]][,1];vst_scc <- GO_vst_scc[[1]][,1]
vst_kcc <- GO_vst_kcc[[1]][,1];vst_gcc <- GO_vst_gcc[[1]][,1]
vst_bicor <- GO_vst_bicor[[1]][,1]
vst_aa <- GO_vst_aa[[1]][,1];vst_ma <- GO_vst_ma[[1]][,1]
vst_mrnet <- GO_vst_mrnet[[1]][,1];vst_clr <- GO_vst_clr[[1]][,1]

total_vst <- cbind(vst_pcc,vst_scc,vst_kcc,vst_gcc,vst_bicor,vst_aa,vst_ma,vst_mrnet,vst_clr)
vst_go_median <- apply(total_vst,2,median)


# CPM
cpm_pcc <- GO_cpm_pcc[[1]][,1];cpm_scc <- GO_cpm_scc[[1]][,1]
cpm_kcc <- GO_cpm_kcc[[1]][,1];cpm_gcc <- GO_cpm_gcc[[1]][,1]
cpm_bicor <- GO_cpm_bicor[[1]][,1]
cpm_aa <- GO_cpm_aa[[1]][,1];cpm_ma <- GO_cpm_ma[[1]][,1]
cpm_mrnet <- GO_cpm_mrnet[[1]][,1];cpm_clr <- GO_cpm_clr[[1]][,1]

total_cpm <- cbind(cpm_pcc,cpm_scc,cpm_kcc,cpm_gcc,cpm_bicor,cpm_aa,cpm_ma,cpm_mrnet,cpm_clr)
cpm_go_median <- apply(total_cpm,2,median)

#RPKM
rpkm_pcc <- GO_rpkm_pcc[[1]][,1];rpkm_scc <- GO_rpkm_scc[[1]][,1]
rpkm_kcc <- GO_rpkm_kcc[[1]][,1];rpkm_gcc <- GO_rpkm_gcc[[1]][,1]
rpkm_bicor <- GO_rpkm_bicor[[1]][,1]
rpkm_aa <- GO_rpkm_aa[[1]][,1];rpkm_ma <- GO_rpkm_ma[[1]][,1]
rpkm_mrnet <- GO_rpkm_mrnet[[1]][,1];rpkm_clr <- GO_rpkm_clr[[1]][,1]

total_rpkm <- cbind(rpkm_pcc,rpkm_scc,rpkm_kcc,rpkm_gcc,rpkm_bicor,rpkm_aa,rpkm_ma,rpkm_mrnet,rpkm_clr)
rpkm_go_median <- apply(total_rpkm,2,median)

median_go_threeMethod <- rbind(vst_go_median,cpm_go_median,rpkm_go_median)
colnames(median_go_threeMethod) <- c("pcc","scc","kcc","gcc","bicor",
                                     "aa","ma","mrnet","clr")
write.table(median_go_threeMethod,file="AUCGOmedian_threeNormAllConst.txt",sep="\t",quote=F)

rm(list=ls(pattern="^cpm_"));rm(list=ls(pattern="^vst"));rm(list=ls(pattern="^rpkm"))
rm(list=ls(pattern="^GO_vst"));rm(list=ls(pattern="^GO_rpkm_"));rm(list=ls(pattern="^GO_cpm_"))
rm(list=ls(pattern="^GO_rc"))


#############################
# for different size comparision
#PPPTY
##PCC
p_12 <- unlist(unname(auc_12_pcc));p_36 <- unlist(unname(auc_36_pcc))
p_65 <- unlist(unname(auc_65_pcc));p_108 <- unlist(unname(auc_108_pcc))
p_270 <- unlist(unname(auc_270_pcc));p_404 <- unlist(unname(auc_404_pcc))
p_1266 <- unlist(unname(auc_cpm_pcc));
p_six <- unlist(unname(auc_aggSixAgg_pcc));
p_fif <- unlist(unname(auc_agg_pcc_noRank));p_pr <- unlist(unname(auc_pr_norm_pcc))

PPPTY_PCC <- cbind(p_12,p_36,p_65,p_108,p_270,p_404,p_1266,p_six,p_fif,p_pr)
pcc_pppty_median <- apply(PPPTY_PCC,2,median)

#MRNET
p_12 <- unlist(unname(auc_12_mrnet));p_36 <- unlist(unname(auc_36_mrnet))
p_65 <- unlist(unname(auc_65_mrnet));p_108 <- unlist(unname(auc_108_mrnet))
p_270 <- unlist(unname(auc_270_mrnet));p_404 <- unlist(unname(auc_404_mrnet))
p_1266 <- unlist(unname(auc_cpm_mrnet));
p_six <- unlist(unname(auc_aggSixAgg_mrnet));
p_fif <- unlist(unname(auc_agg_mrnet_noRank));p_pr <- unlist(unname(auc_pr_norm_mrnet))

PPPTY_mrnet <- cbind(p_12,p_36,p_65,p_108,p_270,p_404,p_1266,p_six,p_fif,p_pr)
mrnet_pppty_median <- apply(PPPTY_mrnet,2,median)

#clr
p_12 <- unlist(unname(auc_12_clr));p_36 <- unlist(unname(auc_36_clr))
p_65 <- unlist(unname(auc_65_clr));p_108 <- unlist(unname(auc_108_clr))
p_270 <- unlist(unname(auc_270_clr));p_404 <- unlist(unname(auc_404_clr))
p_1266 <- unlist(unname(auc_cpm_clr));
p_six <- unlist(unname(auc_aggSixAgg_clr));
p_fif <- unlist(unname(auc_agg_clr_noRank));p_pr <- unlist(unname(auc_pr_norm_clr))

PPPTY_clr <- cbind(p_12,p_36,p_65,p_108,p_270,p_404,p_1266,p_six,p_fif,p_pr)
clr_pppty_median <- apply(PPPTY_clr,2,median)

median_pppty <- rbind(pcc_pppty_median,mrnet_pppty_median,clr_pppty_median)
write.table(median_pppty,file="AUCPPPTYmedian_SizeNtwk.txt",sep="\t",quote=F)
rm(list=ls(pattern="^p_"));rm(list=ls(pattern="^auc_"));rm(list=ls(pattern="^PPPTY"))

#############################
# GO

#PCC
p_12 <- unname(GO_12_pcc[[1]][,1]);p_36 <- unname(GO_36_pcc[[1]][,1]);p_65 <- unname(GO_65_pcc[[1]][,1])
p_108 <- unname(GO_108_pcc[[1]][,1]);p_270 <- unname(GO_270_pcc[[1]][,1]);p_404 <- unname(GO_404_pcc[[1]][,1])
p_1266 <- unname(GO_cpm_pcc[[1]][,1]);p_six <- unname(GO_six_pcc[[1]][,1]);
p_fif <- unname(GO_agg_pcc_noRank[[1]][,1]);p_pr <- unname(GO_pr_pcc[[1]][,1])

total_pcc <- cbind(p_12,p_36,p_65,p_108,p_270,p_404,p_1266,p_six,p_fif)
pcc_go_median <- apply(total_pcc,2,median)
pcc_go_median <- c(pcc_go_median,median(p_pr))

#MRNET
p_12 <- GO_12_mrnet[[1]][,1];p_36 <- GO_36_mrnet[[1]][,1];p_65 <- GO_65_mrnet[[1]][,1]
p_108 <- GO_108_mrnet[[1]][,1];p_270 <- GO_270_mrnet[[1]][,1];p_404 <- GO_404_mrnet[[1]][,1]
p_1266 <- GO_cpm_mrnet[[1]][,1];p_six <- GO_six_mrnet[[1]][,1];
p_fif <- GO_agg_mrnet_noRank[[1]][,1];p_pr <- GO_pr_mrnet[[1]][,1]

total_mrnet <- cbind(p_12,p_36,p_65,p_108,p_270,p_404,p_1266,p_six,p_fif)
mrnet_go_median <- apply(total_mrnet,2,median)
mrnet_go_median <- c(mrnet_go_median,median(p_pr))

#CLR
p_12 <- GO_12_clr[[1]][,1];p_36 <- GO_36_clr[[1]][,1];p_65 <- GO_65_clr[[1]][,1]
p_108 <- GO_108_clr[[1]][,1];p_270 <- GO_270_clr[[1]][,1];p_404 <- GO_404_clr[[1]][,1]
p_1266 <- GO_cpm_clr[[1]][,1];p_six <- GO_six_clr[[1]][,1];
p_fif <- GO_agg_clr_noRank[[1]][,1];p_pr <- GO_pr_clr[[1]][,1]

total_clr <- cbind(p_12,p_36,p_65,p_108,p_270,p_404,p_1266,p_six,p_fif)
clr_go_median <- apply(total_clr,2,median)
clr_go_median <- c(clr_go_median,median(p_pr))


median_go <- rbind(pcc_go_median,mrnet_go_median,clr_go_median)
write.table(median_go,file="AUCGOmedian_SizeNtwk.txt",sep="\t",quote=F)
rm(list=ls(pattern="^p_"));rm(list=ls(pattern="^auc_"));rm(list=ls(pattern="^PPPTY"))
