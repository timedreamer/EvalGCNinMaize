
# load six experiments PPPTY and GO evaluation data.
#PPPTY_six_individual
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP/aggregate_ntwk/AUC_cpm12.RData")
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP/aggregate_ntwk/AUC_cpm36.RData")
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP/aggregate_ntwk/AUC_cpm65.RData")
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP/aggregate_ntwk/AUC_cpm108.RData")
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP/aggregate_ntwk/AUC_cpm270.RData")
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP/aggregate_ntwk/AUC_cpm404.RData")
#GO_six_individual
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/GO_AUROC/aggregate_ntwk/GO_AUROC_SixAggregate_ntwk.RData")


#PPPTY_1266
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP/total_ntwk/AUC_cpm.RData")
rm(auc_cpm_aa,auc_cpm_ma,auc_cpm_bicor,auc_cpm_scc,auc_cpm_kcc,auc_cpm_gcc)
#GO_1266
load("~/Co-expression/Normalization result/Before_Oct2016/AUC/GO_cpm.RData")
rm(GO_cpm_scc,GO_cpm_kcc,GO_cpm_gcc,GO_cpm_aa,GO_cpm_ma,GO_cpm_bicor)

#PPPTY_sixAggregate
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP/aggregate_ntwk/agg_cpm_pcc_PPPTY.RData")
auc_aggSixAgg_pcc <- auc_agg_test 
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP/aggregate_ntwk/agg_cpm_mrnet_PPPTY.RData")
auc_aggSixAgg_mrnet <- auc_agg_test
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP/aggregate_ntwk/agg_cpm_clr_PPPTY.RData")
auc_aggSixAgg_clr <- auc_agg_test;rm(auc_agg_test)
#GO_sixAggregate
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/GO_AUROC/aggregate_ntwk/agg_cpm_pcc_GO.RData")
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/GO_AUROC/aggregate_ntwk/agg_cpm_clr_GO.RData")
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/GO_AUROC/aggregate_ntwk/agg_cpm_mrnet_GO.RData")

#PPPTY_15Aggregate
load("~/Co-expression/aggregate_ntwk/withExp14_total15Exps/agg_cpm_pccAll_PPPTY_noRank.RData")
load("~/Co-expression/aggregate_ntwk/withExp14_total15Exps/agg_cpm_mrnetAll_PPPTY_noRank.RData")
load("~/Co-expression/aggregate_ntwk/withExp14_total15Exps/agg_cpm_clrAll_PPPTY_noRank.RData")
#GO_15Aggregate
load("~/Co-expression/aggregate_ntwk/withExp14_total15Exps/agg_cpmClrAll_GO_noRank.RData")
load("~/Co-expression/aggregate_ntwk/withExp14_total15Exps/agg_cpmPccAll_GO_noRank.RData")
load("~/Co-expression/aggregate_ntwk/withExp14_total15Exps/agg_cpmMrnetAll_GO_noRank.RData")


#PPPTY_protein
load("~/Co-expression/Protein/protein_auc_PPPTY.RData")
#GO_protein
load("~/Co-expression/Protein/GO_AUROC_protein.RData")


##For PPPTY
#PCC
p_12 <- unlist(unname(auc_12_pcc));names(p_12) <- rep("12",1721)
p_36 <- unlist(unname(auc_36_pcc));names(p_36) <- rep("36",1721)
p_65 <- unlist(unname(auc_65_pcc));names(p_65) <- rep("65",1721)
p_108 <- unlist(unname(auc_108_pcc));names(p_108) <- rep("108",1721)
p_270 <- unlist(unname(auc_270_pcc));names(p_270) <- rep("270",1721)
p_404 <- unlist(unname(auc_404_pcc));names(p_404) <- rep("404",1721)
p_1266 <- unlist(unname(auc_cpm_pcc));names(p_1266) <- rep("1266",1721)
p_six <- unlist(unname(auc_aggSixAgg_pcc));names(p_six) <- rep("six",1721)
p_fif <- unlist(unname(auc_agg_pcc_noRank));names(p_fif) <- rep("all_15",1721)
p_pr <- unlist(unname(auc_pr_norm_pcc));names(p_pr) <- rep("protein",2402)


total_pcc <- c(p_12,p_36,p_65,p_108,p_270,p_404,p_1266,p_six,p_fif,p_pr)
total_name_pcc <- c(names(p_12),names(p_36),names(p_65),names(p_108),names(p_270),
                    names(p_404),names(p_1266),names(p_six),names(p_fif),names(p_pr))
method_pcc <- rep("pcc",17891) # 1721*9+2402

#MRNET
p_12 <- unlist(unname(auc_12_mrnet));names(p_12) <- rep("12",1721)
p_36 <- unlist(unname(auc_36_mrnet));names(p_36) <- rep("36",1721)
p_65 <- unlist(unname(auc_65_mrnet));names(p_65) <- rep("65",1721)
p_108 <- unlist(unname(auc_108_mrnet));names(p_108) <- rep("108",1721)
p_270 <- unlist(unname(auc_270_mrnet));names(p_270) <- rep("270",1721)
p_404 <- unlist(unname(auc_404_mrnet));names(p_404) <- rep("404",1721)
p_1266 <- unlist(unname(auc_cpm_mrnet));names(p_1266) <- rep("1266",1721)
p_six <- unlist(unname(auc_aggSixAgg_mrnet));names(p_six) <- rep("six",1721)
p_fif <- unlist(unname(auc_agg_mrnet_noRank));names(p_fif) <- rep("all_15",1721)
p_pr <- unlist(unname(auc_pr_norm_mrnet));names(p_pr) <- rep("protein",2402)


total_mrnet <- c(p_12,p_36,p_65,p_108,p_270,p_404,p_1266,p_six,p_fif,p_pr)
total_name_mrnet <- c(names(p_12),names(p_36),names(p_65),names(p_108),names(p_270),
                    names(p_404),names(p_1266),names(p_six),names(p_fif),names(p_pr))
method_mrnet <- rep("mrnet",17891)

#CLR
p_12 <- unlist(unname(auc_12_clr));names(p_12) <- rep("12",1721)
p_36 <- unlist(unname(auc_36_clr));names(p_36) <- rep("36",1721)
p_65 <- unlist(unname(auc_65_clr));names(p_65) <- rep("65",1721)
p_108 <- unlist(unname(auc_108_clr));names(p_108) <- rep("108",1721)
p_270 <- unlist(unname(auc_270_clr));names(p_270) <- rep("270",1721)
p_404 <- unlist(unname(auc_404_clr));names(p_404) <- rep("404",1721)
p_1266 <- unlist(unname(auc_cpm_clr));names(p_1266) <- rep("1266",1721)
p_six <- unlist(unname(auc_aggSixAgg_clr));names(p_six) <- rep("six",1721)
p_fif <- unlist(unname(auc_agg_clr_noRank));names(p_fif) <- rep("all_15",1721)
p_pr <- unlist(unname(auc_pr_norm_clr));names(p_pr) <- rep("protein",2402)


total_clr <- c(p_12,p_36,p_65,p_108,p_270,p_404,p_1266,p_six,p_fif,p_pr)
total_name_clr <- c(names(p_12),names(p_36),names(p_65),names(p_108),names(p_270),
                      names(p_404),names(p_1266),names(p_six),names(p_fif),names(p_pr))
method_clr <- rep("clr",17891)

# Add together
total_value <- c(total_pcc,total_mrnet,total_clr)
total_netname <- c(total_name_pcc,total_name_mrnet,total_name_clr)
total_constructMethod <- c(method_pcc,method_mrnet,method_clr)


v2 <- data.frame(total_value,total_netname,total_constructMethod)

# Two-way anova and pairwise comparison
v2$total_constructMethod <- factor(v2$total_constructMethod,levels = c("pcc","mrnet","clr"))
v2$total_netname <- factor(v2$total_netname,
                             level = c("12","36","65","108","270","404","1266",
                                       "six","all_15","protein"))
# plot is the same as plot singe var.
par(mfrow=c(1,2)) # export pdf 16*9
plot(total_value ~ total_netname*total_constructMethod,data=v2,ylab = "AUC",xlab="var")

anova(lm(total_value ~ total_netname*total_constructMethod,data=v2)) # show anova table directly


pairwise.t.test(total_value, total_netname, p.adjust="bonferroni")
pairwise.t.test(total_value, total_constructMethod, p.adjust="bonferroni")

pairwise.wilcox.test(total_value, total_netname, p.adjust="bonferroni")
pairwise.wilcox.test(total_value, total_constructMethod, p.adjust="bonferroni")

##########################################################################################################
##########################################################################################################
##########################################################################################################

# GO

#PCC
p_12 <- GO_12_pcc[[1]][,1];p_36 <- GO_36_pcc[[1]][,1];p_65 <- GO_65_pcc[[1]][,1]
p_108 <- GO_108_pcc[[1]][,1];p_270 <- GO_270_pcc[[1]][,1];p_404 <- GO_404_pcc[[1]][,1]
p_1266 <- GO_cpm_pcc[[1]][,1];p_six <- GO_pcc_agg_ntwk[[1]][,1];
p_fif <- GO_agg_pcc_noRank[[1]][,1];p_pr <- GO_pr_pcc[[1]][,1]

total_pcc <- c(p_12,p_36,p_65,p_108,p_270,p_404,p_1266,p_six,p_fif,p_pr)
total_name_pcc <- c(rep("12",277),rep("36",277),rep("65",277),rep("108",277),rep("270",277),
                   rep("404",277),rep("1266",277),
                   rep("six",277),rep("all_15",277),rep("protein",307))
method_pcc <- rep("pcc",2800) #277*9+307

#MRNET
p_12 <- GO_12_mrnet[[1]][,1];p_36 <- GO_36_mrnet[[1]][,1];p_65 <- GO_65_mrnet[[1]][,1]
p_108 <- GO_108_mrnet[[1]][,1];p_270 <- GO_270_mrnet[[1]][,1];p_404 <- GO_404_mrnet[[1]][,1]
p_1266 <- GO_cpm_mrnet[[1]][,1];p_six <- GO_mrnet_agg_ntwk[[1]][,1];
p_fif <- GO_agg_mrnet_noRank[[1]][,1];p_pr <- GO_pr_mrnet[[1]][,1]

total_mrnet <- c(p_12,p_36,p_65,p_108,p_270,p_404,p_1266,p_six,p_fif,p_pr)
total_name_mrnet <- c(rep("12",277),rep("36",277),rep("65",277),rep("108",277),rep("270",277),
                    rep("404",277),rep("1266",277),
                    rep("six",277),rep("all_15",277),rep("protein",307))
method_mrnet <- rep("mrnet",2800) #277*9+307

#CLR
p_12 <- GO_12_clr[[1]][,1];p_36 <- GO_36_clr[[1]][,1];p_65 <- GO_65_clr[[1]][,1]
p_108 <- GO_108_clr[[1]][,1];p_270 <- GO_270_clr[[1]][,1];p_404 <- GO_404_clr[[1]][,1]
p_1266 <- GO_cpm_clr[[1]][,1];p_six <- GO_clr_agg_ntwk[[1]][,1];
p_fif <- GO_agg_clr_noRank[[1]][,1];p_pr <- GO_pr_clr[[1]][,1]

total_clr <- c(p_12,p_36,p_65,p_108,p_270,p_404,p_1266,p_six,p_fif,p_pr)
total_name_clr <- c(rep("12",277),rep("36",277),rep("65",277),rep("108",277),rep("270",277),
                    rep("404",277),rep("1266",277),rep("six",277),
                    rep("all_15",277),rep("protein",307))
method_clr <- rep("clr",2800) #277*9+307

# Add together
total_value <- c(total_pcc,total_mrnet,total_clr)
total_netname <- c(total_name_pcc,total_name_mrnet,total_name_clr)
total_constructMethod <- c(method_pcc,method_mrnet,method_clr)

v2 <- data.frame(total_value,total_netname,total_constructMethod)

# Two-way anova and pairwise comparison
v2$total_constructMethod <- factor(v2$total_constructMethod,levels = c("pcc","mrnet","clr"))
v2$total_netname <- factor(v2$total_netname,
                           level = c("12","36","65","108","270","404","1266",
                                     "six","all_15","protein"))
# plot is the same as plot singe var.
par(mfrow=c(1,2)) # export pdf 16*9
plot(total_value ~ total_netname*total_constructMethod,data=v2,ylab = "AUC",xlab="var")

anova(lm(total_value ~ total_netname*total_constructMethod,data=v2)) # show anova table directly


pairwise.t.test(total_value, total_netname, p.adjust="bonferroni")
pairwise.t.test(total_value, total_constructMethod, p.adjust="bonferroni")

pairwise.wilcox.test(total_value, total_netname, p.adjust="bonferroni")
pairwise.wilcox.test(total_value, total_constructMethod, p.adjust="bonferroni")


