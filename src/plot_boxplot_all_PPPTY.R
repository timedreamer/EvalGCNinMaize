# for maize I decided to ignore all_fifteeen but use fif_rank instead

#PPPTY_six_individual
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP/aggregate_ntwk/AUC_cpm12.RData")
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP/aggregate_ntwk/AUC_cpm36.RData")
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP/aggregate_ntwk/AUC_cpm65.RData")
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP/aggregate_ntwk/AUC_cpm108.RData")
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP/aggregate_ntwk/AUC_cpm270.RData")
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP/aggregate_ntwk/AUC_cpm404.RData")
#PPPTY_1266
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP/total_ntwk/AUC_cpm.RData")
rm(auc_cpm_aa,auc_cpm_ma,auc_cpm_bicor,auc_cpm_scc,auc_cpm_kcc,auc_cpm_gcc)
#PPPTY_sixAggregate
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP/aggregate_ntwk/agg_cpm_pcc_PPPTY.RData")
auc_aggSixAgg_pcc <- auc_agg_test 
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP/aggregate_ntwk/agg_cpm_mrnet_PPPTY.RData")
auc_aggSixAgg_mrnet <- auc_agg_test
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP/aggregate_ntwk/agg_cpm_clr_PPPTY.RData")
auc_aggSixAgg_clr <- auc_agg_test;rm(auc_agg_test)
#PPPTY_15Aggregate
# load("~/Co-expression/aggregate_ntwk/withExp14_total15Exps/agg_cpm_pccAll_PPPTY_noRank.RData")
# load("~/Co-expression/aggregate_ntwk/withExp14_total15Exps/agg_cpm_mrnetAll_PPPTY_noRank.RData")
# load("~/Co-expression/aggregate_ntwk/withExp14_total15Exps/agg_cpm_clrAll_PPPTY_noRank.RData")
#PPPTY_protein
load("~/Co-expression/Protein/protein_auc_PPPTY.RData")
#PPPTY_fifteen_ranked
load("~/Co-expression/aggregate_ntwk/MaizeRank/agg_rank_maize_PPPTYResult.RData")



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
# p_fif <- unlist(unname(auc_agg_pcc_noRank));names(p_fif) <- rep("fifteen",1721)
p_fif_rank <- unlist(unname(auc_aggRank_pcc_maize));names(p_fif_rank) <- rep("fif_rank",1721)
p_pr <- unlist(unname(auc_pr_norm_pcc));names(p_pr) <- rep("protein",2402)


total_pcc <- c(p_12,p_36,p_65,p_108,p_270,p_404,p_1266,p_six,p_fif_rank,p_pr)
total_name_pcc <- c(names(p_12),names(p_36),names(p_65),names(p_108),names(p_270),
                    names(p_404),names(p_1266),names(p_six),names(p_fif_rank),names(p_pr))
total_name_pcc <- factor(total_name_pcc,level = c("12","36","65","108","270","404","1266",
                                                  "six","fif_rank","protein"))
method_pcc <- rep("pcc",17891) # 1721*10 + 2402
mean_pcc <- c(mean(p_12),mean(p_36),mean(p_65),mean(p_108),mean(p_270),mean(p_404),
              mean(p_1266),mean(p_six),mean(p_fif_rank),mean(p_pr))

#MRNET
p_12 <- unlist(unname(auc_12_mrnet));names(p_12) <- rep("12",1721)
p_36 <- unlist(unname(auc_36_mrnet));names(p_36) <- rep("36",1721)
p_65 <- unlist(unname(auc_65_mrnet));names(p_65) <- rep("65",1721)
p_108 <- unlist(unname(auc_108_mrnet));names(p_108) <- rep("108",1721)
p_270 <- unlist(unname(auc_270_mrnet));names(p_270) <- rep("270",1721)
p_404 <- unlist(unname(auc_404_mrnet));names(p_404) <- rep("404",1721)
p_1266 <- unlist(unname(auc_cpm_mrnet));names(p_1266) <- rep("1266",1721)
p_six <- unlist(unname(auc_aggSixAgg_mrnet));names(p_six) <- rep("six",1721)
# p_fif <- unlist(unname(auc_agg_mrnet_noRank));names(p_fif) <- rep("fifteen",1721)
p_fif_rank <- unlist(unname(auc_aggRank_mrnet_maize));names(p_fif_rank) <- rep("fif_rank",1721)
p_pr <- unlist(unname(auc_pr_norm_mrnet));names(p_pr) <- rep("protein",2402)


total_mrnet <- c(p_12,p_36,p_65,p_108,p_270,p_404,p_1266,p_six,p_fif_rank,p_pr)
total_name_mrnet <- c(names(p_12),names(p_36),names(p_65),names(p_108),names(p_270),
                      names(p_404),names(p_1266),names(p_six),names(p_fif_rank),names(p_pr))
total_name_mrnet <- factor(total_name_mrnet,level = c("12","36","65","108","270","404","1266",
                                                      "six","fif_rank","protein"))
method_mrnet <- rep("mrnet",17891)
mean_mrnet <- c(mean(p_12),mean(p_36),mean(p_65),mean(p_108),mean(p_270),mean(p_404),
                mean(p_1266),mean(p_six),mean(p_fif_rank),mean(p_pr))

#CLR
p_12 <- unlist(unname(auc_12_clr));names(p_12) <- rep("12",1721)
p_36 <- unlist(unname(auc_36_clr));names(p_36) <- rep("36",1721)
p_65 <- unlist(unname(auc_65_clr));names(p_65) <- rep("65",1721)
p_108 <- unlist(unname(auc_108_clr));names(p_108) <- rep("108",1721)
p_270 <- unlist(unname(auc_270_clr));names(p_270) <- rep("270",1721)
p_404 <- unlist(unname(auc_404_clr));names(p_404) <- rep("404",1721)
p_1266 <- unlist(unname(auc_cpm_clr));names(p_1266) <- rep("1266",1721)
p_six <- unlist(unname(auc_aggSixAgg_clr));names(p_six) <- rep("six",1721)
# p_fif <- unlist(unname(auc_agg_clr_noRank));names(p_fif) <- rep("fifteen",1721)
p_fif_rank <- unlist(unname(auc_aggRank_clr_maize));names(p_fif_rank) <- rep("fif_rank",1721)
p_pr <- unlist(unname(auc_pr_norm_clr));names(p_pr) <- rep("protein",2402)


total_clr <- c(p_12,p_36,p_65,p_108,p_270,p_404,p_1266,p_six,p_fif_rank,p_pr)
total_name_clr <- c(names(p_12),names(p_36),names(p_65),names(p_108),names(p_270),
                    names(p_404),names(p_1266),names(p_six),names(p_fif_rank),names(p_pr))
total_name_clr <- factor(total_name_clr,level = c("12","36","65","108","270","404","1266",
                                                  "six","fif_rank","protein"))
method_clr <- rep("clr",17891)
mean_clr <- c(mean(p_12),mean(p_36),mean(p_65),mean(p_108),mean(p_270),mean(p_404),
              mean(p_1266),mean(p_six),mean(p_fif_rank),mean(p_pr))

# one-factor boxplot. export svg as 800*500
par(mar=c(9,5,4,1),mfrow=c(1,3))
boxplot(total_pcc~total_name_pcc,main="PPPTY_PCC",las=2,ylim=c(0.4,0.75),
        color=c("white"),
        outcol="grey",whiskcol="grey",staplecol="grey",cex.axis=2,
        names=c("12","36","65","108","270","404","1266","six","fif_rank","protein"));box(lwd=2)
points(mean_pcc,col="black",pch=8)

boxplot(total_mrnet~total_name_mrnet,main="PPPTY_MRNET",las=2,ylim=c(0.4,0.75),
        outcol="grey",whiskcol="grey",staplecol="grey",cex.axis=2,
        names=c("12","36","65","108","270","404","1266","six","fif_rank","protein"));box(lwd=2)
points(mean_mrnet,col="black",pch=8)


boxplot(total_clr~total_name_clr,main="PPPTY_CLR",las=2,ylim=c(0.4,0.75),
        outcol="grey",whiskcol="grey",staplecol="grey",cex.axis=2,
        names=c("12","36","65","108","270","404","1266","six","fif_rank","protein"));box(lwd=2)
points(mean_clr,col="black",pch=8)


# Plot everything together in boxplot .export as svg 600*600
total_value <- c(total_pcc,total_mrnet,total_clr)
total_netname <- rep(c(rep("12",1721),rep("36",1721),rep("65",1721),rep("108",1721),rep("270",1721),
                       rep("404",1721),rep("1266",1721),rep("six",1721),
                       rep("fif_rank",1721),rep("protein",2402)),3)
total_constructMethod <- c(method_pcc,method_mrnet,method_clr)


total_netname <- factor(total_netname,
                        c("12","36","65","108","270","404","1266","six","fif_rank","protein"))
total_constructMethod <- factor(total_constructMethod,c("pcc","mrnet","clr"))


par(mfrow=c(1,1),mar=c(7,4,4,2))
boxplot(total_value~total_constructMethod*total_netname,las=2,yaxt="n",ylim=c(0.4,0.75),
        col=rep(c(rep("white",3),rep("grey",3)),5),
        outcol="grey",whiskcol="grey",staplecol="grey",main="allNtwk_PPPTY")
box(lwd=2);axis(2,cex.axis=1.5,las=1)

boxplot(total_value~total_constructMethod*total_netname,las=2,yaxt="n",ylim=c(0.4,0.75),
        col=c("white","grey","grey34"),
        outcol="grey",whiskcol="grey",staplecol="grey",main="allNtwk_PPPTY")
box(lwd=2);axis(2,cex.axis=1.5,las=1)


##########################################################################
# Pairwise wilconxon test
pairwise.wilcox.test(total_pcc,total_name_pcc,p.adjust.method = "b",correct=F)
pairwise.wilcox.test(total_mrnet,total_name_mrnet,p.adjust.method = "b",correct=F)
pairwise.wilcox.test(total_clr,total_name_clr,p.adjust.method = "b",correct=F)


### how about each size
#1266
pcc_1266 <- unlist(unname(auc_cpm_pcc));mrnet_1266 <- unlist(unname(auc_cpm_mrnet))
clr_1266 <- unlist(unname(auc_cpm_clr))
total_1266 <- c(pcc_1266,mrnet_1266,clr_1266)
totalName_1266 <- c(rep('pcc',1721),rep('mrnet',1721),rep('clr',1721))
pairwise.wilcox.test(total_1266,totalName_1266,p.adjust.method = 'b',correct=F)

#six
p_six <- unlist(unname(auc_aggSixAgg_pcc));m_six <- unlist(unname(auc_aggSixAgg_mrnet))
c_six <- unlist(unname(auc_aggSixAgg_clr))
total_six <- c(p_six,m_six,c_six)
totalName_six <- c(rep('pcc',1721),rep('mrnet',1721),rep('clr',1721))
pairwise.wilcox.test(total_six,totalName_six,p.adjust.method = 'b',correct=F)

# fifteen
p_fif <- unlist(unname(auc_agg_pcc_noRank));m_fif <- unlist(unname(auc_agg_mrnet_noRank))
c_fif <- unlist(unname(auc_agg_clr_noRank))
total_fif <- c(p_fif,m_fif,c_fif)
totalName_fif <- c(rep('pcc',1721),rep('mrnet',1721),rep('clr',1721))
pairwise.wilcox.test(total_fif,totalName_fif,p.adjust.method = 'b',correct=F)

# fif_rank
p_fif_rank <- unlist(unname(auc_aggRank_pcc_maize));m_fif_rank <- unlist(unname(auc_aggRank_mrnet_maize))
c_fif_rank <- unlist(unname(auc_aggRank_clr_maize))
total_fif_rank <- c(p_fif_rank,m_fif_rank,c_fif_rank)
totalName_fif_rank <- c(rep('pcc',1721),rep('mrnet',1721),rep('clr',1721))
pairwise.wilcox.test(total_fif_rank,totalName_fif_rank,p.adjust.method = 'b',correct=F)