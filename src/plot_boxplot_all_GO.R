

# load data
#GO_six_individual
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/GO_AUROC/aggregate_ntwk/GO_AUROC_SixAggregate_ntwk.RData")
#GO_1266
load("~/Co-expression/Normalization result/Before_Oct2016/AUC/GO_cpm.RData")
rm(GO_cpm_scc,GO_cpm_kcc,GO_cpm_gcc,GO_cpm_aa,GO_cpm_ma,GO_cpm_bicor)

#GO_sixAggregate
load("D:\\Users\\jhuang\\Documents\\Co-expression\\aggregate_ntwk\\six\\agg_cpm_pcc_GO_noRank.RData")
GO_six_pcc <- GO_agg_pcc_noRank;rm(GO_agg_pcc_noRank)
load("D:\\Users\\jhuang\\Documents\\Co-expression\\aggregate_ntwk\\six\\agg_cpm_mrnet_GO_noRank.RData")
GO_six_mrnet <- GO_agg_mrnet_noRank ;rm(GO_agg_mrnet_noRank)
load("D:\\Users\\jhuang\\Documents\\Co-expression\\aggregate_ntwk\\six\\agg_cpm_clr_GO_noRank.RData")
GO_six_clr <- GO_agg_clr_noRank ;rm(GO_agg_clr_noRank)

#GO_15Aggregate
load("~/Co-expression/aggregate_ntwk/withExp14_total15Exps/agg_cpmClrAll_GO_noRank.RData")
load("~/Co-expression/aggregate_ntwk/withExp14_total15Exps/agg_cpmPccAll_GO_noRank.RData")
load("~/Co-expression/aggregate_ntwk/withExp14_total15Exps/agg_cpmMrnetAll_GO_noRank.RData")

# GO_15_ranked
load("~/Co-expression/aggregate_ntwk/MaizeRank/GO_evaluation/GO_agg_rankPcc_Maize.RData")
load("~/Co-expression/aggregate_ntwk/MaizeRank/GO_evaluation/GO_agg_rankMrnet_Maize.RData")
load("~/Co-expression/aggregate_ntwk/MaizeRank/GO_evaluation/GO_agg_rankClr_Maize.RData")

#GO_protein
load("~/Co-expression/Protein/GO_AUROC_protein.RData")



#PCC
p_12 <- GO_12_pcc[[1]][,1];p_36 <- GO_36_pcc[[1]][,1];p_65 <- GO_65_pcc[[1]][,1]
p_108 <- GO_108_pcc[[1]][,1];p_270 <- GO_270_pcc[[1]][,1];p_404 <- GO_404_pcc[[1]][,1]
p_1266 <- GO_cpm_pcc[[1]][,1];p_six <- GO_six_pcc[[1]][,1];
p_fif <- GO_agg_pcc_noRank[[1]][,1];p_fif_rank <- GO_agg_rankPccMaize[[1]][,1]
p_pr <- GO_pr_pcc[[1]][,1]

total_pcc <- c(p_12,p_36,p_65,p_108,p_270,p_404,p_1266,p_six,p_fif,p_fif_rank,p_pr)
total_name_pcc <- c(rep("12",277),rep("36",277),rep("65",277),rep("108",277),rep("270",277),
                    rep("404",277),rep("1266",277),
                    rep("six",277),rep("fifteen",277),rep("fif_rank",277),rep("protein",307))
total_name_pcc <- factor(total_name_pcc,level = c("12","36","65","108","270","404","1266",
                                                  "six","fifteen","fif_rank","protein"))
method_pcc <- rep("pcc",3077) #277*10+307
mean_pcc <- c(mean(p_12),mean(p_36),mean(p_65),mean(p_108),mean(p_270),mean(p_404),
              mean(p_1266),mean(p_six),mean(p_fif),mean(p_fif_rank),mean(p_pr))

#MRNET
p_12 <- GO_12_mrnet[[1]][,1];p_36 <- GO_36_mrnet[[1]][,1];p_65 <- GO_65_mrnet[[1]][,1]
p_108 <- GO_108_mrnet[[1]][,1];p_270 <- GO_270_mrnet[[1]][,1];p_404 <- GO_404_mrnet[[1]][,1]
p_1266 <- GO_cpm_mrnet[[1]][,1];p_six <- GO_six_mrnet[[1]][,1];
p_fif <- GO_agg_mrnet_noRank[[1]][,1];p_fif_rank <- GO_agg_rankMrnetMaize[[1]][,1]
p_pr <- GO_pr_mrnet[[1]][,1]

total_mrnet <- c(p_12,p_36,p_65,p_108,p_270,p_404,p_1266,p_six,p_fif,p_fif_rank,p_pr)
total_name_mrnet <- c(rep("12",277),rep("36",277),rep("65",277),rep("108",277),rep("270",277),
                      rep("404",277),rep("1266",277),
                      rep("six",277),rep("fifteen",277),rep("fif_rank",277),rep("protein",307))
total_name_mrnet <- factor(total_name_mrnet,level = c("12","36","65","108","270","404","1266",
                                                  "six","fifteen","fif_rank","protein"))
method_mrnet <- rep("mrnet",3077) 
mean_mrnet <- c(mean(p_12),mean(p_36),mean(p_65),mean(p_108),mean(p_270),mean(p_404),
                mean(p_1266),mean(p_six),mean(p_fif),mean(p_fif_rank),mean(p_pr))

#CLR
p_12 <- GO_12_clr[[1]][,1];p_36 <- GO_36_clr[[1]][,1];p_65 <- GO_65_clr[[1]][,1]
p_108 <- GO_108_clr[[1]][,1];p_270 <- GO_270_clr[[1]][,1];p_404 <- GO_404_clr[[1]][,1]
p_1266 <- GO_cpm_clr[[1]][,1];p_six <- GO_six_clr[[1]][,1];
p_fif <- GO_agg_clr_noRank[[1]][,1];p_fif_rank <- GO_agg_rankClrMaize[[1]][,1]
p_pr <- GO_pr_clr[[1]][,1]

total_clr <- c(p_12,p_36,p_65,p_108,p_270,p_404,p_1266,p_six,p_fif,p_fif_rank,p_pr)
total_name_clr <- c(rep("12",277),rep("36",277),rep("65",277),rep("108",277),rep("270",277),
                    rep("404",277),rep("1266",277),rep("six",277),
                    rep("fifteen",277),rep("fif_rank",277),rep("protein",307))
total_name_clr <- factor(total_name_clr,level = c("12","36","65","108","270","404","1266",
                                                      "six","fifteen","fif_rank","protein"))
method_clr <- rep("clr",3077) 
mean_clr <- c(mean(p_12),mean(p_36),mean(p_65),mean(p_108),mean(p_270),mean(p_404),
              mean(p_1266),mean(p_six),mean(p_fif),mean(p_fif_rank),mean(p_pr))




# one-factor boxplot. export 8*6
par(mar=c(9,4,4,1),mfrow=c(1,3))
boxplot(total_pcc~total_name_pcc,main="GO_PCC",las=2,
        color=c("white"),
        outcol="grey",whiskcol="grey",staplecol="grey",cex.axis=2,
        names=c("12","36","65","108","270","404","1266","six","fifteen","fif_rank","protein"));box(lwd=2)
points(mean_pcc,col="black",pch=8)

boxplot(total_mrnet~total_name_mrnet,main="GO_MRNET",las=2,
        outcol="grey",whiskcol="grey",staplecol="grey",cex.axis=2,
        names=c("12","36","65","108","270","404","1266","six","fifteen","fif_rank","protein"));box(lwd=2)
points(mean_mrnet,col="black",pch=8)


boxplot(total_clr~total_name_clr,main="GO_CLR",las=2,
        outcol="grey",whiskcol="grey",staplecol="grey",cex.axis=2,
        names=c("12","36","65","108","270","404","1266","six","fifteen","fif_rank","protein"));box(lwd=2)
points(mean_clr,col="black",pch=8)


# Plot everything together in boxplot .export as 6*6
total_value <- c(total_pcc,total_mrnet,total_clr)
total_netname <- rep(c(rep("12",277),rep("36",277),rep("65",277),rep("108",277),rep("270",277),
                   rep("404",277),rep("1266",277),rep("six",277),
                   rep("fifteen",277),rep("fif_rank",277),rep("protein",307)),3)
total_constructMethod <- c(method_pcc,method_mrnet,method_clr)


total_netname <- factor(total_netname,
                           c("12","36","65","108","270","404","1266","six","fifteen","fif_rank","protein"))
total_constructMethod <- factor(total_constructMethod,c("pcc","mrnet","clr"))


par(mfrow=c(1,1),mar=c(7,4,4,2))
boxplot(total_value~total_constructMethod*total_netname,las=2,yaxt="n",
        col=rep(c(rep("white",3),rep("grey",3)),5),
        outcol="grey",whiskcol="grey",staplecol="grey",main="allNtwk")
box(lwd=2);axis(2,cex.axis=1.5,las=1)

boxplot(total_value~total_constructMethod*total_netname,las=2,yaxt="n",
        col=c("white","grey","grey34"),
        outcol="grey",whiskcol="grey",staplecol="grey",main="allNtwk")
box(lwd=2);axis(2,cex.axis=1.5,las=1)



######################################################################
#pairwise wilconxon test
pairwise.wilcox.test(total_pcc,total_name_pcc,p.adjust.method = "b",correct=F)
pairwise.wilcox.test(total_mrnet,total_name_mrnet,p.adjust.method = "b",correct=F)
pairwise.wilcox.test(total_clr,total_name_clr,p.adjust.method = "b",correct=F)


# compare 1266, six and fifteen networks among three methods.
p_pcc_1266 <- GO_cpm_pcc[[1]][,1];p_pcc_six <- GO_pcc_agg_ntwk[[1]][,1];
p_pcc_fif <- GO_agg_pcc_noRank[[1]][,1];p_pcc_fif_rank <- GO_agg_rankPccMaize[[1]][,1]

p_mrnet_1266 <- GO_cpm_mrnet[[1]][,1];p_mrnet_six <- GO_mrnet_agg_ntwk[[1]][,1];
p_mrnet_fif <- GO_agg_mrnet_noRank[[1]][,1];p_mrnet_fif_rank <- GO_agg_rankMrnetMaize[[1]][,1]

p_clr_1266 <- GO_cpm_clr[[1]][,1];p_clr_six <- GO_clr_agg_ntwk[[1]][,1];
p_clr_fif <- GO_agg_clr_noRank[[1]][,1];p_clr_fif_rank <- GO_agg_rankClrMaize[[1]][,1]

all_1266 <- c(p_pcc_1266,p_mrnet_1266,p_clr_1266)
all_1266_name <- c(rep("pcc",277),rep("mrnet",277),rep("clr",277))

all_six <- c(p_pcc_six,p_mrnet_six,p_clr_six)
all_six_name <- c(rep("pcc",277),rep("mrnet",277),rep("clr",277))

all_fif <- c(p_pcc_fif,p_mrnet_fif,p_clr_fif)
all_fif_name <- c(rep("pcc",277),rep("mrnet",277),rep("clr",277))

all_fif_rank <- c(p_pcc_fif_rank,p_mrnet_fif_rank,p_clr_fif_rank)
all_fif_rank_name <-  c(rep("pcc",277),rep("mrnet",277),rep("clr",277))

pairwise.wilcox.test(all_1266,all_1266_name,p.adjust.method = 'b',paired = F)
pairwise.wilcox.test(all_six,all_six_name,p.adjust.method = 'b',paired = F)
pairwise.wilcox.test(all_fif,all_fif_name,p.adjust.method = 'b',paired = F)
pairwise.wilcox.test(all_fif_rank,all_fif_rank_name,p.adjust.method = 'b',paired = F)

# 
