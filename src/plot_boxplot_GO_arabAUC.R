#


#load data
# 1363 GO
load("~/Co-expression/arabidopsis/1363/arabidopsis_GO_clr.RData") # GO_arab_clr
load("~/Co-expression/arabidopsis/1363/arabidopsis_GO_pcc.RData") # GO_arab_pcc
load("~/Co-expression/arabidopsis/1363/arabidopsis_GO_mrnet.RData") # GO_arab_mrnet

# all GO
load("~/Co-expression/arabidopsis/Text_arab/arab_all_pcc/GO_arab_agg_pcc.RData") # GO_arab_agg_pcc
load("~/Co-expression/arabidopsis/Text_arab/arab_all_clr/GO_arab_agg_clr.RData") # GO_arab_agg_clr
load("~/Co-expression/arabidopsis/Text_arab/arab_all_mrnet/GO_arab_agg_mrnet.RData") # GO_arab_agg_mrnet

# allRank GO
load("~/Co-expression/arabidopsis/agg_rank/GO_arab_aggRank_clr.RData") # GO_ara_aggRank_clr
load("~/Co-expression/arabidopsis/agg_rank/GO_arab_aggRank_pcc.RData")  # GO_arab_aggRank_pcc
load("~/Co-expression/arabidopsis/agg_rank/GO_arab_aggRank_mrnet.RData") # GO_arab_aggRank_mrnet


#PCC
p_1363 <- GO_arab_pcc[[1]][,1];p_all <- GO_arab_agg_pcc[[1]][,1];p_allRank <- GO_arab_aggRank_pcc[[1]][,1]
total_pcc <- c(p_1363,p_all,p_allRank)
total_name_pcc <- c(rep("1363",445),rep("all",445),rep("all_rank",445))
total_name_pcc <- factor(total_name_pcc,level = c("1363","all","all_rank"))
method_pcc <- rep("pcc",1335) #
mean_pcc <- c(mean(p_1363),mean(p_all),mean(p_allRank))

#mrnet
p_1363 <- GO_arab_mrnet[[1]][,1];p_all <- GO_arab_agg_mrnet[[1]][,1];p_allRank <- GO_arab_aggRank_mrnet[[1]][,1]
total_mrnet <- c(p_1363,p_all,p_allRank)
total_name_mrnet <- c(rep("1363",445),rep("all",445),rep("all_rank",445))
total_name_mrnet <- factor(total_name_mrnet,level = c("1363","all","all_rank"))
method_mrnet <- rep("mrnet",1335) #
mean_mrnet <- c(mean(p_1363),mean(p_all),mean(p_allRank))

#clr
p_1363 <- GO_arab_clr[[1]][,1];p_all <- GO_arab_agg_clr[[1]][,1];p_allRank <- GO_arab_aggRank_clr[[1]][,1]
total_clr <- c(p_1363,p_all,p_allRank)
total_name_clr <- c(rep("1363",445),rep("all",445),rep("all_rank",445))
total_name_clr <- factor(total_name_clr,level = c("1363","all","all_rank"))
method_clr <- rep("clr",1335) #
mean_clr <- c(mean(p_1363),mean(p_all),mean(p_allRank))



# one-factor boxplot. export 8*6
par(mar=c(9,4,4,1),mfrow=c(1,3))
boxplot(total_pcc~total_name_pcc,main="GO_PCC",las=2,ylim=c(0.4,1),
        color=c("white"),
        outcol="grey",whiskcol="grey",staplecol="grey",cex.axis=2,
        names=c("1363","all","all_rank"));box(lwd=2)
points(mean_pcc,col="black",pch=8)

boxplot(total_mrnet~total_name_mrnet,main="GO_mrnet",las=2,ylim=c(0.4,1),
        color=c("white"),
        outcol="grey",whiskcol="grey",staplecol="grey",cex.axis=2,
        names=c("1363","all","all_rank"));box(lwd=2)
points(mean_mrnet,col="black",pch=8)

boxplot(total_clr~total_name_clr,main="GO_clr",las=2,ylim=c(0.4,1),
        color=c("white"),
        outcol="grey",whiskcol="grey",staplecol="grey",cex.axis=2,
        names=c("1363","all","all_rank"));box(lwd=2)
points(mean_clr,col="black",pch=8)

# Plot everything together in boxplot .export as 6*6
total_value <- c(total_pcc,total_mrnet,total_clr)
total_netname <- rep(c(rep("1363",445),rep("all",445),rep("all_rank",445)),3)
total_constructMethod <- c(method_pcc,method_mrnet,method_clr)


total_netname <- factor(total_netname,c("1363","all","all_rank"))
total_constructMethod <- factor(total_constructMethod,c("pcc","mrnet","clr"))


par(mfrow=c(1,1),mar=c(7,4,4,2))
boxplot(total_value~total_constructMethod*total_netname,las=2,yaxt="n",
        col=rep(c(rep("white",3),rep("grey",3)),5),ylim=c(0.4,1),
        outcol="grey",whiskcol="grey",staplecol="grey",main="allNtwk")
box(lwd=2);axis(2,cex.axis=1.5,las=1)

boxplot(total_value~total_constructMethod*total_netname,las=2,yaxt="n",
        col=c("white","grey","grey34"),ylim=c(0.4,1),
        outcol="grey",whiskcol="grey",staplecol="grey",main="allNtwk")
box(lwd=2);axis(2,cex.axis=1.5,las=1)


######################################################################
#pairwise wilconxon test
pairwise.wilcox.test(total_pcc,total_name_pcc,p.adjust.method = "b",correct=F)
pairwise.wilcox.test(total_mrnet,total_name_mrnet,p.adjust.method = "b",correct=F)
pairwise.wilcox.test(total_clr,total_name_clr,p.adjust.method = "b",correct=F)


# compare 1266, six and fifteen networks among three methods.
p_pcc_1363 <- GO_arab_pcc[[1]][,1];p_pcc_all <- GO_arab_agg_pcc[[1]][,1]
p_pcc_all_rank <- GO_arab_aggRank_pcc[[1]][,1]

p_mrnet_1363 <- GO_arab_mrnet[[1]][,1];p_mrnet_all <- GO_arab_agg_mrnet[[1]][,1]
p_mrnet_all_rank <- GO_arab_aggRank_mrnet[[1]][,1]

p_clr_1363 <- GO_arab_clr[[1]][,1];p_clr_all <- GO_arab_agg_clr[[1]][,1]
p_clr_all_rank <- GO_arab_aggRank_clr[[1]][,1]

total_1363 <- c(p_pcc_1363,p_mrnet_1363,p_clr_1363)
all_1266_name <- c(rep("pcc",445),rep("mrnet",445),rep("clr",445))

total_all <- c(p_pcc_all,p_mrnet_all,p_clr_all)
all_six_name <- c(rep("pcc",445),rep("mrnet",445),rep("clr",445))

total_all_rank <- c(p_pcc_all_rank,p_mrnet_all_rank,p_clr_all_rank)
all_fif_name <- c(rep("pcc",445),rep("mrnet",445),rep("clr",445))


pairwise.wilcox.test(total_1363,all_1266_name,p.adjust.method = 'b',paired = F,correct=F)
pairwise.wilcox.test(total_all,all_six_name,p.adjust.method = 'b',paired = F,correct=F)
pairwise.wilcox.test(total_all_rank,all_fif_name,p.adjust.method = 'b',paired = F,correct=F)
