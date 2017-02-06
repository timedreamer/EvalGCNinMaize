
# for arabidopsis, use all_rank only, disgard all.

# load Data
#1363
load("~/Co-expression/arabidopsis/1363/auc_PP_3Methods.RData")
pcc_1363 <- auc_ara_pcc;mrnet_1363 <- auc_ara_mrnet;clr_1363 <- auc_ara_clr

#all
# load("~/Co-expression/arabidopsis/Text_arab/auc_AggPP_3Methods.RData")

# all_rank
load("~/Co-expression/arabidopsis/agg_rank/arab_PPPTY_aggRank.RData")

#PCC
p_1363 <- unlist(unname(pcc_1363));names(p_1363) <- rep("1363",3131)
# p_all <- unlist(unname(auc_ara_pcc));names(p_all) <- rep("all",3131)
p_all_rank <- unlist(unname(auc_ara_pccRank));names(p_all_rank) <- rep("all_rank",3131)

total_pcc <- c(p_1363,p_all_rank)
total_name_pcc <- c(names(p_1363),names(p_all_rank))
total_name_pcc <- factor(total_name_pcc,level = c("1363","all_rank"))
method_pcc <- rep("pcc",6262) 
mean_pcc <- c(mean(p_1363),mean(p_all_rank))

#mrnet
p_1363 <- unlist(unname(mrnet_1363));names(p_1363) <- rep("1363",3131)
# p_all <- unlist(unname(auc_ara_mrnet));names(p_all) <- rep("all",3131)
p_all_rank <- unlist(unname(auc_ara_mrnetRank));names(p_all_rank) <- rep("all_rank",3131)

total_mrnet <- c(p_1363,p_all_rank)
total_name_mrnet <- c(names(p_1363),names(p_all_rank))
total_name_mrnet <- factor(total_name_mrnet,level = c("1363","all_rank"))
method_mrnet <- rep("mrnet",6262) 
mean_mrnet <- c(mean(p_1363),mean(p_all_rank))

#clr
p_1363 <- unlist(unname(clr_1363));names(p_1363) <- rep("1363",3131)
# p_all <- unlist(unname(auc_ara_clr));names(p_all) <- rep("all",3131)
p_all_rank <- unlist(unname(auc_ara_clrRank));names(p_all_rank) <- rep("all_rank",3131)

total_clr <- c(p_1363,p_all_rank)
total_name_clr <- c(names(p_1363),names(p_all_rank))
total_name_clr <- factor(total_name_clr,level = c("1363","all_rank"))
method_clr <- rep("clr",6262) 
mean_clr <- c(mean(p_1363),mean(p_all_rank))


# one-factor boxplot. export as svg, 500*500
par(mar=c(9,4,4,1),mfrow=c(1,3))
boxplot(total_pcc~total_name_pcc,col=c("white","grey"),
        main="PP_PCC",las=2,ylim=c(0.4,0.8),
        outcol="grey",whiskcol="grey",staplecol="grey",cex.axis=2,
        names=c("1363","all_rank"));box(lwd=2)
points(mean_pcc,col="black",pch=8)
abline(h=0.560484,lty=2)

boxplot(total_mrnet~total_name_mrnet,main="PP_MRNET",las=2,ylim=c(0.4,0.8),
        col=c("white","grey"),
        outcol="grey",whiskcol="grey",staplecol="grey",cex.axis=2,
        names=c("1363","all_rank"));box(lwd=2)
points(mean_mrnet,col="black",pch=8)
abline(h=0.521885,lty=2)

boxplot(total_clr~total_name_clr,main="PP_CLR",las=2,ylim=c(0.4,0.8),
        col=c("white","grey"),
        outcol="grey",whiskcol="grey",staplecol="grey",cex.axis=2,
        names=c("1363","all_rank"));box(lwd=2)
points(mean_clr,col="black",pch=8)
abline(h=0.522989,lty=2)

# Plot everything together in boxplot .export as svg, 500*600
total_value <- c(total_pcc,total_mrnet,total_clr)
total_netname <- rep(c(rep("1363",3131),rep("all_rank",3131)),3)
total_constructMethod <- c(method_pcc,method_mrnet,method_clr)


total_netname <- factor(total_netname,c("1363","all_rank"))
total_constructMethod <- factor(total_constructMethod,c("pcc","mrnet","clr"))
mean_value_all <- c(0.5963725,0.6111365, 0.6144819,0.6135576,0.6186971,0.6204410)
# 
# par(mfrow=c(1,1),mar=c(7,4,4,2))
# boxplot(total_value~total_constructMethod*total_netname,las=2,yaxt="n",
#         col=rep(c(rep("white",3),rep("grey",3)),5),ylim=c(0.4,0.8),
#         outcol="grey",whiskcol="grey",staplecol="grey",main="allNtwk")
# box(lwd=2);axis(2,cex.axis=1.5,las=1)

par(mfrow=c(1,1),mar=c(9,4,4,2))
boxplot(total_value~total_constructMethod*total_netname,las=2,yaxt="n",
        col=c("white","grey","grey34"),ylim=c(0.4,0.8),cex.axis=1.5,
        outcol="grey",whiskcol="grey",staplecol="grey",main="allNtwk")
box(lwd=2);axis(2,cex.axis=1.5,las=1)
points(mean_value_all,col="black",pch=8)

######################################################################
#pairwise wilconxon test
pairwise.wilcox.test(total_pcc,total_name_pcc,p.adjust.method = "b",correct=F)
pairwise.wilcox.test(total_mrnet,total_name_mrnet,p.adjust.method = "b",correct=F)
pairwise.wilcox.test(total_clr,total_name_clr,p.adjust.method = "b",correct=F)


# compare 1363, all and all_rank networks among three methods.
p_pcc_1363 <- unlist(unname(pcc_1363));p_pcc_all <- unlist(unname(auc_ara_pcc))
p_pcc_all_rank <- unlist(unname(auc_ara_pccRank))

p_mrnet_1363 <- unlist(unname(mrnet_1363));p_mrnet_all <- unlist(unname(auc_ara_mrnet))
p_mrnet_all_rank <- unlist(unname(auc_ara_mrnetRank))

p_clr_1363 <- unlist(unname(clr_1363));p_clr_all <- unlist(unname(auc_ara_clr))
p_clr_all_rank <- unlist(unname(auc_ara_clrRank))


total_1363 <- c(p_pcc_1363,p_mrnet_1363,p_clr_1363)
all_1266_name <- c(rep("pcc",3131),rep("mrnet",3131),rep("clr",3131))

total_all <- c(p_pcc_all,p_mrnet_all,p_clr_all)
all_six_name <- c(rep("pcc",3131),rep("mrnet",3131),rep("clr",3131))

total_all_rank <- c(p_pcc_all_rank,p_mrnet_all_rank,p_clr_all_rank)
all_fif_name <- c(rep("pcc",3131),rep("mrnet",3131),rep("clr",3131))


pairwise.wilcox.test(total_1363,all_1266_name,p.adjust.method = 'b',paired = F,correct=F)
pairwise.wilcox.test(total_all,all_six_name,p.adjust.method = 'b',paired = F,correct=F)
pairwise.wilcox.test(total_all_rank,all_fif_name,p.adjust.method = 'b',paired = F,correct=F)
