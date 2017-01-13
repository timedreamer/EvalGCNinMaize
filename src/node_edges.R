#######################################################
#Calculate how many nodes and edged within each co-expression network.
# threshold can be changed by change calSumR function 0.7/0.8/0.9.
# Draw coefficiency distribution.
# Part I is for five correlation method. Part II is for MI.
#
setwd("D:\\Users\\jhuang\\Documents\\Co-expression\\Normalization result/")
library(poweRlaw)

threshold <- quantile(vst_pcc,0.9956235)

calsumR <- function(x){
  passThreshold <- length(which(x>threshold))
  return(passThreshold)
}

#VST
passT1 <- apply(vst_pcc,2,calsumR);passT1 <- passT1[which(passT1>1)]
passT2 <- apply(vst_scc,2,calsumR);passT2 <- passT2[which(passT2>1)]
passT3 <- apply(vst_kcc,2,calsumR);passT3 <- passT3[which(passT3>1)]
passT4 <- apply(vst_gcc,2,calsumR);passT4 <- passT4[which(passT4>1)]
passT5 <- apply(vst_bicor,2,calsumR);passT5 <- passT5[which(passT5>1)]


##################################################################################################
##################################################################################################
#export figure as 10*10.
setwd("D:\\Users\\jhuang\\Documents\\Co-expression\\Normalization result\\FromNewScript_Oct2016")

plcalc <- function(ntwk,top1M=0.9956235){
  threshold <- quantile(ntwk,top1M)
  calsumR <- function(x){
    passThreshold <- length(which(x>threshold))
    return(passThreshold)
  }
  passT1 <- apply(ntwk,2,calsumR);
  passT1 <- passT1[which(passT1>1)]
  m_pl <- displ$new(passT1);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
  bs_p <- bootstrap_p(m_pl,threads = 12,seed=115,no_of_sims = 1000)
  plot(m_pl,main=deparse(substitute(ntwk)),xlab="",ylab="",cex.lab=1.5);
  lines(m_pl,col=2,lwd=2);
  legend("bottomleft",legend = paste0("p-value: ",bs_p$p),bty="n",cex=1.5)
}

load("./vst_20161027/vst_pcc_scc.RData");load("./vst_20161027/vst_kcc.RData")
load("./vst_20161027/vst_bicor.RData");load("vst_20161027/vst_four_MI.RData")

par(mfrow=c(3,3),mar=c(6,6,3,1))

plcalc(vst_pcc);plcalc(vst_scc);plcalc(vst_kcc);plcalc(vst_gcc);
plcalc(vst_bicor);plcalc(vst_aa);plcalc(vst_ma);
plcalc(vst_mrnet);plcalc(vst_clr)

rm(vst_pcc,vst_scc,vst_kcc,vst_gcc,vst_aa,vst_ma,vst_mrnet,vst_clr,vst_bicor)




#CPM
load("./cpm_20161024/cpm_pcc_scc.RData");load("cpm_20161024/cpm_kcc.RData")
load("./cpm_20161024/cpm_bicor.RData");load("cpm_20161024/cpm_four_MI.RData")

par(mfrow=c(3,3),mar=c(6,6,3,1))

plcalc(cpm_pcc);plcalc(cpm_scc);plcalc(cpm_kcc);plcalc(cpm_gcc);
plcalc(cpm_bicor);plcalc(cpm_aa);plcalc(cpm_ma);
plcalc(cpm_mrnet);plcalc(cpm_clr)

rm(cpm_pcc,cpm_scc,cpm_kcc,cpm_gcc,cpm_aa,cpm_ma,cpm_mrnet,cpm_clr,cpm_bicor)


#RPKM
load("rpkm_20161027/rpkm_pcc_scc.RData");load("rpkm_20161027/rpkm_kcc.RData")
load("rpkm_20161027/rpkm_bicor.RData");load("rpkm_20161027/rpkm_four_MI.RData")

par(mfrow=c(3,3),mar=c(6,6,3,1))

plcalc(rpkm_pcc);plcalc(rpkm_scc);plcalc(rpkm_kcc);plcalc(rpkm_gcc);
plcalc(rpkm_bicor);plcalc(rpkm_aa);plcalc(rpkm_ma);
plcalc(rpkm_mrnet);plcalc(rpkm_clr)

rm(rpkm_pcc,rpkm_scc,rpkm_kcc,rpkm_gcc,rpkm_aa,rpkm_ma,rpkm_mrnet,rpkm_clr,rpkm_bicor)

#RC
load("rc_20161027/rc_pcc_scc.RData");load("rc_20161027/rc_kcc.RData")
load("rc_20161027/rc_bicor.RData");load("rc_20161027/rc_four_MI.RData")

par(mfrow=c(3,3),mar=c(6,6,3,1))

plcalc(rc_pcc);plcalc(rc_scc);plcalc(rc_kcc);plcalc(rc_gcc);
plcalc(rc_bicor);plcalc(rc_aa);plcalc(rc_ma);
plcalc(rc_mrnet);plcalc(rc_clr)

rm(rc_pcc,rc_scc,rc_kcc,rc_gcc,rc_aa,rc_ma,rc_mrnet,rc_clr,rc_bicor)

#################################################################
#################################################################

threshold <- quantile(vst_pcc,top1M)
calsumR <- function(x){
  passThreshold <- length(which(x>threshold))
  return(passThreshold)
}

passT1 <- apply(vst_pcc,2,calsumR);
passT1 <- passT1[which(passT1>1)]
m_pl <- displ$new(passT1);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
bs_p <- bootstrap_p(m_pl,threads = 12,seed=115)
plot(m_pl,main=deparse(substitute(ntwk)),xlab="",ylab="",cex.lab=1.5);
legend("bottomleft",legend = paste0("p-value: ",bs_p$p),bty="n",cex=1.5)
lines(m_pl,col=2,lwd=2)





#total node
vst_node <- c(length(passT1),length(passT2),length(passT3),length(passT4),length(passT5))
#total edges
vst_edges <- c(sum(passT1),sum(passT2),sum(passT3),sum(passT4),sum(passT5))

par(mfrow=c(2,3))
.hist1 <- hist(passT1,main="vst_pcc")
.hist2 <- hist(passT2,main="vst_scc")
.hist3 <- hist(passT3,main="vst_kcc")
.hist4 <- hist(passT4,main="vst_gcc")
.hist5 <- hist(passT5,main="vst_bicor")
par(mfrow=c(2,3))
with(.hist1,plot(mids,counts,main="vst_pcc"))
with(.hist2,plot(mids,counts,main="vst_scc"))
with(.hist3,plot(mids,counts,main="vst_kcc"))
with(.hist4,plot(mids,counts,main="vst_gcc"))
with(.hist5,plot(mids,counts,main="vst_bicor"))
max(passT1)

par(mfrow =c(2,3))
m_pl <- displ$new(passT1);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="vst_pcc",ylab="CDF");lines(m_pl,col=2,lwd=2)
#bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p # 0

m_pl <- displ$new(passT2);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="vst_scc",ylab="CDF");lines(m_pl,col=2,lwd=2)
#bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p # 0

m_pl <- displ$new(passT3);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="vst_kcc",ylab="CDF");lines(m_pl,col=2,lwd=2)
#bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p # 0

m_pl <- displ$new(passT4);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="vst_gcc",ylab="CDF");lines(m_pl,col=2,lwd=2)
#bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p # 0
m_pl <- displ$new(passT5);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="vst_bicor",ylab="CDF");lines(m_pl,col=2,lwd=2)
#bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p # 0


#CPM
passT1 <- apply(cpm_pcc,2,calsumR);passT1 <- passT1[which(passT1>1)]
passT2 <- apply(cpm_scc,2,calsumR);passT2 <- passT2[which(passT2>1)]
passT3 <- apply(cpm_kcc,2,calsumR);passT3 <- passT3[which(passT3>1)]
passT4 <- apply(cpm_gcc,2,calsumR);passT4 <- passT4[which(passT4>1)]
passT5 <- apply(cpm_bicor,2,calsumR);passT5 <- passT5[which(passT5>1)]

#total node
cpm_node <- c(length(passT1),length(passT2),length(passT3),length(passT4),length(passT5))
#total edges
cpm_edges <- c(sum(passT1),sum(passT2),sum(passT3),sum(passT4),sum(passT5))

par(mfrow=c(2,3))
.hist1 <- hist(passT1,main="cpm_pcc")
.hist2 <- hist(passT2,main="cpm_scc")
.hist3 <- hist(passT3,main="cpm_kcc")
.hist4 <- hist(passT4,main="cpm_gcc")
.hist5 <- hist(passT5,main="cpm_bicor")
par(mfrow=c(2,3))
with(.hist1,plot(mids,counts,main="cpm_pcc"))
with(.hist2,plot(mids,counts,main="cpm_scc"))
with(.hist3,plot(mids,counts,main="cpm_kcc"))
with(.hist4,plot(mids,counts,main="cpm_gcc"))
with(.hist5,plot(mids,counts,main="cpm_bicor"))
max(passT1)

#RawCount
passT1 <- apply(rc_pcc,2,calsumR);passT1 <- passT1[which(passT1>1)]
passT2 <- apply(rc_scc,2,calsumR);passT2 <- passT2[which(passT2>1)]
passT3 <- apply(rawcount_kcc,2,calsumR);passT3 <- passT3[which(passT3>1)]
passT4 <- apply(rawcount_gcc,2,calsumR);passT4 <- passT4[which(passT4>1)]
passT5 <- apply(rc_bicor,2,calsumR);passT5 <- passT5[which(passT5>1)]

#total node
rc_node <- c(length(passT1),length(passT2),length(passT3),length(passT4),length(passT5))
#total edges
rc_edges <- c(sum(passT1),sum(passT2),sum(passT3),sum(passT4),sum(passT5))

par(mfrow=c(2,3))
.hist1 <- hist(passT1,main="rc_pcc")
.hist2 <- hist(passT2,main="rc_scc")
.hist3 <- hist(passT3,main="rc_kcc")
.hist4 <- hist(passT4,main="rc_gcc")
.hist5 <- hist(passT5,main="rc_bicor")
par(mfrow=c(2,3))
with(.hist1,plot(mids,counts,main="rc_pcc"))
with(.hist2,plot(mids,counts,main="rc_scc"))
with(.hist3,plot(mids,counts,main="rc_kcc"))
with(.hist4,plot(mids,counts,main="rc_gcc"))
with(.hist5,plot(mids,counts,main="rc_bicor"))
max(passT1)

par(mfrow =c(2,3))
m_pl <- displ$new(passT1);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="rc_pcc",ylab="CDF");lines(m_pl,col=2,lwd=2)
#bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p # 0

m_pl <- displ$new(passT2);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="rc_scc",ylab="CDF");lines(m_pl,col=2,lwd=2)
#bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p # 0

m_pl <- displ$new(passT3);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="rc_kcc",ylab="CDF");lines(m_pl,col=2,lwd=2)
#bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p # 0

m_pl <- displ$new(passT4);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="rc_gcc",ylab="CDF");lines(m_pl,col=2,lwd=2)
#bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p # 0
m_pl <- displ$new(passT5);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="rc_bicor",ylab="CDF");lines(m_pl,col=2,lwd=2)
#bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p # 0

#RPKM
passT1 <- apply(rpkm_pcc,2,calsumR);passT1 <- passT1[which(passT1>1)]
passT2 <- apply(rpkm_scc,2,calsumR);passT2 <- passT2[which(passT2>1)]
passT3 <- apply(rpkm_kcc,2,calsumR);passT3 <- passT3[which(passT3>1)]
passT4 <- apply(rpkm_gcc,2,calsumR);passT4 <- passT4[which(passT4>1)]
passT5 <- apply(rpkm_bicor,2,calsumR);passT5 <- passT5[which(passT5>1)]

#total node
rpkm_node <- c(length(passT1),length(passT2),length(passT3),length(passT4),length(passT5))
#total edges
rpkm_edges <- c(sum(passT1),sum(passT2),sum(passT3),sum(passT4),sum(passT5))

par(mfrow=c(2,3))
.hist1 <- hist(passT1,main="rpkm_pcc")
.hist2 <- hist(passT2,main="rpkm_scc")
.hist3 <- hist(passT3,main="rpkm_kcc")
.hist4 <- hist(passT4,main="rpkm_gcc")
.hist5 <- hist(passT5,main="rpkm_bicor")
par(mfrow=c(2,3))
with(.hist1,plot(mids,counts,main="rpkm_pcc"))
with(.hist2,plot(mids,counts,main="rpkm_scc"))
with(.hist3,plot(mids,counts,main="rpkm_kcc"))
with(.hist4,plot(mids,counts,main="rpkm_gcc"))
with(.hist5,plot(mids,counts,main="rpkm_bicor"))
max(passT1)

#################################################################

#Draw histgram distribution for all normalization/correlation method
par(mfrow=c(4,5))
hist(vst_pcc,probability = T,xlim=c(-1,1),main="vst_pcc");curve(dnorm(x, mean=mean(vst_pcc), sd=sd(vst_pcc)), add=TRUE,col=2,lwd=2)
hist(vst_scc,probability = T,xlim=c(-1,1),main="vst_scc");curve(dnorm(x, mean=mean(vst_scc), sd=sd(vst_scc)), add=TRUE,col=2,lwd=2)
hist(vst_kcc,probability = T,xlim=c(-1,1),main="vst_kcc");curve(dnorm(x, mean=mean(vst_kcc), sd=sd(vst_kcc)), add=TRUE,col=2,lwd=2)
hist(vst_gcc,probability = T,xlim=c(-1,1),main="vst_gcc");curve(dnorm(x, mean=mean(vst_gcc), sd=sd(vst_gcc)), add=TRUE,col=2,lwd=2)
hist(vst_bicor,probability = T,xlim=c(-1,1),main="vst_bicor_distribution");curve(dnorm(x, mean=mean(vst_bicor), sd=sd(vst_bicor)), add=TRUE,col=2,lwd=2)

hist(cpm_pcc,probability = T,xlim=c(-1,1),main="cpm_pcc");curve(dnorm(x, mean=mean(cpm_pcc), sd=sd(cpm_pcc)), add=TRUE,col=2,lwd=2)
hist(cpm_scc,probability = T,xlim=c(-1,1),main="cpm_scc");curve(dnorm(x, mean=mean(cpm_scc), sd=sd(cpm_scc)), add=TRUE,col=2,lwd=2)
hist(cpm_kcc,probability = T,xlim=c(-1,1),main="cpm_kcc");curve(dnorm(x, mean=mean(cpm_kcc), sd=sd(cpm_kcc)), add=TRUE,col=2,lwd=2)
hist(cpm_gcc,probability = T,xlim=c(-1,1),main="cpm_gcc");curve(dnorm(x, mean=mean(cpm_gcc), sd=sd(cpm_gcc)), add=TRUE,col=2,lwd=2)
hist(cpm_bicor,probability = T,xlim=c(-1,1),main="cpm_bicor");curve(dnorm(x, mean=mean(cpm_bicor), sd=sd(cpm_bicor)), add=TRUE,col=2,lwd=2)

hist(rc_pcc,probability = T,xlim=c(-1,1),main="rc_pcc");curve(dnorm(x, mean=mean(rc_pcc), sd=sd(rc_pcc)), add=TRUE,col=2,lwd=2)
hist(rc_scc,probability = T,xlim=c(-1,1),main="rc_scc");curve(dnorm(x, mean=mean(rc_scc), sd=sd(rc_scc)), add=TRUE,col=2,lwd=2)
hist(rawcount_kcc,probability = T,xlim=c(-1,1),main="rc_kcc");curve(dnorm(x, mean=mean(rawcount_kcc), sd=sd(rawcount_kcc)), add=TRUE,col=2,lwd=2)
hist(rawcount_gcc,probability = T,xlim=c(-1,1),main="rc_gcc");curve(dnorm(x, mean=mean(rawcount_gcc), sd=sd(rawcount_gcc)), add=TRUE,col=2,lwd=2)
hist(rc_bicor,probability = T,xlim=c(-1,1),main="rc_bicor");curve(dnorm(x, mean=mean(rc_bicor), sd=sd(rc_bicor)), add=TRUE,col=2,lwd=2)

hist(rpkm_pcc,probability = T,xlim=c(-1,1),main="rpkm_pcc");curve(dnorm(x, mean=mean(rpkm_pcc), sd=sd(rpkm_pcc)), add=TRUE,col=2,lwd=2)
hist(rpkm_scc,probability = T,xlim=c(-1,1),main="rpkm_scc");curve(dnorm(x, mean=mean(rpkm_scc), sd=sd(rpkm_scc)), add=TRUE,col=2,lwd=2)
hist(rpkm_kcc,probability = T,xlim=c(-1,1),main="rpkm_kcc");curve(dnorm(x, mean=mean(rpkm_kcc), sd=sd(rpkm_kcc)), add=TRUE,col=2,lwd=2)
hist(rpkm_gcc,probability = T,xlim=c(-1,1),main="rpkm_gcc");curve(dnorm(x, mean=mean(rpkm_gcc), sd=sd(rpkm_gcc)), add=TRUE,col=2,lwd=2)
hist(rpkm_bicor,probability = T,xlim=c(-1,1),main="rpkm_bicor");curve(dnorm(x, mean=mean(rpkm_bicor), sd=sd(rpkm_bicor)), add=TRUE,col=2,lwd=2)
########Try power law distribution



##################################PART II###############################################
# node and edges for MI method
library(poweRlaw)
# find 80%/90% quantile of aaracene method. all zero
quantile(vst_aaracne,0.9);quantile(cpm_aaracne,0.9);quantile(rc_aaracne,0.9);quantile(rpkm_aaracne,0.9)


# set calsumR threshold to 0
passT1 <- apply(vst_aaracne,2,calsumR);passT1 <- passT1[which(passT1>1)]
passT2 <- apply(cpm_aaracne,2,calsumR);passT2 <- passT2[which(passT2>1)]
passT3 <- apply(rc_aaracne,2,calsumR);passT3 <- passT3[which(passT3>1)]
passT4 <- apply(rpkm_aaracne,2,calsumR);passT4 <- passT4[which(passT4>1)]


#total node
aa_node <- c(length(passT1),length(passT2),length(passT3),length(passT4))
#total edges
aa_edges <- c(sum(passT1),sum(passT2),sum(passT3),sum(passT4))

# draw powerlaw graph and calculate p value.
# p-value >0.05/0.01 powerlaw
par(mfrow =c(2,2))
m_pl <- displ$new(passT1);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="vst_aa",ylab="CDF");lines(m_pl,col=2,lwd=2)
bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p #0.06

m_pl <- displ$new(passT2);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="cpm_aa",ylab="CDF");lines(m_pl,col=2,lwd=2)
bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p #0.02

m_pl <- displ$new(passT3);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="rc_aa",ylab="CDF");lines(m_pl,col=2,lwd=2)
bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p #0

m_pl <- displ$new(passT4);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="rpkm_aa",ylab="CDF");lines(m_pl,col=2,lwd=2)
bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p #0.02



###MARACENE
# find 80%/90% quantile of maracene method. all zero
quantile(vst_maracne,0.9);quantile(cpm_maracne,0.9);quantile(rc_maracne,0.9);quantile(rpkm_maracne,0.9)

# set calsumR threshold to 0
passT1 <- apply(vst_maracne,2,calsumR);passT1 <- passT1[which(passT1>1)]
passT2 <- apply(cpm_maracne,2,calsumR);passT2 <- passT2[which(passT2>1)]
passT3 <- apply(rc_maracne,2,calsumR);passT3 <- passT3[which(passT3>1)]
passT4 <- apply(rpkm_maracne,2,calsumR);passT4 <- passT4[which(passT4>1)]


#total node
ma_node <- c(length(passT1),length(passT2),length(passT3),length(passT4))
#total edges
ma_edges <- c(sum(passT1),sum(passT2),sum(passT3),sum(passT4))

ma_node;ma_edges
# draw powerlaw graph and calculate p value.
# p-value >0.05/0.01 powerlaw
par(mfrow =c(2,2))
m_pl <- displ$new(passT1);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="vst_ma",ylab="CDF");lines(m_pl,col=2,lwd=2)
bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p # 0

m_pl <- displ$new(passT2);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="cpm_ma",ylab="CDF");lines(m_pl,col=2,lwd=2)
bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p # 0.02

m_pl <- displ$new(passT3);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="rc_ma",ylab="CDF");lines(m_pl,col=2,lwd=2)
bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p # 0

m_pl <- displ$new(passT4);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="rpkm_ma",ylab="CDF");lines(m_pl,col=2,lwd=2)
bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p # 0.05

###CLR
# find 90% quantile of maracene method. 3.1,3.1,2.83,3.1
quantile(vst_clr,0.99);quantile(cpm_clr,0.99);quantile(rc_clr,0.99);quantile(rpkm_clr,0.99)
#Then try 99%. 9.14,9.15,8.92,9.15


# set calsumR threshold to 0
passT1 <- apply(vst_clr,2,calsumR);passT1 <- passT1[which(passT1>1)]
passT2 <- apply(cpm_clr,2,calsumR);passT2 <- passT2[which(passT2>1)]
passT3 <- apply(rc_clr,2,calsumR);passT3 <- passT3[which(passT3>1)]
passT4 <- apply(rpkm_clr,2,calsumR);passT4 <- passT4[which(passT4>1)]


#total node
clr_node <- c(length(passT1),length(passT2),length(passT3),length(passT4))
#total edges
clr_edges <- c(sum(passT1),sum(passT2),sum(passT3),sum(passT4))

clr_node;clr_edges
# draw powerlaw graph and calculate p value.
# p-value >0.05/0.01 powerlaw
par(mfrow =c(2,2))
m_pl <- displ$new(passT1);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="vst_clr",ylab="CDF");lines(m_pl,col=2,lwd=2)
bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p # 0

m_pl <- displ$new(passT2);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="cpm_clr",ylab="CDF");lines(m_pl,col=2,lwd=2)
bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p # 0

m_pl <- displ$new(passT3);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="rc_clr",ylab="CDF");lines(m_pl,col=2,lwd=2)
bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p # 0

m_pl <- displ$new(passT4);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="rpkm_clr",ylab="CDF");lines(m_pl,col=2,lwd=2)
bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p # 0

###MRNET
# find 90% quantile of maracene method.0.092,0.092,0.108,0.092
quantile(vst_mrnet,0.99);quantile(cpm_mrnet,0.99);quantile(rc_mrnet,0.99);quantile(rpkm_mrnet,0.99)
#then choose 99%. 0.153,0.153,0.165,0.153


# set calsumR threshold to 0
passT1 <- apply(vst_mrnet,2,calsumR);passT1 <- passT1[which(passT1>1)]
passT2 <- apply(cpm_mrnet,2,calsumR);passT2 <- passT2[which(passT2>1)]
passT3 <- apply(rc_mrnet,2,calsumR);passT3 <- passT3[which(passT3>1)]
passT4 <- apply(rpkm_mrnet,2,calsumR);passT4 <- passT4[which(passT4>1)]


#total node
mrnet_node <- c(length(passT1),length(passT2),length(passT3),length(passT4))
#total edges
mrnet_edges <- c(sum(passT1),sum(passT2),sum(passT3),sum(passT4))

mrnet_node;mrnet_edges
# draw powerlaw graph and calculate p value.
# p-value >0.05/0.01 powerlaw
par(mfrow =c(2,2))
m_pl <- displ$new(passT1);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="vst_mrnet",ylab="CDF");lines(m_pl,col=2,lwd=2)
bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p # 0

m_pl <- displ$new(passT2);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="cpm_mrnet",ylab="CDF");lines(m_pl,col=2,lwd=2)
#bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p # 0

m_pl <- displ$new(passT3);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="rc_mrnet",ylab="CDF");lines(m_pl,col=2,lwd=2)
#bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p # 0

m_pl <- displ$new(passT4);est_pl <- estimate_xmin(m_pl);m_pl$setXmin(est_pl)
plot(m_pl,main="rpkm_mrnet",ylab="CDF");lines(m_pl,col=2,lwd=2)
#bs_p <- bootstrap_p(m_pl,threads = 10);bs_p$p # 0


#########################
#Draw node and edges
eg_pcc_vst <- c(1307754,227642,22162);eg_pcc_cpm <- c(1307751,242048,27032)
eg_pcc_rc <- c(50656865,14121756,531530);eg_pcc_rpkm <- c(1409595,265786,24552)

eg_scc_vst <- c(1384474,242048,27032);eg_scc_cpm <- c(1384474,255214,27776)
eg_scc_rc <- c(56304374,17013678,690364);eg_scc_rpkm <- c(1401796,255214,27772)

eg_kcc_vst <- c(29545,16322,15116);eg_kcc_cpm <- c(32386,16308,15116)
eg_kcc_rc <- c(2170349,33678,15126);eg_kcc_rpkm <- c(32370,16306,15116)

eg_gcc_vst <- c(2429413,470488,38976);eg_gcc_cpm <- c(2426413,467312,39302)
eg_gcc_rc <- c(66177987,22248960,1154376);eg_gcc_rpkm <- c(2343299,467306,39298)

eg_bi_vst <- c(1155467,183476,20584);eg_bi_cpm <- c(1254346,221928,24804)
eg_bi_rc <- c(52591528,15225462,576968);eg_bi_rpkm <- c(1254346,221928,24804)

par(mfrow=c(2,3))
plot(log2(eg_pcc_vst),type="b",col="1",ylim=c(10,30),main="edges_pcc",xaxt="n",xlab="PCC",ylab="log2(edges)")
axis(1, at=1:3, labels=c(0.7,0.8,0.9))
legend("topright",c("VST","CPM","RC","RPKM"),lwd=c(2,2,2),col=c(1,2,3,4),cex=0.6)
lines(log2(eg_pcc_cpm),type="b",col="2")
lines(log2(eg_pcc_rpkm),type="b",col="3",lwd=2)
lines(log2(eg_pcc_rc),type="b",col="4")

plot(log2(eg_scc_vst),type="b",col="1",ylim=c(10,30),main="edges_scc",xaxt="n",xlab="SCC",ylab="log2(edges)")
axis(1, at=1:3, labels=c(0.7,0.8,0.9))
legend("topright",c("VST","CPM","RC","RPKM"),lwd=c(2,2,2),col=c(1,2,3,4),cex=0.6)
lines(log2(eg_scc_cpm),type="b",col="2")
lines(log2(eg_scc_rpkm),type="b",col="3",lwd=2)
lines(log2(eg_scc_rc),type="b",col="4")

plot(log2(eg_kcc_vst),type="b",col="1",ylim=c(10,30),main="edges_kcc",xaxt="n",xlab="KCC",ylab="log2(edges)")
axis(1, at=1:3, labels=c(0.7,0.8,0.9))
legend("topright",c("VST","CPM","RC","RPKM"),lwd=c(2,2,2),col=c(1,2,3,4),cex=0.6)
lines(log2(eg_kcc_cpm),type="b",col="2")
lines(log2(eg_kcc_rpkm),type="b",col="3",lwd=2)
lines(log2(eg_kcc_rc),type="b",col="4")

plot(log2(eg_gcc_vst),type="b",col="1",ylim=c(10,30),main="edges_gcc",xaxt="n",xlab="GCC",ylab="log2(edges)")
axis(1, at=1:3, labels=c(0.7,0.8,0.9))
legend("topright",c("VST","CPM","RC","RPKM"),lwd=c(2,2,2),col=c(1,2,3,4),cex=0.6)
lines(log2(eg_gcc_cpm),type="b",col="2")
lines(log2(eg_gcc_rpkm),type="b",col="3",lwd=2)
lines(log2(eg_gcc_rc),type="b",col="4")

plot(log2(eg_bi_vst),type="b",col="1",ylim=c(10,30),main="edges_bi",xaxt="n",xlab="Bicor",ylab="log2(edges)")
axis(1, at=1:3, labels=c(0.7,0.8,0.9))
legend("topright",c("VST","CPM","RC","RPKM"),lwd=c(2,2,2),col=c(1,2,3,4),cex=0.6)
lines(log2(eg_bi_cpm),type="b",col="2")
lines(log2(eg_bi_rpkm),type="b",col="3",lwd=2)
lines(log2(eg_bi_rc),type="b",col="4")
