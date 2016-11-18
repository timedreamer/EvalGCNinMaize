# plot correlation matrix distribution. Export as pdf 10*10

setwd("D:\\Users\\jhuang\\Documents\\Co-expression\\Normalization result\\FromNewScript_Oct2016")

# load data
load("./vst_20161027/vst_pcc_scc.RData");load("./vst_20161027/vst_kcc.RData")
load("./vst_20161027/vst_bicor.RData");load("vst_20161027/vst_four_MI.RData")
#load(GCC) from old folder

load("./cpm_20161024/cpm_pcc_scc.RData");load("cpm_20161024/cpm_kcc.RData")
load("./cpm_20161024/cpm_bicor.RData");load("cpm_20161024/cpm_four_MI.RData")
load("~/Co-expression/Normalization result/Before_Oct2016/cpm/cpm_gcc.RData")

load("rpkm_20161027/rpkm_pcc_scc.RData");load("rpkm_20161027/rpkm_kcc.RData")
load("rpkm_20161027/rpkm_bicor.RData");load("rpkm_20161027/rpkm_four_MI.RData")
load("~/Co-expression/Normalization result/Before_Oct2016/rpkm/rpkm_gcc.RData")

load("rc_20161027/rc_pcc_scc.RData");load("rc_20161027/rc_kcc.RData")
load("rc_20161027/rc_bicor.RData");load("rc_20161027/rc_four_MI.RData")
#load(GCC)


########################################################################################
################################ VST ###################################################
########################################################################################


min(vst_aa);max(vst_aa);min(vst_ma);max(vst_ma); # decide the ylim
min(vst_mrnet);max(vst_mrnet);min(vst_clr);max(vst_clr)

#Draw histgram distribution for all normalization/correlation method
###################VST##############################
par(mfrow=c(3,3))
hist(vst_pcc,probability = T,xlim=c(-1,1),main="VST_PCC",cex.main=2,xlab="",cex.axis=1.5);
curve(dnorm(x, mean=mean(vst_pcc), sd=sd(vst_pcc)), add=TRUE,col=2,lwd=2)

hist(vst_scc,probability = T,xlim=c(-1,1),main="VST_SCC",cex.main=2,xlab="",cex.axis=1.5);
curve(dnorm(x, mean=mean(vst_scc), sd=sd(vst_scc)), add=TRUE,col=2,lwd=2)

hist(vst_kcc,probability = T,xlim=c(-1,1),main="VST_KCC",cex.main=2,xlab="",cex.axis=1.5);
curve(dnorm(x, mean=mean(vst_kcc), sd=sd(vst_kcc)), add=TRUE,col=2,lwd=2)

hist(vst_gcc,probability = T,xlim=c(-1,1),main="VST_GCC",cex.main=2,xlab="",cex.axis=1.5);
curve(dnorm(x, mean=mean(vst_gcc), sd=sd(vst_gcc)), add=TRUE,col=2,lwd=2)

hist(vst_bicor,probability = T,xlim=c(-1,1),main="VST_Bicor",cex.main=2,xlab="",cex.axis=1.5);
curve(dnorm(x, mean=mean(vst_bicor), sd=sd(vst_bicor)), add=TRUE,col=2,lwd=2)

hist(vst_aa,probability = T,xlim=c(0,1.784),main="VST_AA",cex.main=2,xlab="",cex.axis=1.5);

hist(vst_ma,probability = T,xlim=c(0,1.784),main="VST_MA",cex.main=2,xlab="",cex.axis=1.5);

hist(vst_mrnet,probability = T,xlim=c(0,1.784),main="VST_MRNET",cex.main=2,xlab="",cex.axis=1.5);

hist(vst_clr,probability = T,xlim=c(0,416),main="VST_CLR",cex.main=2,xlab="",cex.axis=1.5);

rm(vst_pcc,vst_scc,vst_kcc,vst_gcc,vst_aa,vst_ma,vst_mrnet,vst_clr,vst_bicor)

########################################################################################
################################ CPM ###################################################
########################################################################################

min(cpm_aa);max(cpm_aa);min(cpm_ma);max(cpm_ma);
min(cpm_mrnet);max(cpm_mrnet);min(cpm_clr);max(cpm_clr)

par(mfrow=c(3,3))
hist(cpm_pcc,probability = T,xlim=c(-1,1),main="cpm_PCC",cex.main=2,xlab="",cex.axis=1.5);
curve(dnorm(x, mean=mean(cpm_pcc), sd=sd(cpm_pcc)), add=TRUE,col=2,lwd=2)

hist(cpm_scc,probability = T,xlim=c(-1,1),main="cpm_SCC",cex.main=2,xlab="",cex.axis=1.5);
curve(dnorm(x, mean=mean(cpm_scc), sd=sd(cpm_scc)), add=TRUE,col=2,lwd=2)

hist(cpm_kcc,probability = T,xlim=c(-1,1),main="cpm_KCC",cex.main=2,xlab="",cex.axis=1.5);
curve(dnorm(x, mean=mean(cpm_kcc), sd=sd(cpm_kcc)), add=TRUE,col=2,lwd=2)

hist(cpm_gcc,probability = T,xlim=c(-1,1),main="cpm_GCC",cex.main=2,xlab="",cex.axis=1.5);
curve(dnorm(x, mean=mean(cpm_gcc), sd=sd(cpm_gcc)), add=TRUE,col=2,lwd=2)

hist(cpm_bicor,probability = T,xlim=c(-1,1),main="cpm_Bicor",cex.main=2,xlab="",cex.axis=1.5);
curve(dnorm(x, mean=mean(cpm_bicor), sd=sd(cpm_bicor)), add=TRUE,col=2,lwd=2)

hist(cpm_aa,probability = T,xlim=c(0,1.789),main="cpm_AA",cex.main=2,xlab="",cex.axis=1.5);

hist(cpm_ma,probability = T,xlim=c(0,1.789),main="cpm_MA",cex.main=2,xlab="",cex.axis=1.5);

hist(cpm_mrnet,probability = T,xlim=c(0,1.789),main="cpm_MRNET",cex.main=2,xlab="",cex.axis=1.5);

hist(cpm_clr,probability = T,xlim=c(0,435),main="cpm_CLR",cex.main=2,xlab="",cex.axis=1.5);

rm(cpm_pcc,cpm_scc,cpm_kcc,cpm_gcc,cpm_aa,cpm_ma,cpm_mrnet,cpm_clr,cpm_bicor)



########################################################################################
################################ RPKM ##################################################
########################################################################################

min(rpkm_aa);max(rpkm_aa);min(rpkm_ma);max(rpkm_ma);
min(rpkm_mrnet);max(rpkm_mrnet);min(rpkm_clr);max(rpkm_clr)



par(mfrow=c(3,3))
hist(rpkm_pcc,probability = T,xlim=c(-1,1),main="rpkm_PCC",cex.main=2,xlab="",cex.axis=1.5);
curve(dnorm(x, mean=mean(rpkm_pcc), sd=sd(rpkm_pcc)), add=TRUE,col=2,lwd=2)

hist(rpkm_scc,probability = T,xlim=c(-1,1),main="rpkm_SCC",cex.main=2,xlab="",cex.axis=1.5);
curve(dnorm(x, mean=mean(rpkm_scc), sd=sd(rpkm_scc)), add=TRUE,col=2,lwd=2)

hist(rpkm_kcc,probability = T,xlim=c(-1,1),main="rpkm_KCC",cex.main=2,xlab="",cex.axis=1.5);
curve(dnorm(x, mean=mean(rpkm_kcc), sd=sd(rpkm_kcc)), add=TRUE,col=2,lwd=2)

hist(rpkm_gcc,probability = T,xlim=c(-1,1),main="rpkm_GCC",cex.main=2,xlab="",cex.axis=1.5);
curve(dnorm(x, mean=mean(rpkm_gcc), sd=sd(rpkm_gcc)), add=TRUE,col=2,lwd=2)

hist(rpkm_bicor,probability = T,xlim=c(-1,1),main="rpkm_Bicor",cex.main=2,xlab="",cex.axis=1.5);
curve(dnorm(x, mean=mean(rpkm_bicor), sd=sd(rpkm_bicor)), add=TRUE,col=2,lwd=2)

hist(rpkm_aa,probability = T,xlim=c(0,1.787),main="rpkm_AA",cex.main=2,xlab="",cex.axis=1.5);

hist(rpkm_ma,probability = T,xlim=c(0,1.787),main="rpkm_MA",cex.main=2,xlab="",cex.axis=1.5);

hist(rpkm_mrnet,probability = T,xlim=c(0,1.787),main="rpkm_MRNET",cex.main=2,xlab="",cex.axis=1.5);

hist(rpkm_clr,probability = T,xlim=c(0,437.9),main="rpkm_CLR",cex.main=2,xlab="",cex.axis=1.5);

rm(rpkm_pcc,rpkm_scc,rpkm_kcc,rpkm_gcc,rpkm_aa,rpkm_ma,rpkm_mrnet,rpkm_clr,rpkm_bicor)


########################################################################################
################################ RC  ###################################################
########################################################################################

min(rc_aa);max(rc_aa);min(rc_ma);max(rc_ma);
min(rc_mrnet);max(rc_mrnet);min(rc_clr);max(rc_clr)

par(mfrow=c(3,3))
hist(rc_pcc,probability = T,xlim=c(-1,1),main="rc_PCC",cex.main=2,xlab="",cex.axis=1.5);
curve(dnorm(x, mean=mean(rc_pcc), sd=sd(rc_pcc)), add=TRUE,col=2,lwd=2)

hist(rc_scc,probability = T,xlim=c(-1,1),main="rc_SCC",cex.main=2,xlab="",cex.axis=1.5);
curve(dnorm(x, mean=mean(rc_scc), sd=sd(rc_scc)), add=TRUE,col=2,lwd=2)

hist(rc_kcc,probability = T,xlim=c(-1,1),main="rc_KCC",cex.main=2,xlab="",cex.axis=1.5);
curve(dnorm(x, mean=mean(rc_kcc), sd=sd(rc_kcc)), add=TRUE,col=2,lwd=2)

hist(rc_gcc,probability = T,xlim=c(-1,1),main="rc_GCC",cex.main=2,xlab="",cex.axis=1.5);
curve(dnorm(x, mean=mean(rc_gcc), sd=sd(rc_gcc)), add=TRUE,col=2,lwd=2)

hist(rc_bicor,probability = T,xlim=c(-1,1),main="rc_Bicor",cex.main=2,xlab="",cex.axis=1.5);
curve(dnorm(x, mean=mean(rc_bicor), sd=sd(rc_bicor)), add=TRUE,col=2,lwd=2)

hist(rc_aa,probability = T,xlim=c(0,2),main="rc_AA",cex.main=2,xlab="",cex.axis=1.5);

hist(rc_ma,probability = T,xlim=c(0,2),main="rc_MA",cex.main=2,xlab="",cex.axis=1.5);

hist(rc_mrnet,probability = T,xlim=c(0,2),main="rc_MRNET",cex.main=2,xlab="",cex.axis=1.5);

hist(rc_clr,probability = T,xlim=c(0,423.8),main="rc_CLR",cex.main=2,xlab="",cex.axis=1.5);

rm(rc_pcc,rc_scc,rc_kcc,rc_gcc,rc_aa,rc_ma,rc_mrnet,rc_clr,rc_bicor)
