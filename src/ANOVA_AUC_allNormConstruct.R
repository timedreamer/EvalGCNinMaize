#KIN3077 Windows
load("~/Co-expression/Normalization result/Before_Oct2016/AUC/GO_cpm.RData")
load("~/Co-expression/Normalization result/Before_Oct2016/AUC/GO_vst.RData")
load("~/Co-expression/Normalization result/Before_Oct2016/AUC/GO_rc.RData")
load("~/Co-expression/Normalization result/Before_Oct2016/AUC/GO_rpkm.RData")

load("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP/total_ntwk/AUC_vst.RData")
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP/total_ntwk/AUC_rc.RData")
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP/total_ntwk/AUC_cpm.RData")
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/AUC_PTY_PP/total_ntwk/AUC_rpkm.RData")


####################
# For PPPTY


##build dataframe with three columns. First, AUROC value. Second, Construction method
# Third, Normalization method.

#VST
vst_pcc <- unlist(unname(auc_vst_pcc));names(vst_pcc) <- rep("pcc",1721)
vst_scc <- unlist(unname(auc_vst_scc));names(vst_scc) <- rep("scc",1721)
vst_kcc <- unlist(unname(auc_vst_kcc));names(vst_kcc) <- rep("kcc",1721)
vst_gcc <- unlist(unname(auc_vst_gcc));names(vst_gcc) <- rep("gcc",1721)
vst_bicor <- unlist(unname(auc_vst_bicor));names(vst_bicor) <- rep("bicor",1721)

vst_aa <- unlist(unname(auc_vst_aa));names(vst_aa) <- rep("aa",1721)
vst_ma <- unlist(unname(auc_vst_ma));names(vst_ma) <- rep("ma",1721)
vst_mrnet <- unlist(unname(auc_vst_mrnet));names(vst_mrnet) <- rep("mrnet",1721)
vst_clr <- unlist(unname(auc_vst_clr));names(vst_clr) <- rep("clr",1721)

total_vst <- c(vst_pcc,vst_scc,vst_kcc,vst_gcc,vst_bicor,vst_aa,vst_ma,vst_mrnet,vst_clr)
total_name_vst <- c(names(vst_pcc),names(vst_scc),names(vst_kcc),names(vst_gcc),names(vst_bicor),
               names(vst_aa),names(vst_ma),names(vst_mrnet),names(vst_clr))
norm_vst <- rep("vst",15489)


#CPM
cpm_pcc <- unlist(unname(auc_cpm_pcc));names(cpm_pcc) <- rep("pcc",1721)
cpm_scc <- unlist(unname(auc_cpm_scc));names(cpm_scc) <- rep("scc",1721)
cpm_kcc <- unlist(unname(auc_cpm_kcc));names(cpm_kcc) <- rep("kcc",1721)
cpm_gcc <- unlist(unname(auc_cpm_gcc));names(cpm_gcc) <- rep("gcc",1721)
cpm_bicor <- unlist(unname(auc_cpm_bicor));names(cpm_bicor) <- rep("bicor",1721)

cpm_aa <- unlist(unname(auc_cpm_aa));names(cpm_aa) <- rep("aa",1721)
cpm_ma <- unlist(unname(auc_cpm_ma));names(cpm_ma) <- rep("ma",1721)
cpm_mrnet <- unlist(unname(auc_cpm_mrnet));names(cpm_mrnet) <- rep("mrnet",1721)
cpm_clr <- unlist(unname(auc_cpm_clr));names(cpm_clr) <- rep("clr",1721)


total_cpm <- c(cpm_pcc,cpm_scc,cpm_kcc,cpm_gcc,cpm_bicor,cpm_aa,cpm_ma,cpm_mrnet,cpm_clr)
total_name_cpm <- c(names(cpm_pcc),names(cpm_scc),names(cpm_kcc),names(cpm_gcc),names(cpm_bicor),
                names(cpm_aa),names(cpm_ma),names(cpm_mrnet),names(cpm_clr))
norm_cpm <- rep("cpm",15489)


#RC
rc_pcc <- unlist(unname(auc_rc_pcc));names(rc_pcc) <- rep("pcc",1721)
rc_scc <- unlist(unname(auc_rc_scc));names(rc_scc) <- rep("scc",1721)
rc_kcc <- unlist(unname(auc_rc_kcc));names(rc_kcc) <- rep("kcc",1721)
rc_gcc <- unlist(unname(auc_rc_gcc));names(rc_gcc) <- rep("gcc",1721)
rc_bicor <- unlist(unname(auc_rc_bicor));names(rc_bicor) <- rep("bicor",1721)

rc_aa <- unlist(unname(auc_rc_aa));names(rc_aa) <- rep("aa",1721)
rc_ma <- unlist(unname(auc_rc_ma));names(rc_ma) <- rep("ma",1721)
rc_mrnet <- unlist(unname(auc_rc_mrnet));names(rc_mrnet) <- rep("mrnet",1721)
rc_clr <- unlist(unname(auc_rc_clr));names(rc_clr) <- rep("clr",1721)


total_rc <- c(rc_pcc,rc_scc,rc_kcc,rc_gcc,rc_bicor,rc_aa,rc_ma,rc_mrnet,rc_clr)
total_name_rc <- c(names(rc_pcc),names(rc_scc),names(rc_kcc),names(rc_gcc),names(rc_bicor),
                names(rc_aa),names(rc_ma),names(rc_mrnet),names(rc_clr))
norm_rc <- rep("rc",15489)


#RPKM
rpkm_pcc <- unlist(unname(auc_rpkm_pcc));names(rpkm_pcc) <- rep("pcc",1721)
rpkm_scc <- unlist(unname(auc_rpkm_scc));names(rpkm_scc) <- rep("scc",1721)
rpkm_kcc <- unlist(unname(auc_rpkm_kcc));names(rpkm_kcc) <- rep("kcc",1721)
rpkm_gcc <- unlist(unname(auc_rpkm_gcc));names(rpkm_gcc) <- rep("gcc",1721)
rpkm_bicor <- unlist(unname(auc_rpkm_bicor));names(rpkm_bicor) <- rep("bicor",1721)

rpkm_aa <- unlist(unname(auc_rpkm_aa));names(rpkm_aa) <- rep("aa",1721)
rpkm_ma <- unlist(unname(auc_rpkm_ma));names(rpkm_ma) <- rep("ma",1721)
rpkm_mrnet <- unlist(unname(auc_rpkm_mrnet));names(rpkm_mrnet) <- rep("mrnet",1721)
rpkm_clr <- unlist(unname(auc_rpkm_clr));names(rpkm_clr) <- rep("clr",1721)


total_rpkm <- c(rpkm_pcc,rpkm_scc,rpkm_kcc,rpkm_gcc,rpkm_bicor,rpkm_aa,rpkm_ma,rpkm_mrnet,rpkm_clr)
total_name_rpkm <- c(names(rpkm_pcc),names(rpkm_scc),names(rpkm_kcc),names(rpkm_gcc),names(rpkm_bicor),
                names(rpkm_aa),names(rpkm_ma),names(rpkm_mrnet),names(rpkm_clr))
norm_rpkm <- rep("rpkm",15489)


# Add together. with RC
#total_value <- c(total_vst,total_cpm,total_rc,total_rpkm)
#total_construct <- c(total_name_vst,total_name_cpm,total_name_rc,total_name_rpkm)

# Add together without RC.
total_value <- c(total_vst,total_cpm,total_rpkm)
total_construct <- c(total_name_vst,total_name_cpm,total_name_rpkm)
total_normMethod <- c(norm_vst,norm_cpm,norm_rpkm)

v2 <- data.frame(total_value,total_construct,total_normMethod)


# Two-way anova and pairwise comparison
v2$total_normMethod <- factor(v2$total_normMethod,levels = c("vst","cpm","rpkm"))
v2$total_construct <- factor(v2$total_construct,
                          level = c("pcc","scc","kcc","gcc","bicor","aa","ma",
                                    "mrnet","clr"))
# plot is the same as plot singe var.
par(mfrow=c(1,2)) # export pdf 15*9
plot(total_value ~ total_construct*total_normMethod, data=v2,ylab = "AUC",xlab="var")

######

a2 <- aov(total_value ~ total_construct*total_normMethod,data=v2)
res<-a2$residuals
hist(res,main="Histogram of residuals",xlab="Residuals")
summary(a2)
par(mar=c(5,5,4,2))
# this plot is the same as what I drew for the mean. redundent.
interaction.plot(total_netname,total_constructMethod,total_value,type="b",
                 col=c("red","green","blue"),lwd=2,bty="n",cex.axis=1.2,
                 cex.lab=1.2) 
box(lwd=2)
tky2 <- TukeyHSD(a2)
tky2

# # unused code. Use tukeyHSD instead.
# anova(lm(total_value ~ total_construct*total_normMethod,data=v2)) # show anova table directly
# 
# pairwise.t.test(total_value, total_construct, p.adjust="bonferroni")
# pairwise.t.test(total_value, total_normMethod, p.adjust="bonferroni")
# 
# pairwise.wilcox.test(total_value, total_construct, p.adjust="bonferroni")
# pairwise.wilcox.test(total_value, total_normMethod, p.adjust="bonferroni")


# # Onw-way anova, did not use
# v1_vst <- data.frame(total_vst,total_name)
# 
# plot(total_vst ~ total_name, data=v1_vst)
# v1_result <-  aov(total_vst ~ total_name, data=v1_vst)
# summary(v1_result) # diff
# 
# pairwise.t.test(total_vst, total_name, p.adjust="bonferroni")
# TukeyHSD(v1_result, conf.level = 0.95)
##########################################


###########For Gene Ontology
vst_pcc <- GO_vst_pcc[[1]][,1];vst_scc <- GO_vst_scc[[1]][,1]
vst_kcc <- GO_vst_kcc[[1]][,1];vst_gcc <- GO_vst_gcc[[1]][,1]
vst_bicor <- GO_vst_bicor[[1]][,1]
vst_aa <- GO_vst_aa[[1]][,1];vst_ma <- GO_vst_ma[[1]][,1]
vst_mrnet <- GO_vst_mrnet[[1]][,1];vst_clr <- GO_vst_clr[[1]][,1]

total_vst <- c(vst_pcc,vst_scc,vst_kcc,vst_gcc,vst_bicor,vst_aa,vst_ma,vst_mrnet,vst_clr)
total_name_vst <- c(rep("pcc",277),rep("scc",277),rep("kcc",277),rep("gcc",277),rep("bicor",277),
                    rep("aa",277),rep("ma",277),rep("mrnet",277),rep("clr",277))
norm_vst <- rep("vst",2493)


cpm_pcc <- GO_cpm_pcc[[1]][,1];cpm_scc <- GO_cpm_scc[[1]][,1]
cpm_kcc <- GO_cpm_kcc[[1]][,1];cpm_gcc <- GO_cpm_gcc[[1]][,1]
cpm_bicor <- GO_cpm_bicor[[1]][,1]
cpm_aa <- GO_cpm_aa[[1]][,1];cpm_ma <- GO_cpm_ma[[1]][,1]
cpm_mrnet <- GO_cpm_mrnet[[1]][,1];cpm_clr <- GO_cpm_clr[[1]][,1]

total_cpm <- c(cpm_pcc,cpm_scc,cpm_kcc,cpm_gcc,cpm_bicor,cpm_aa,cpm_ma,cpm_mrnet,cpm_clr)
total_name_cpm <- c(rep("pcc",277),rep("scc",277),rep("kcc",277),rep("gcc",277),rep("bicor",277),
                    rep("aa",277),rep("ma",277),rep("mrnet",277),rep("clr",277))
norm_cpm <- rep("cpm",2493)


rc_pcc <- GO_rc_pcc[[1]][,1];rc_scc <- GO_rc_scc[[1]][,1]
rc_kcc <- GO_rc_kcc[[1]][,1];rc_gcc <- GO_rc_gcc[[1]][,1]
rc_bicor <- GO_rc_bicor[[1]][,1]
rc_aa <- GO_rc_aa[[1]][,1];rc_ma <- GO_rc_ma[[1]][,1]
rc_mrnet <- GO_rc_mrnet[[1]][,1];rc_clr <- GO_rc_clr[[1]][,1]

total_rc <- c(rc_pcc,rc_scc,rc_kcc,rc_gcc,rc_bicor,rc_aa,rc_ma,rc_mrnet,rc_clr)
total_name_rc <-c(rep("pcc",277),rep("scc",277),rep("kcc",277),rep("gcc",277),rep("bicor",277),
                  rep("aa",277),rep("ma",277),rep("mrnet",277),rep("clr",277))
norm_rc <- rep("rc",2493)


rpkm_pcc <- GO_rpkm_pcc[[1]][,1];rpkm_scc <- GO_rpkm_scc[[1]][,1]
rpkm_kcc <- GO_rpkm_kcc[[1]][,1];rpkm_gcc <- GO_rpkm_gcc[[1]][,1]
rpkm_bicor <- GO_rpkm_bicor[[1]][,1]
rpkm_aa <- GO_rpkm_aa[[1]][,1];rpkm_ma <- GO_rpkm_ma[[1]][,1]
rpkm_mrnet <- GO_rpkm_mrnet[[1]][,1];rpkm_clr <- GO_rpkm_clr[[1]][,1]

total_rpkm <- c(rpkm_pcc,rpkm_scc,rpkm_kcc,rpkm_gcc,rpkm_bicor,rpkm_aa,rpkm_ma,rpkm_mrnet,rpkm_clr)
total_name_rpkm <- c(rep("pcc",277),rep("scc",277),rep("kcc",277),rep("gcc",277),rep("bicor",277),
                     rep("aa",277),rep("ma",277),rep("mrnet",277),rep("clr",277))
norm_rpkm <- rep("rpkm",2493)


# Add together. wth RC
# total_value <- c(total_vst,total_cpm,total_rc,total_rpkm)
# total_construct <- c(total_name_vst,total_name_cpm,total_name_rc,total_name_rpkm)
# total_normMethod <- c(norm_vst,norm_cpm,norm_rc,norm_rpkm)

# Add together without RC.
total_value <- c(total_vst,total_cpm,total_rpkm)
total_construct <- c(total_name_vst,total_name_cpm,total_name_rpkm)
total_normMethod <- c(norm_vst,norm_cpm,norm_rpkm)


v3 <- data.frame(total_value,total_construct,total_normMethod)


a3 <- aov(total_value ~ total_construct*total_normMethod,data=v3)
res<-a3$residuals
hist(res,main="Histogram of residuals",xlab="Residuals")
summary(a3)
par(mar=c(5,5,4,2))
# this plot is the same as what I drew for the mean. redundent.
interaction.plot(total_netname,total_constructMethod,total_value,type="b",
                 col=c("red","green","blue"),lwd=2,bty="n",cex.axis=1.2,
                 cex.lab=1.2) 
box(lwd=2)
tky3 <- TukeyHSD(a3)
tky3


# # Two-way anova and pairwise comparison
# v3$total_normMethod <- factor(v3$total_normMethod,levels = c("vst","rc","cpm","rpkm"))
# v3$total_construct <- factor(v3$total_construct,
#                              level = c("pcc","scc","kcc","gcc","bicor","aa","ma",
#                                        "mrnet","clr"))
# # plot is the same as plot singe var.
# par(mfrow=c(1,2)) # export pdf 15*9
# plot(total_value ~ total_construct*total_normMethod, data=v3,ylab = "AUC",xlab="var")
# 
# anova(lm(total_value ~ total_construct*total_normMethod,data=v3)) # show anova table directly
# 
# pairwise.t.test(total_value, total_construct, p.adjust="bonferroni")
# pairwise.t.test(total_value, total_normMethod, p.adjust="bonferroni")
# 
# pairwise.wilcox.test(total_value, total_construct, p.adjust="bonferroni")
# pairwise.wilcox.test(total_value, total_normMethod, p.adjust="bonferroni")
# 
