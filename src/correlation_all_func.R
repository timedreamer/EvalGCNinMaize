# This is the file for calcuate diffrent correlation method.

#######################################################################################
# correlation matrix builder for different normalization method.
# PCC,SCC are built using cor()
# GCC use rsgcc
# KCC use pcapp
# bicor use WGCNA
# AA, MA, MRNET and CLR use parmigene.

# adjacencymatrix needs rowname and colname in order to calculate.

library(rsgcc)
library(pcaPP)
library(parmigene)
library(WGCNA)



calc_corr_all <- function(x){
  if (deparse(substitute(x)) == "vst_norm"){
    vst_pcc <- cor(t(x),method = "pearson") # 224s
    vst_scc <- cor(t(x),method="s") # 228s
    save(vst_scc,vst_pcc,file="vst_pcc_scc.RData");rm(vst_scc,vst_pcc)
    vst_gcc <- adjacencymatrix(x,method = "GCC",cpus =12) # 2hours
    vst_kcc <- cor.fk(t(x))
    save(vst_gcc,vst_kcc,file="vst_gcc_kcc.RData");rm(vst_gcc,vst_kcc)
    vst_bicor <- bicor(t(x),nThreads = 12)
    save(vst_bicor,file="vst_bicor.RData");rm(vst_bicor)
    vst_MI <- knnmi.all(x,k=3)
    vst_aa <-aracne.a(vst_MI);vst_ma <- aracne.m(vst_MI)
    vst_mrnet <- mrnet(vst_MI);vst_clr <- clr(vst_MI)
    save(vst_aa,vst_ma,vst_mrnet,vst_clr,file="vst_four_MI.RData")
    save(vst_MI,file="vst_MI_raw.RData");
    rm(vst_aa,vst_ma,vst_mrnet,vst_clr,vst_MI)
  } else if (deparse(substitute(x)) == "cpm_norm"){
    cpm_pcc <- cor(t(x),method = "pearson") # 224s
    cpm_scc <- cor(t(x),method="s") # 228s
    save(cpm_scc,cpm_pcc,file="cpm_pcc_scc.RData");rm(cpm_scc,cpm_pcc)
    cpm_gcc <- adjacencymatrix(x,method = "GCC",cpus =12) # 2hours
    cpm_kcc <- cor.fk(t(x))
    save(cpm_gcc,cpm_kcc,file="cpm_gcc_kcc.RData");rm(cpm_gcc,cpm_kcc)
    cpm_bicor <- bicor(t(x),nThreads = 12)
    save(cpm_bicor,file="cpm_bicor.RData");rm(cpm_bicor)
    cpm_MI <- knnmi.all(x,k=3)
    cpm_aa <-aracne.a(cpm_MI);cpm_ma <- aracne.m(cpm_MI)
    cpm_mrnet <- mrnet(cpm_MI);cpm_clr <- clr(cpm_MI)
    save(cpm_aa,cpm_ma,cpm_mrnet,cpm_clr,file="cpm_four_MI.RData")
    save(cpm_MI,file="cpm_MI_raw.RData");
    rm(cpm_aa,cpm_ma,cpm_mrnet,cpm_clr,cpm_MI)
  } else if (deparse(substitute(x)) =="rc_norm"){
    rc_pcc <- cor(t(x),method = "pearson") # 224s
    rc_scc <- cor(t(x),method="s") # 228s
    save(rc_scc,rc_pcc,file="rc_pcc_scc.RData");rm(rc_scc,rc_pcc)
    rc_gcc <- adjacencymatrix(x,method = "GCC",cpus =12) # 2hours
    rc_kcc <- cor.fk(t(x))
    save(rc_gcc,rc_kcc,file="rc_gcc_kcc.RData");rm(rc_gcc,rc_kcc)
    rc_bicor <- bicor(t(x),nThreads = 12)
    save(rc_bicor,file="rc_bicor.RData");rm(rc_bicor)
    rc_MI <- knnmi.all(x,k=3)
    rc_aa <-aracne.a(rc_MI);rc_ma <- aracne.m(rc_MI)
    rc_mrnet <- mrnet(rc_MI);rc_clr <- clr(rc_MI)
    save(rc_aa,rc_ma,rc_mrnet,rc_clr,file="rc_four_MI.RData")
    save(rc_MI,file="rc_MI_raw.RData");
    rm(rc_aa,rc_ma,rc_mrnet,rc_clr,rc_MI)
  } else if (deparse(substitute(x)) =="rpkm_norm") {
    rpkm_pcc <- cor(t(x),method = "pearson") # 224s
    rpkm_scc <- cor(t(x),method="s") # 228s
    save(rpkm_scc,rpkm_pcc,file="rpkm_pcc_scc.RData");rm(rpkm_scc,rpkm_pcc)
    rpkm_gcc <- adjacencymatrix(x,method = "GCC",cpus =12) # 2hours
    rpkm_kcc <- cor.fk(t(x))
    save(rpkm_gcc,rpkm_kcc,file="rpkm_gcc_kcc.RData");rm(rpkm_gcc,rpkm_kcc)
    rpkm_bicor <- bicor(t(x),nThreads = 12)
    save(rpkm_bicor,file="rpkm_bicor.RData");rm(rpkm_bicor)
    rpkm_MI <- knnmi.all(x,k=3)
    rpkm_aa <-aracne.a(rpkm_MI);rpkm_ma <- aracne.m(rpkm_MI)
    rpkm_mrnet <- mrnet(rpkm_MI);rpkm_clr <- clr(rpkm_MI)
    save(rpkm_aa,rpkm_ma,rpkm_mrnet,rpkm_clr,file="rpkm_four_MI.RData")
    save(rpkm_MI,file="rpkm_MI_raw.RData");
    rm(rpkm_aa,rpkm_ma,rpkm_mrnet,rpkm_clr,rpkm_MI)
  }
}

#read data
setwd("D:\\Users\\jhuang\\Documents\\Co-expression\\Normalization result/")

###################Calculation########################################
vst_norm <- read.delim("VST_result.txt",sep="\t",header=TRUE) # no geneLength column
calc_corr_all(vst_norm)

rc_norm <- read.delim("rawCount_log2_result.txt",sep="\t",header=TRUE) # no geneLength column
calc_corr_all(rc_norm)

cpm_norm <- read.delim("CPM_log2.txt",sep="\t",header=T) # geneLength last column
cpm_norm[1267] <- NULL # need to delete last column
calc_corr_all(cpm_norm)

rpkm_norm <- read.delim("RPKM_log2.txt",sep="\t",header=T) # geneLength last column
rpkm_norm[1267] <- NULL
calc_corr_all(rpkm_norm)


# # calculate ED, using rsgcc
# ptm <- proc.time()
# cpm_ed <- adjacencymatrix(cpm_norm,method = "ED",cpus = 12)
# proc.time() - ptm
# save(cpm_ed,file="cpm_ed.RData")
#
#
# ptm <- proc.time()
# rpkm_ed <- adjacencymatrix(rpkm_ed,method = "ED",cpus = 12)
# proc.time() - ptm
# save(rpkm_ed,file="rpkm_ed.RData")
