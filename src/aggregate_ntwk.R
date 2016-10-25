
library(parmigene)
library(edgeR)

setwd("C:/WORK/Co-expression/")
setwd("D:\\Users\\jhuang\\Documents\\Co-expression")
load("15116_geneNames.RData")
#
raw_data  <- read.delim("ALL_FC_noDuplicateLib_biggerThan5Million_70allignmentRate_1266.txt",
stringsAsFactors = FALSE,check.names = FALSE)



aggregate_cpm <- function(x,output) {

  x <- t(x)
  raw_data_agg <- raw_data[,which(colnames(raw_data) %in% x)]
  raw_data_agg <- cbind(raw_data[,c(1,2)],raw_data_agg)
  raw_data_agg <- raw_data_agg[which(raw_data_agg[,1] %in% ntwk_geneName),]
  gene_name <- unname(raw_data_agg$Geneid)
  y <- DGEList(counts = raw_data_agg[3:ncol(raw_data_agg)],genes = gene_name)
  y <- calcNormFactors(y)
  counts.p.m <- cpm(y,normalized.lib.sizes = TRUE,log = TRUE,prior.count = 1)
  rownames(counts.p.m) <- gene_name
  assign(paste0("cpm",output,"_pcc"),cor(t(counts.p.m),method = "pearson"))
  assign(paste0("cpm",output,"_MI"),knnmi.all(counts.p.m,k=3))
  assign(paste0("cpm",output,"_mrnet"),mrnet(eval(as.name(paste0("cpm",output,"_MI")))))
  assign(paste0("cpm",output,"_clr"),clr(eval(as.name(paste0("cpm",output,"_MI")))))
  #cpm_pcc <- cor(t(counts.p.m),method = "pearson")
  # cpm_MI <- knnmi.all(counts.p.m,k=3)
  # cpm_mrnet <- mrnet(cpm_MI)
  # cpm_clr <- clr(cpm_MI)
  save(list=paste0("cpm",output,"_pcc"),file=paste0(output,"_cpm_pcc.RData"))
  save(list=paste0("cpm",output,"_MI"),file=paste0(output,"_cpm_MI.RData"))
  save(list=paste0("cpm",output,"_mrnet"),file=paste0(output,"_cpm_mrnet.RData"))
  save(list=paste0("cpm",output,"_clr"),file=paste0(output,"_cpm_clr.RData"))
  # save(cpm_MI,file=paste0(output,"_cpm_MI.RData"))
  # save(cpm_mrnet,file=paste0(output,"_cpm_mrnet.RData"))
  # save(cpm_clr,file=paste0(output,"_cpm_clr.RData"))
}


aggre_exp <- read.delim("12.txt",header = F)
aggregate_cpm(aggre_exp,"12")
