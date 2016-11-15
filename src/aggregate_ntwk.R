# This script is to code different size of network and only calculate PCC, MRNET and CLR.
library(parmigene)
library(edgeR)

setwd("D:\\Users\\jhuang\\Documents\\Co-expression")
setwd("/home/bio.local/jhuang/data/aggregate_ntwk")
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
  
  save(list=paste0("cpm",output,"_pcc"),file=paste0(output,"_cpm_pcc.RData"))
  save(list=paste0("cpm",output,"_MI"),file=paste0(output,"_cpm_MI.RData"))
  save(list=paste0("cpm",output,"_mrnet"),file=paste0(output,"_cpm_mrnet.RData"))
  save(list=paste0("cpm",output,"_clr"),file=paste0(output,"_cpm_clr.RData"))
}


aggre_exp <- read.delim("12.txt",header = F)
aggregate_cpm(aggre_exp,"12")

ptm <- proc.time()
aggre_exp <- read.delim("36.txt",header = F)
aggregate_cpm(aggre_exp,"36")
proc.time()  - ptm # 20012 s


ptm <- proc.time()
aggre_exp <- read.delim("65.txt",header = F)
aggregate_cpm(aggre_exp,"65")
proc.time()  - ptm # 22060 s

ptm <- proc.time()
aggre_exp <- read.delim("108.txt",header = F)
aggregate_cpm(aggre_exp,"108")
proc.time()  - ptm # 25620 s

ptm <- proc.time()
aggre_exp <- read.delim("270.txt",header = F)
aggregate_cpm(aggre_exp,"270")
proc.time()  - ptm # 37187 s

ptm <- proc.time()
aggre_exp <- read.delim("404.txt",header = F)
aggregate_cpm(aggre_exp,"404")
proc.time() - ptm # 52316 s

ptm <- proc.time()
aggre_exp <- read.delim("exp7.txt",header = F)
aggregate_cpm(aggre_exp,"exp7")
proc.time()  - ptm #4038s
rm(cpmexp7_pcc,cpmexp7_clr,cpmexp7_mrnet)

ptm <- proc.time()
aggre_exp <- read.delim("exp8.txt",header = F)
aggregate_cpm(aggre_exp,"exp8")
proc.time()  - ptm #3963s
rm(cpmexp8_pcc,cpmexp8_clr,cpmexp8_mrnet)

ptm <- proc.time()
aggre_exp <- read.delim("exp9.txt",header = F)
aggregate_cpm(aggre_exp,"exp9")
proc.time()  - ptm #3184s
rm(cpmexp9_pcc,cpmexp9_clr,cpmexp9_mrnet)

ptm <- proc.time()
aggre_exp <- read.delim("exp10.txt",header = F)
aggregate_cpm(aggre_exp,"exp10")
proc.time()  - ptm #3737s
rm(cpmexp10_pcc,cpmexp10_clr,cpmexp10_mrnet)

ptm <- proc.time()
aggre_exp <- read.delim("exp11.txt",header = F)
aggregate_cpm(aggre_exp,"exp11")
proc.time()  - ptm #3258s
rm(cpmexp11_pcc,cpmexp11_clr,cpmexp11_mrnet)

ptm <- proc.time()
aggre_exp <- read.delim("exp12.txt",header = F)
aggregate_cpm(aggre_exp,"exp12")
proc.time()  - ptm # 3061s
rm(cpmexp11_pcc,cpmexp11_clr,cpmexp11_mrnet)

ptm <- proc.time()
aggre_exp <- read.delim("exp13.txt",header = F)
aggregate_cpm(aggre_exp,"exp13")
proc.time()  - ptm #3023s
rm(cpmexp13_pcc,cpmexp13_clr,cpmexp13_mrnet)


ptm <- proc.time()
aggre_exp <- read.delim("exp14.txt",header = F)
aggregate_cpm(aggre_exp,"exp14")
proc.time()  - ptm #3129s
rm(cpmexp14_pcc,cpmexp14_clr,cpmexp14_mrnet)

ptm <- proc.time()
aggre_exp <- read.delim("exp15.txt",header = F)
aggregate_cpm(aggre_exp,"exp15")
proc.time()  - ptm #2890s
rm(cpmexp15_pcc,cpmexp15_clr,cpmexp15_mrnet)