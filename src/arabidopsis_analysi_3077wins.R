library(rentrez)
library(vioplot)
# Just search one year of SRA samples.
entrez_search(db="sra",
              term="(Arabidopsis[ORGN] OR Arabidopsis thaliana[ORGN]) AND 2015[PDAT] AND RNA-Seq[STRA]",
              retmax=0)

search_year_sra <- function(year, term){
  query <- paste(term, "AND (", year, "[PDAT])")
  entrez_search(db="sra", term=query, retmax=0)$count
}

year <- 2008:2016
sra_arabidopsis <- sapply(year, search_year_sra, term="(Arabidopsis[ORGN] OR Arabidopsis thaliana[ORGN]) AND RNA-Seq[STRA] AND illumina[PLAT]",
                    USE.NAMES=FALSE)

entrez_search(db="gds",
              term="GPL198[ACCN] AND 2016[PDAT]",
              retmax=0)

# 1022, 1655, 1283, 1241, 1100, 1273,1216, 1015, 802
# from 2008 to 2016



setwd("D:\\Users\\jhuang\\Documents\\Co-expression\\arabidopsis")
library(edgeR)
#raw_data <- read.delim("Round50.featureCount", stringsAsFactors = FALSE,check.names = FALSE)
raw_data  <- read.delim("arabi.txt",stringsAsFactors = FALSE,check.names = FALSE,header=T)
colnames(raw_data)[1] <- "gene_name"

colnames(raw_data)[1364]

filter_data <- raw_data[rowSums(raw_data>2)>=1000,]
filter_data[1:5,1:5]

filter_data[,1] <- substr(filter_data[,1],1,9)
noDupli_data <- filter_data[!duplicated(filter_data$gene_name),]

rownames(noDupli_data) <- noDupli_data[,1]
noDupli_data[,1] <- NULL
noDupli_data[1:5,1:5]

save(noDupli_data,file="arabidopsis_expression.RData")

ara_pcc <- cor(t(noDupli_data),method = "p")
save(ara_pcc,file="arabidopsis_pcc.RData")

library (EGAD)

load("~/Co-expression/arabidopsis/arabidopsis_pcc.RData")
load("~/Co-expression/arabidopsis/arabidopsis_mrnet.RData")
load("~/Co-expression/arabidopsis/arabidopsis_clr.RData")

gotable <- read.table("gene_association.tair_parsed.txt",header=T,sep="\t")
gotable[1:5,1:2]

genes <- unique(gotable$ID)
annotationlist <- unique(gotable$GO)

annotataions_arab <- make_annotations(gotable,genes,annotationlist)
annotations_sub <- filter_network(annotataions_arab,flag = 2,min=20,max = 300)
save(annotations_sub, file="go_annotations_arabSub.RData")

rm(gotable,annotataions_arab)

GO_arab_pcc <- run_GBA(ara_pcc,annotations_sub)
GO_arab_mrnet <- run_GBA(ara_mrnet,annotations_sub)
GO_arab_clr <- run_GBA(ara_clr,annotations_sub)

  GO_arab_pcc[[3]];GO_arab_mrnet[[3]];GO_arab_clr[[3]]

save(GO_arab_pcc,file="arabidopsis_GO_pcc.RData")
save(GO_arab_mrnet,file="arabidopsis_GO_mrnet.RData")
save(GO_arab_clr,file="arabidopsis_GO_clr.RData")

# pairwise wilconxin test

arab1363_GOpcc <- GO_arab_pcc[[1]][,1];arab1363_GOmrnet <- GO_arab_mrnet[[1]][,1];
arab1363_GOclr <- GO_arab_clr[[1]][,1];

arab1363_totalGO <- c(arab1363_GOpcc,arab1363_GOmrnet,arab1363_GOclr)
arab1363_totalGO_name <- c(rep("pcc",445),rep("mrnet",445),rep("clr",445))
pairwise.wilcox.test(arab1363_totalGO,arab1363_totalGO_name,p.adjust.method = 'b')

arab_aggPcc <- GO_arab_agg_pcc[[1]][,1];arab_aggMrnet <- GO_arab_agg_mrnet[[1]][,1]
arab_aggclr <- GO_arab_agg_clr[[1]][,1]
arab_aggTotal <- c(arab_aggPcc,arab_aggMrnet,arab_aggclr);
arab_aggTotal_name <- c(rep('pcc',445),rep('mrnet',445),rep('clr',445))
pairwise.wilcox.test(arab_aggTotal,arab_aggTotal_name,p.adjust.method = 'b')

# Convert PP interaction dataset

convert_list <- function(x){
  return(PP_new[grep(x,PP_new$V1),2])
}

PP_SIF <- read.delim("AtPIN_PPI.txt",header=F,sep="\t")
PP_SIF <- PP_SIF[order(PP_SIF$V1),]

length(which(table(PP_SIF$V1)>=5))
fivenum(table(PP_SIF$V1))
PP_new <- PP_SIF[!(as.numeric(PP_SIF$V1)) %in% which(table(PP_SIF$V1)<5),]

PP_new[,1] <- as.character(PP_new[,1]);PP_new[,2] <- as.character(PP_new[,2])
node_name <- unique(PP_new[,1])
node_number <- length(node_name)
cogene_result <- sapply(node_name,convert_list) # this step takes several minutes

rm(PP_SIF,PP_new,convert_list)

save(cogene_result,file= paste0(node_number,"arab_PP_list.RData"))
save(node_name,file=paste0(node_number,"arab_nodeName.RData"))


library(ROCR)
calc_auc <- function(corMatrix,x){
  t1 <- abs(corMatrix[,x])
  m1 <- data.frame(matrix(nrow = 13519,ncol=3))
  m1[,1] <- as.character();m1[,2] <- as.double();m1[,3] <- 0
  m1[,1] <- names(t1);m1[,2] <- t1;m1[which(m1[,1] %in% cogene_result[[x]]),3] <- 1
  if (length(table(m1[,3]))==2){
    pred <- prediction(m1[,2],m1[,3])
    auc <- performance(pred,"auc");auc <- unlist(slot(auc, "y.values"))
  } else {
    auc <- NULL
  }
  return(auc)
}

filter_node_name <- node_name[which(node_name %in% colnames(ara_mrnet))]
cogene_result <- cogene_result[names(cogene_result) %in% filter_node_name]

ara_pcc_new <- (ara_pcc - mean(ara_pcc))/sd(ara_pcc) 
# the pcc correlaition matrix was right shifted a little.
# this cause AUROC higher, like what happend for RC.So I did a normalization.
ara_pcc_ft <- ara_pcc_new[,which(colnames(ara_pcc_new) %in% filter_node_name)] # dim vst_gcc_ft 15116*753                        #
auc_ara_pcc <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=ara_pcc_new);
mean(unlist(auc_ara_pcc));

ara_mrnet_ft <- ara_mrnet[,which(colnames(ara_mrnet) %in% filter_node_name)] # dim vst_gcc_ft 15116*753                        #
auc_ara_mrnet <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=ara_mrnet);
mean(unlist(auc_ara_mrnet));

ara_clr_ft <- ara_clr[,which(colnames(ara_clr) %in% filter_node_name)] # dim vst_gcc_ft 15116*753                        #
auc_ara_clr <- sapply(X = filter_node_name,FUN = calc_auc,corMatrix=ara_clr);
mean(unlist(auc_ara_clr));

save(auc_ara_pcc,auc_ara_mrnet,auc_ara_clr,file="auc_PP_3Methods.RData")

ap <- unlist(auc_ara_pcc);am <- unlist(auc_ara_mrnet);ac <- unlist(auc_ara_clr)


library(vioplot)

al <- c(ap,am,ac);aname <- c(rep("pcc",3131),rep("mrnet",3131),rep("clr",3131))
ap <- data.frame(cbind(al,aname),row.names = c(1:9393))
ap$al <- as.numeric(ap$al)
vioplot(ap, am, ac, names=c("pcc", "mrnet", "clr"), col="gold")
# this test confirmed mrnet and clr are better than pcc.
pairwise.wilcox.test(ap$al,ap$aname,p.adjust.method = "b")





