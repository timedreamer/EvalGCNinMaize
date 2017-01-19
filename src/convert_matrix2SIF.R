# This file is to cnovert correlation matrix file to SIF file.
# For correlation, change diag to 0 to avoid self-connections.
# The resulting SIF can be imported in Cytoscape.

#Maize Final Network
# KIN3077 Windows
library(igraph)
load("~/Co-expression/aggregate_ntwk/MaizeRank/pcc/maizeAgg_Rankpcc.RData") # PCC use agg15 Rank
# MRNET and CLR use 1266 one matrix data
load("~/Co-expression/Normalization result/FromNewScript_Oct2016/cpm_20161024/cpm_four_MI.RData")
rm(cpm_aa,cpm_ma)


# PCC_rank_15experiments
diag(agg_RankntwkPcc) <- 0
#choose threshold.top 1 Million
one_million <- 2e06/(15116*15116)
g  <- graph.adjacency((agg_RankntwkPcc+t(agg_RankntwkPcc))/2,mode = "undirected",weighted=TRUE)
df <- get.data.frame(g)
head(df)
threshold <- quantile(df[,3],1-one_million)
df_pass1 <- subset(df,weight>threshold)
write.table(df_pass1,file = paste0("agg_RankntwkPcc_",round(threshold,3),".SIF"),quote = F,sep = "\t",col.names = T,row.names = F)


# MRNET 1266 maize
diag(cpm_mrnet) <- 0
#choose threshold.top 1 Million
one_million <- 1e06/(73689151)
g  <- graph.adjacency((cpm_mrnet+t(cpm_mrnet))/2,mode = "undirected",weighted=TRUE)
df <- get.data.frame(g)
head(df)
threshold <- quantile(df[,3],1-one_million)
df_pass1 <- subset(df,weight>threshold)
write.table(df_pass1,file = paste0("cpm_mrnet_",round(threshold,3),".SIF"),quote = F,sep = "\t",col.names = T,row.names = F)


# CLR_1266_maize.
diag(cpm_clr) <- 0
#choose threshold.top 1 Million
one_million <- 1e06/(73201950)
g  <- graph.adjacency((cpm_clr+t(cpm_clr))/2,mode = "undirected",weighted=TRUE)
df <- get.data.frame(g)
head(df)
threshold <- quantile(df[,3],1-one_million)
df_pass1 <- subset(df,weight>threshold)
write.table(df_pass1,file = paste0("cpm_clr_",round(threshold,3),".SIF"),quote = F,sep = "\t",col.names = T,row.names = F)


#############################################
##Arabidopsis Final Network
############################################
load("~/Co-expression/finalized_networks/Arabidopsis/arab_agg_ntwkRank_pcc.Rdata")
load("~/Co-expression/finalized_networks/Arabidopsis/arabidopsis_clr.RData")
load("~/Co-expression/finalized_networks/Arabidopsis/arabidopsis_mrnet.RData")


# PCC_rank_44experiments
diag(agg_ntwk_pccRank) <- 0
#choose threshold.top 1 Million
g  <- graph.adjacency((agg_ntwk_pccRank+t(agg_ntwk_pccRank))/2,mode = "undirected",weighted=TRUE)
df <- get.data.frame(g)
head(df)
one_million <- 1e06/length(df[,3])
threshold <- quantile(df[,3],1-one_million)
df_pass1 <- subset(df,weight>threshold)
write.table(df_pass1,file = paste0("Arab_agg_ntwk_pccRank_",round(threshold,3),".SIF"),quote = F,sep = "\t",col.names = T,row.names = F)


# MRNET 1363 single network
diag(ara_mrnet) <- 0
#choose threshold.top 1 Million
g  <- graph.adjacency((ara_mrnet+t(ara_mrnet))/2,mode = "undirected",weighted=TRUE)
df <- get.data.frame(g)
head(df)
one_million <- 1e06/length(df[,3])
threshold <- quantile(df[,3],1-one_million)
df_pass1 <- subset(df,weight>threshold)
write.table(df_pass1,file = paste0("Arab_ara_mrnet_",round(threshold,3),".SIF"),quote = F,sep = "\t",col.names = T,row.names = F)

# CLR 1363 single network
diag(ara_clr) <- 0
#choose threshold.top 1 Million
g  <- graph.adjacency((ara_clr+t(ara_clr))/2,mode = "undirected",weighted=TRUE)
df <- get.data.frame(g)
head(df)
one_million <- 1e06/length(df[,3])
threshold <- quantile(df[,3],1-one_million)
df_pass1 <- subset(df,weight>threshold)
write.table(df_pass1,file = paste0("Arab_ara_clr_",round(threshold,3),".SIF"),quote = F,sep = "\t",col.names = T,row.names = F)
