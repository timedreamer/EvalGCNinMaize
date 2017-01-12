# This file is to cnovert correlation matrix file to SIF file.
# For correlation, change diag to 0 to avoid self-connections.
# The resulting SIF can be imported in Cytoscape.

# KIN3077 Windows
library(igraph)
load("~/Co-expression/aggregate_ntwk/agg15_ntwkPCC.Rdata")
load("~/Co-expression/aggregate_ntwk/agg15_ntwkMrnet.Rdata")
load("~/Co-expression/aggregate_ntwk/agg15_ntwkClr.Rdata")

# CPM_all_fifteen_maize
# use absolute value
agg15_ntwkPcc <- abs(agg15_ntwkPcc)
diag(agg15_ntwkPcc) <- 0

#choose threshold.top 1 Million 
one_million <- 1e06/(15116*15116)

threshold <- quantile(agg15_ntwkPcc,1-one_million)

ptm <- proc.time()
g  <- graph.adjacency((agg15_ntwkPcc+t(agg15_ntwkPcc))/2,mode = "undirected",weighted=TRUE)
df <- get.data.frame(g)
head(df)
df_pass1 <- subset(df,weight>threshold)
write.table(df_pass1,file = paste0("agg15_ntwkPcc_",round(threshold,3),".SIF"),quote = F,sep = "\t",col.names = T,row.names = F)
proc.time() - ptm # For threshold top 1million, 80.73s


# MRNET_all_fifteen_maize
# use absolute value
agg15_ntwkMrnet <- abs(agg15_ntwkMrnet)
diag(agg15_ntwkMrnet) <- 0

#choose threshold.top 1 Million 
one_million <- 1e06/(15116*15116)

threshold <- quantile(agg15_ntwkMrnet,1-one_million)

ptm <- proc.time()
g  <- graph.adjacency((agg15_ntwkMrnet+t(agg15_ntwkMrnet))/2,mode = "undirected",weighted=TRUE)
df <- get.data.frame(g)
head(df)
df_pass1 <- subset(df,weight>threshold)
write.table(df_pass1,file = paste0("agg15_ntwkMrnet_",round(threshold,3),".SIF"),quote = F,sep = "\t",col.names = T,row.names = F)
proc.time() - ptm # For threshold top 1million, 80.45s



# CLR_all_fifteen_maize
# use absolute value
agg15_ntwkClr <- abs(agg15_ntwkClr)
diag(agg15_ntwkClr) <- 0

#choose threshold.top 1 Million 
one_million <- 1e06/(15116*15116)

threshold <- quantile(agg15_ntwkClr,1-one_million)

ptm <- proc.time()
g  <- graph.adjacency((agg15_ntwkClr+t(agg15_ntwkClr))/2,mode = "undirected",weighted=TRUE)
df <- get.data.frame(g)
head(df)
df_pass1 <- subset(df,weight>threshold)
write.table(df_pass1,file = paste0("agg15_ntwkClr_",round(threshold,3),".SIF"),quote = F,sep = "\t",col.names = T,row.names = F)
proc.time() - ptm # For threshold top 1million, 80.85s
