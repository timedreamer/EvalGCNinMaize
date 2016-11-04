# This file is to cnovert correlation matrix file to SIF file.
# For correlation, change diag to 0 to avoid self-connections.
# The resulting SIF can be imported in Cytoscape.


library(igraph)

# use absolute value
cpm_pcc <- abs(cpm_pcc)
diag(cpm_pcc) <- 0

#choose threshold. For example, top 0.1%, 0.01%.
threshold <- quantile(cpm_pcc,0.999)

ptm <- proc.time()
g  <- graph.adjacency((cpm_pcc+t(cpm_pcc))/2,mode = "undirected",weighted=TRUE)
df <- get.data.frame(g)
head(df)
df_pass1 <- subset(df,weight>threshold)
write.table(df_pass1,file = paste0("cpm_pcc_",round(threshold,3),".SIF"),quote = F,sep = "\t",col.names = T,row.names = F)
proc.time() - ptm # For threshold=0.99 took ~85s.


# CLR
cpm_clr <- abs(cpm_clr)
diag(cpm_clr) <- 0

#choose threshold. For example, top 0.1%, 0.01%.
threshold <- quantile(cpm_clr,0.999)

ptm <- proc.time()
g  <- graph.adjacency((cpm_clr+t(cpm_clr))/2,mode = "undirected",weighted=TRUE)
df <- get.data.frame(g)
head(df)
df_pass1 <- subset(df,weight>threshold)
write.table(df_pass1,file = paste0("cpm_clr_",round(threshold,3),".SIF"),quote = F,sep = "\t",col.names = T,row.names = F)
proc.time() - ptm # 0.999 threshold 18.519
