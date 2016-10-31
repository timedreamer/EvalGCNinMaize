# this file is to explore those low values after normalization.
# When plot density, there is a small peak at cpm(log2) and vst.
# At rawcount(log2), there are multiple small peaks; while rpkm(log2) is not obviouse.
# I wonder are those low values from specific genes or libraries?

# calculate how many genes above threshold.
low_exp <- function(x,thd){
  l <- length(which(x < threshold))
  return(l)

}

########CPM######################
setwd("C:\\WORK\\EvalGCNinMaize\\data")

cpm_norm <- read.delim("CPM_log2.txt",sep="\t",header=T) # geneLength last column
cpm_norm[1267] <- NULL # need to delete last column

# plot density of gene expression as a whole. small peak can be observed.
plot(density(as.matrix(cpm_norm[1:1266])),main="CPM_allGeneExprDistribution")
abline(v=-3.7,col="red") # detect where is the small peak

# plot average expresion. MARGIN can be 1(genes) or 2(SRA)
avg_expr_cpm <- apply(X = unname(cpm_norm),MARGIN = 1,FUN = mean)
plot(density(avg_expr_cpm))


# set the threshold based on abline
threshold <- -3.7

# For two margines, calculate number of low expression elements.
low_exp_sra <- apply(cpm_norm,2,low_exp,thd = threshold)
low_exp_gene <- apply(cpm_norm,1,low_exp,thd= threshold)

# Explore low expression elements.
fivenum(low_exp_sra)
sum(low_exp_sra)
plot(density(low_exp_sra))
low_exp_sra[low_exp_sra>1000]
sum(low_exp_sra[low_exp_sra>1000])
# conclusion: for SRA, low elements are widely spreaded. However, eight libaries contribute to
# over 25% of low elements. Closer look at those libraries show they are pollen tissues.


#some genes has higher number of low elements, but highest is 211 which is still considerable low.
# saved those genes and did GO enrichment in agriGO, nothing sig returned,
fivenum(low_exp_gene)
sum(low_exp_gene)
plot(density(low_exp_gene))
low_exp_gene[which(low_exp_gene>100)]
sum(low_exp_gene[low_exp_gene>100])


t <-  names(low_exp_gene[which(low_exp_gene>100)])
write.table(t,file="t_low_exp_gene.txt",sep="\t",col.names = F,row.names = F,quote = F)
# no sig GO can be retrieved.

rm(cpm_norm,avg_expr_cpm,low_exp_gene,low_exp_sra,threshold)



#################################################################################################
##How about VST? # Almost the same as CPM
vst_norm <- read.delim("VST_result.txt",sep="\t",header=TRUE) # no geneLength column

# plot density of gene expression as a whole. small peak can be observed.
plot(density(as.matrix(vst_norm)),main="VST_allGeneExprDistribution")
abline(v=-4.2,col="red") # detect where is the small peak

# plot average expresion. MARGIN can be 1(genes) or 2(SRA)
avg_expr_vst <- apply(X = unname(vst_norm),MARGIN = 1,FUN = mean)
plot(density(avg_expr_vst))


# set the threshold based on abline
threshold <- -4.2

# For two margines, calculate number of low expression elements.
low_exp_sra <- apply(vst_norm,2,low_exp,thd = threshold)
low_exp_gene <- apply(vst_norm,1,low_exp,thd= threshold)

# Explore low expression elements.
fivenum(low_exp_sra)
sum(low_exp_sra)
plot(density(low_exp_sra))
low_exp_sra[low_exp_sra>1000]
sum(low_exp_sra[low_exp_sra>1000])

rm(vst_norm,avg_expr_vst,low_exp_gene,low_exp_sra,threshold)






###############################################################################################
##Raw Count is also very similar with CPM and VST, though have higher number of low expression
# elements.
rc_norm <- read.delim("rawCount_log2_result.txt",sep="\t",header=TRUE) # no geneLength column
# plot density of gene expression as a whole. small peak can be observed.
plot(density(as.matrix(rc_norm)),main="RC_allGeneExprDistribution")
abline(v=3.2,col="red") # detect where is the small peak

# plot average expresion. MARGIN can be 1(genes) or 2(SRA)
avg_expr_rc <- apply(X = unname(rc_norm),MARGIN = 1,FUN = mean)
plot(density(avg_expr_rc))


# set the threshold based on abline
threshold <- 3.2

# For two margines, calculate number of low expression elements.
low_exp_sra <- apply(rc_norm,2,low_exp,thd = threshold)
low_exp_gene <- apply(rc_norm,1,low_exp,thd= threshold)

# Explore low expression elements.
fivenum(low_exp_sra)
sum(low_exp_sra)
plot(density(low_exp_sra))
low_exp_sra[low_exp_sra>1000]
sum(low_exp_sra[low_exp_sra>1000])
