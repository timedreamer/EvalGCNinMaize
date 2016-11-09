# PCA analysis for RNA-Seq samples in GCN building.
# follow this http://www-huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html

setwd("C:\\WORK\\EvalGCNinMaize\\data")

cpm_norm <- read.delim("CPM_log2.txt",sep="\t",header=T) # geneLength last column
cpm_norm[1267] <- NULL # need to delete last column
cpm_norm[1:5,1:5]

# Draw by default R plot
fit <- prcomp(t(cpm_norm),scale. = F,retx =T)
# this is how calculate the percentage of each principle.
percentVar <- round(100*fit$sdev^2/sum(fit$sdev^2),1)
plot(fit$x[,1],fit$x[,2],xlab=paste0("PC1, VarExp:", round(percentVar[1],4)),
     ylab = paste0("PC2, VarExp:", round(percentVar[2],4)))
plot(fit) # show variance explained.
# an naive way to show what those samples are.
pc1_higher100 <- names(fit$x[,1][which(fit$x[,1] > 100)])


# Better to draw using ggplot2
library(ggfortify)
setwd("D:\\Dropbox\\Project_related\\Library_SRA_info")
# the order of libs are the same with cpm_norm
libs <- read.table("ALL_FILE_INCLUDE_1266_lib_info.txt",header=T,sep="\t")

# Color-blind friendly palette. https://goo.gl/6rX2j8
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")

# Color different tissue
p <- autoplot(prcomp(t(cpm_norm),scale. = F,retx =T), data = libs, colour = 'Tissue_short')
p+scale_color_manual(values=cbPalette) + labs(x = "PC1(28%)",y="PC2(14.6%)") + 
  theme(axis.title=element_text(size="22"))

# Color different genotypes
p_genotype <- autoplot(prcomp(t(cpm_norm),scale. = F,retx =T), data = libs, colour = 'Genotype_short')
p_genotype + labs(x = "PC1(28%)",y="PC2(14.6%)") + theme(axis.title=element_text(size="22"))

# Plot no color.
p_noColor <- autoplot(prcomp(t(cpm_norm),scale. = F,retx = T))
p_noColor + labs(x = "PC1(28%)",y="PC2(14.6%)") + theme(axis.title=element_text(size="22"))

