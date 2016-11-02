# PCA analysis for RNA-Seq samples in GCN building.
# follow this http://www-huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html

setwd("C:\\WORK\\EvalGCNinMaize\\data")

cpm_norm <- read.delim("CPM_log2.txt",sep="\t",header=T) # geneLength last column
cpm_norm[1267] <- NULL # need to delete last column

fit <- prcomp(t(cpm_norm),scale. = F,retx =T)

percentVar <- round(100*fit$sdev^2/sum(fit$sdev^2),1)

plot(fit$x[,1],fit$x[,2],xlab=paste0("PC1, VarExp:", round(percentVar[1],4)),
     ylab = paste0("PC2, VarExp:", round(percentVar[2],4)))


plot(fit) # show variance explained.


# an naive way to show what those samples are.
pc1_higher100 <- names(fit$x[,1][which(fit$x[,1] > 100)])
