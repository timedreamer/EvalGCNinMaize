# 


# load data for HDA101 vst_gcc, and HDA101 target.
library(ROCR)
getwd()
setwd("C:\\WORK\\Co-expression\\ChIPSeq")
hda_target <- read.delim("HDA101_targets.txt",header=F,sep='\t')

setwd("C:\\WORK\\Co-expression\\ChIPSeq\\HDA101")
upGenes <- read.delim("HDA101_increasedExpressGenes.txt",header=T,sep="\t")
name3 <- upGenes$Gene




setwd("C:\\WORK\\EvalGCNinMaize\\data\\hda101_rank")
hda101_vst_gcc <- read.delim("hda101.rpkmclr",header=F,sep='\t')
# result without sorting is the same as sorting.

rownames(hda101_vst_gcc) <- hda101_vst_gcc[,1]
name1 <- rownames(hda101_vst_gcc)
hda101_vst_gcc[,3] <- 0

rownames(hda_target) <- hda_target[,1]
name2 <- rownames(hda_target)


int_name <- intersect(name1,name2)

tar_up <- intersect(name2,name3)

# hda101_vst_gcc[which(name1 %in% int_name),3] <- 1
# colnames(hda101_vst_gcc) <- c("Gene","Score","Target") 

hda101_vst_gcc[which(name1 %in% tar_up),3] <- 1
colnames(hda101_vst_gcc) <- c("Gene","Score","Target")

pred <- prediction(hda101_vst_gcc$Score,hda101_vst_gcc$Target)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
# plot(perf,main="hda101_rpkmClr")
# abline(0,1,col="red")
# table(hda101_vst_gcc$Target)
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
auc <- round(auc, digits = 3)
1-auc
# auct <- paste(c("AUC = "),auc,sep="")
# legend(0.65,0.3,auct,border="white",cex=1.4,box.col = "white")


##################################################################################################################
##################################################################################################################
# Calculate ROC for other result from other database.

setwd("C:\\WORK\\Co-expression\\ChIPSeq\\HDA101")

cob <- read.del











