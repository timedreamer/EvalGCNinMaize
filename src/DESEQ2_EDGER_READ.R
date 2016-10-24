#This code file is to calculate normalized gene count from the original gene count file. Also,
#do a quick check on the result expression level.
#Install two packages
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

#change workingdir
getwd()
setwd("C:\\WORK\\Co-expression")
setwd("D:\\Users\\jhuang\\Documents\\Co-expression")

# EdgeR normalization. TMM. Then CPM and RPKM.
####################################################################################################################################

#read data.
library(edgeR)
#raw_data <- read.delim("Round50.featureCount", stringsAsFactors = FALSE,check.names = FALSE)
raw_data  <- read.delim("ALL_FC_noDuplicateLib_biggerThan5Million_70allignmentRate_1266.txt",stringsAsFactors = FALSE,check.names = FALSE)

# delete unused columns
raw_data[2:5] <- NULL
ncol(raw_data) # first column is the geneid

# read to DGEList. Data filteration
gene_name <- unname(raw_data$Geneid)
ncol(raw_data) # 1268
y <- DGEList(counts = raw_data[3:1268],genes = gene_name)
keep <- rowSums(cpm(y)>2) >= 1000  # only keep some expressed data
table(keep)
y <- y[keep,,keep.lib.size=FALSE]
head(y$genes)
# Do the normalization
y <- calcNormFactors(y)
sample_name <- rownames(y$samples)

# bargraph show library size
lib_size <- sort(y$samples$lib.size,decreasing = FALSE)
mp <- barplot(lib_size*1e-6,ylab="Library size(millions)") # bargraph show library size


abline(h=5,col="red") # this is lib size after remove a lot genes.


# Get gene length of filtered genes that can be used to count RPKM.
gene_name <- as.vector(y$genes$genes) # updated gene names
id <- as.vector(unname(raw_data$Geneid))
geneLength <- raw_data[id %in% gene_name,2]
summary(geneLength)

# Calculate cpm and rpkm. log2 transformed
counts.p.m <- cpm(y,normalized.lib.sizes = TRUE,log = TRUE,prior.count = 1)
counts.p.m <- cbind(counts.p.m,geneLength) # last column is the geneLength
rownames(counts.p.m) <- gene_name


y$genes$Length <- geneLength # need geneLength info for calculating RPKM.
reads.p.m <- rpkm(y,normalized.lib.sizes = TRUE,log = TRUE,prior.count = 1)
reads.p.m <- cbind(reads.p.m,geneLength) # last column is the gene length
rownames(reads.p.m) <- gene_name

#save normalized counts
write.table(x = counts.p.m, file = "CPM_log2.txt",sep = "\t",quote = FALSE)
write.table(x = reads.p.m, file = "RPKM_log2.txt",sep = "\t",quote = FALSE)


###################################################################################
# some basic info about expression data.
# look at gene length distribution
geneLength_total <- raw_data$Length
summary(geneLength);hist(geneLength)
summary(geneLength_total)
par(mfrow=c(1,2))
hist(geneLength_total,main="total_39479")
hist(geneLength,main="afterFilter_18221")  # it seems I filter some short genes.

# gene length and expression level relationship.Use RPKM data.
par(mfrow=c(2,2))
reads.p.m <- as.data.frame(reads.p.m)
summary(reads.p.m[1:1266]) # last column is the gene length
max(reads.p.m[1:1266]);min(reads.p.m[1:1266])
plot(density(as.matrix(reads.p.m[1:1266])),main="RPKM_allGeneExprDistribution") # all gene expression distribution
avg_expr_rpkm <- apply(X = unname(reads.p.m[1:1266]),MARGIN = 1,FUN = mean)
plot(reads.p.m$geneLength,avg_expr_rpkm,main="RPKM_avgExpr vs GeneLength")
plot(density(avg_expr_rpkm),main="RPKM_density_avgExpr")
summary(avg_expr_rpkm)
##only choose gene length from 400- 3000.
exp_f1 <- reads.p.m[which(reads.p.m$geneLength>400),]
exp_f1 <- exp_f1[which(exp_f1$geneLength<3000),]
avg_expr_rpkm_range <- apply(X = unname(exp_f1[1:1266]),MARGIN = 1,FUN = mean)
length(avg_expr_rpkm_range)
plot(exp_f1$geneLength,avg_expr_rpkm_range,main="RPKM_avgExpr_400_3000")
plot(density(avg_expr_rpkm_range),main="RPKM_density_avgExpr_400_3000")

# gene length and expression level. Use CPM data.
counts.p.m <- as.data.frame(counts.p.m)
summary(counts.p.m[1:1266])
max(counts.p.m[1:1266]);min(counts.p.m[1:1266])
plot(density(as.matrix(counts.p.m[1:1266])),main="CPM_allGeneExprDistribution") # all gene expression distribution
avg_expr_cpm <- apply(X = unname(counts.p.m[1:1266]),MARGIN = 1,FUN = mean)
plot(counts.p.m$geneLength,avg_expr_cpm,main="CPM_avgExpr vs GeneLength")
plot(density(avg_expr_cpm),main="CPM_density_avgExpr")
summary(avg_expr_cpm)
##only choose gene length from 400- 3000.
exp_f2 <- counts.p.m[which(counts.p.m$geneLength>400),]
exp_f2 <- exp_f2[which(exp_f2$geneLength<3000),]
avg_expr_cpm_range <- apply(X = unname(exp_f2[1:1266]),MARGIN = 1,FUN = mean)
plot(exp_f2$geneLength,avg_expr_cpm_range,main="CPM_avgExpr_400_3000")
plot(density(avg_expr_cpm_range),main="CPM_density_avgExpr_400_3000")

# what are those highly expressed genes.
avg_exp <- apply(X = reads.p.m,MARGIN = 1,FUN = mean)
avg_exp_high <- avg_exp[which(avg_exp>8)]
write.table(avg_exp_high,file = "highexprGenes.txt",quote = FALSE,sep="\t")

# clean. Keep gene_name. for DESEQ2 use.
rm(raw_data,raw_data1,geneLength,id,keep,y)
rm(exp_f1,exp_f2,avg_exp,avg_expr_rpkm,avg_expr_rpkm_range,avg_expr_cpm,avg_expr_cpm_range,avg_exp_high)
rm(avg_expr,geneLength_total,geneLength_F1)
rm(counts.p.m,reads.p.m)



#DESEQ2_normalization as well as log2 Rawcount
##################################################################################################
library(DESeq2)
setwd("C:\\WORK\\Co-expression")
setwd("D:\\Users\\jhuang\\Documents\\Co-expression")

raw_data  <- read.delim("ALL_FC_noDuplicateLib_biggerThan5Million_70allignmentRate_1266.txt",stringsAsFactors = FALSE,check.names = FALSE)
raw_data[2:6] <- NULL
rownames(raw_data) <- raw_data$Geneid
raw_data[1] <- NULL

# filter data based on EdgeR. gene_name variable is from edgeR
# draw RawCount distribution
filter_data <- raw_data[rownames(raw_data) %in%  gene_name,]
plot(density(as.matrix(log2(filter_data+1))),main="RawCount_allGeneExprDistribution") # raw_count density
avg_expr_rc <- apply(X = unname(log2(filter_data[1:1266])),MARGIN = 1,FUN = mean)
plot(geneLength,avg_expr_rc,main="rawCount_avgExpr vs GeneLength") # gene Length and avgExpr rawCount
plot(density(avg_expr_rc),main="rawCount_density_avgExpr") # density avgExpr of rawCount

raw_log <- log2(filter_data+1)
write.table(x = raw_log,file = "rawCount_log2_result.txt",sep = "\t",quote = FALSE)

#save colnames which are file names into seperate files named "colData".
#Use colData as input for colData paramenter of DESeq
cdInput <- read.delim("colData.txt",stringsAsFactors = FALSE,check.names = FALSE,header=FALSE)
rownames(cdInput) <- cdInput$V1
dds <- DESeqDataSetFromMatrix(countData = filter_data ,colData = cdInput,design= ~1)
str(dds)

#Calculate VST and then save to seperate file
system.time(vst <- varianceStabilizingTransformation(dds,blind = TRUE)) # took 90min
head(assay(vst),3)
write.table(x = assay(vst),file = "VST_result.txt",sep = "\t",quote = FALSE)


#Draw sample tree
sampleTree = hclust(dist(t(assay(vst))), method = "average");
par(mfrow=c(1,1));sizeGrWindow(12,9)
par(cex = 0.6,mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering", sub="", xlab="",cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)


# check expresseion level
par(mfrow=c(2,2))
vst_expr <- assay(vst);vst_expr <- cbind(vst_expr,geneLength) # last column is the geneLength
plot(density(vst_expr[,1:1266]),main="VST_allGeneExprDistribution")
avg_expr_vst <- apply(X = unname(vst_expr[,1:1266]),MARGIN = 1,FUN = mean) # VST avgExpr
plot(vst_expr[,1266],avg_expr_vst,main="VST_avgExpr vs GeneLength",xlab="geneLength") # gene Length and avgExpr rawCount
plot(density(avg_expr_vst),main="VST_density_avgExpr") # density avgExpr of VST

##only choose gene length from 400- 3000.
exp_f3 <- vst_expr[which(vst_expr[,1267]>400),]  # col 1267 is the geneLength
exp_f3 <- exp_f3[which(exp_f3[,1267]<3000),]
avg_expr7 <- apply(X = unname(exp_f3[,1:1266]),MARGIN = 1,FUN = mean)
plot(exp_f3[,1267],avg_expr7,main="VST_avgExpr_400_3000",xlab="geneLength")

#clean
rm(dds,vst,cdInput,filter_data,result_vst,sampleTree,sampleDistMatrix)
rm(d,exp_f3,vst_expr)
rm(avg_expr_rc,avg_expr_vst,avg_expr7,gene_name)
