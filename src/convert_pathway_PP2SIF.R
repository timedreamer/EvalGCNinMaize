################################################################################
# This file is to convert Protein-Protein interaction file and Pathyway file that I parsed
# Before into SIF files, which later to be used to create lists of genes with their
# co-expressed genes.

library(igraph)


##########Process Protein-Protein interaction file##############################
# final result is "PP_newtwork_cmp.SIF"
load("~/Co-expression/highConfidentPPData_geneLevel.RData")
highPPI_final[1:5,1:5]
dim(highPPI_final)

# logic is highPPI_final has 3272 columns. first to extract all unique elements from each column
# this will include one NA. 2. get the name of the column.3. for all elements in that column
# paste the column name with elements and sep by tab.4. final file is PP_newtwork.
PP_newtwork <- matrix()
for (j in c(1:3272)){
  uq_elm <- unique(highPPI_final[,j])
  cname_nf <- colnames(highPPI_final)[j]
  for (i in c(1:length(uq_elm))){
    temp1 <- paste(cname_nf,uq_elm[i],sep = "\t")
    PP_newtwork <- rbind(PP_newtwork,temp1)
  }
}

#write the file,in order to read and sep by two columns.393613 edges.
write.table(PP_newtwork,file="PP_newtwork.SIF",sep="\t",quote = F,row.names = F) # this is middle file

# read the file and now I have two columns.
PP_SIF <- read.delim("PP_newtwork.SIF",header=F,sep="\t")

# remove all NA in the file.
PP_SIF_cmp <- PP_SIF[complete.cases(PP_SIF),]
PP_SIF_cmp <- PP_SIF_cmp[-1,] # remove first row
#write final file.
write.table(PP_SIF_cmp,file="PP_newtwork_cmp.SIF",sep="\t",quote = F,row.names = F,col.names = F)



##########Process pathway file####################################################
# final file is "Pathway_network_cpm.SIF"

# For pathway information, need to find all two element vector within in each pathway and used as SIF
# use combn() function to get all combination of two. The result is "pty_sif" has 363913 edges.
setwd("D:\\Users\\jhuang\\Documents\\Co-expression")
pty <- read.delim("pathway_namedChanged.txt",sep = "\t",stringsAsFactors = FALSE)
dim(pty)

pty_sif <- matrix(nrow=1,ncol = 2)
for (i in c(1:413)){
  p_col <- pty[complete.cases(pty[,i]),i]
  pc_sif <- t(combn(p_col,m=2))
  pty_sif <- rbind(pty_sif,pc_sif)
}

pty_sif<- pty_sif[complete.cases(pty_sif),] # this just delete first row which is NA when I created.

write.table(pty_sif,file="Pathway_network_cpm.SIF",sep="\t",quote=F,row.names = F,col.names = F)
