##Change the format of the Pathway and Protein Protein interaction data. My purpose is to create a list
# with gene node and their co-expressed genes which included in the pathway or PP.
# for this purpose, I combined pathway_SIF and PP_SIF which created from previous code. Then delete genes
# that have less than 5 nodes connnected. Last use conver_list to
# extract all co-expressed genes in list.
getwd()
setwd("D:/Users/jhuang/Documents/Co-expression")

convert_list <- function(x){
  return(P_new_combine[grep(x,P_new_combine$V1),2])
}


PP_SIF <- read.delim("PP_newtwork_cmp.SIF",header=F,sep="\t")
PTY_SIF <- read.delim("Pathway_network_cpm.SIF",header=F,sep="\t")

P_combine <- rbind(PP_SIF,PTY_SIF)
P_combine <- P_combine[order(P_combine$V1),]

# choose Genes that have more than 71 genes connected. 1777 genes left
length(which(table(P_combine$V1)>=5))
fivenum(table(P_combine$V1))
P_new_combine <- P_combine[!(as.numeric(P_combine$V1)) %in% which(table(P_combine$V1)<5),]

# convert two columns into character
P_new_combine[,1] <- as.character(P_new_combine[,1]);P_new_combine[,2] <- as.character(P_new_combine[,2])

node_name <- unique(P_new_combine[,1])
node_number <- length(node_name)

cogene_result <- sapply(node_name,convert_list) # this step takes several minutes
rm(P_combine,P_new_combine,PP_SIF,PTY_SIF,convert_list)
save(cogene_result,file= paste0(node_number,"_PP_PTY_combine_list.RData"))
save(node_name,file=paste0(node_number,"_nodeName.RData"))
