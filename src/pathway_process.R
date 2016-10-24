# this file to evaluate correlation network based on pathway information from MaizeCyc
# Code tested.
getwd()
setwd("C:\\WORK\\Co-expression\\GeneOntology_pathway")

# read pathway information

pathway_raw <- read.delim("pathways.col_final_orginzed.txt",sep="\t",stringsAsFactors=FALSE,header =TRUE)
colnames(pathway_raw) <- pathway_raw[1,]
pathway_raw <- pathway_raw[-1,] # delete first row which was pathway name

# use loop to change transcript number to gene number
for (i in 1: 270){
  for (j in 1:413){
    if (startsWith(pathway_raw[[i,j]],"GRMZM")) {
      pathway_raw[[i,j]] <- substr(pathway_raw[[i,j]],1,13)
    } else if (startsWith(pathway_raw[[i,j]],"AC")){
      pathway_raw[[i,j]] <- sub(pattern = "P",replacement = "",x = pathway_raw[[i,j]])
    } else{
      pathway_raw[[i,j]] <- "NA"
    }
  }
}

write.table(pathway_raw, file="pathway_namedChanged.txt",sep="\t",quote=FALSE,row.names = FALSE)
