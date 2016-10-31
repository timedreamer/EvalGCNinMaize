# This file is to find what are the highly connected genes(hub genes) in each network when choose
# a threshold.


#load network



# this function is to calculate how many genes connected.
calsumR <- function(x){
  passThreshold <- length(which(x>threshold))
  return(passThreshold)
}

# Here is the funciton what has three parameters.
# x is the correlation attributes(cpm_pcc)
# thd is the threhold for the correlation matrix(0.9/0.99)
# ctf is the cutoff value for choosing top percentage of highly connected genes(0.9)

highConnectGenes <- function(x,thd,ctf){

  threshold <- quantile(x,thd)
  passT1 <- apply(x,2,calsumR);passT1 <- passT1[which(passT1>1)]
  cutoff <- quantile(passT1,ctf)
  passT2 <- passT1[which(passT1>=cutoff)]
  highGeneName <- attributes(passT2)
  write.table(highGeneName$names,file=paste0(deparse(substitute(x)),"_",round(threshold,3),"_highConnetGene.txt"),
              sep="\t",row.names = F,quote = F,col.names = F)

}


highConnectGenes(cpm_pcc,0.999,0.9)
highConnectGenes(cpm_mrnet,0.999,0.9)
