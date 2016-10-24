#This file is to change the format of high confident protein-protein interaction
#in maize. Easier for me to use both Excel and R.

# A better way would be like what I processed the SIF files,
# convet two columns data to list. I will modify when I have time.

getwd()
setwd("C:\\WORK\\Co-expression\\Co-expression_Network")
# The table I downloaded from http://comp-sysbio.org/ppim/
# this is high confident prediction.
hPPI <- read.table("highConfidentPPIs.txt", header=FALSE)
t_hPPI <- t(hPPI)
colnames(t_hPPI) <- t_hPPI[1,]
head(t_hPPI)

# remove the query gene, the drop here is essential!
query_name <- t_hPPI[-1,,drop=F]
# convert it to list. With query gene as list name,
# and proteins combined with it as list content.
list_pp <- sapply(unique(colnames(query_name)), function(x) unname(unlist(query_name[,colnames(query_name)==x])))

# This function I found online that can write list into tables.
# However this one is not perfect, because the list content is still
# only seperated by space.
fnlist <- function(x, fil){ z <- deparse(substitute(x))
                            cat(z, "\n", file=fil)
                            nams=names(x)
                            for (i in seq_along(x) ){ cat(nams[i], "\t",  x[[i]], "\r",
                                                          file=fil, append=TRUE) }
}

fnlist(list_pp,'protein_high_prediction.txt')

# Because excel has a limitation of characters in a cell, so when use execel to split data,
# got very complicated. But bascially I did it manually.

##########################################################################################
# Now read the file that I manually curated. Code tested.
highPP_curated <- read.table("highPPI.csv.txt",sep="\t",header=F,na.strings =c(""," ","NA") )

# Change the transcript level to gene level.
highPP_curated <- apply(highPP_curated,2,function(x) sub("FGP","FG",x)) # change genes with "AC"
highPP_curated <- apply(highPP_curated,2,function(x) sub("_P\\d+","",x))  # change genes with "GR"
f[1:5,1:5] #check data


# take a look at the data.
highPPI_final <- t(highPP_curated)
colnames(highPPI_final) <- highPPI_final[1,]
colnames(highPPI_final)[1:5]
highPPI_final<- highPPI_final[-1,,drop=F]  # new_f is the file that I need with column name the bait protein name,
                           # rest of the column are the proteins that interact with it.
save(highPPI_final,file="highConfidentPPData_geneLevel.RBD")
# This new_f file will be the high confident prediction of protein-p interaction.
# the column name is the query gene number, in total 3272.
# rows are proteins that interact with it. Because I need to keep it as a matrix,
# in order to calculate sum/mean etc, so fill blank cell with "NA".

non_na <- apply(highPPI_final,2,function(x) length(which(!is.na(x))))
table(non_na) # most of them just have one protein interaction.
hist(non_na,main="total highPP",xlab="number of interacted proteins")
sum(non_na) #155845
max(non_na)
fivenum(non_na)
d <- density(non_na)
hist(non_na)
plot(d)
