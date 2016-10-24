getwd()
setwd("C:\\WORK\\Co-expression\\Co-expression_Network")
# The table I downloaded from http://comp-sysbio.org/ppim/
# this is high confident prediction.
hPPI <- read.table("highConfidentPPIs.txt", header=FALSE)
t_hPPI <- t(hPPI)
colnames(t_hPPI) <- t_hPPI[1,]
head(t_hPPI)

# remove the query gene, the drop here is essential!
n <- t_hPPI[-1,,drop=F]
# convert it to list. With query gene as list name,
# and proteins combined with it as list content.
ttt <- sapply(unique(colnames(n)), function(x) unname(unlist(n[,colnames(n)==x])))
ttt[1]
ttt[2]

# This function I found online that can write list into tables.
# However this one is not perfect, because the list content is still
# only seperated by space.
fnlist <- function(x, fil){ z <- deparse(substitute(x))
                            cat(z, "\n", file=fil)
                            nams=names(x) 
                            for (i in seq_along(x) ){ cat(nams[i], "\t",  x[[i]], "\r", 
                                                          file=fil, append=TRUE) }
}

fnlist(ttt,'protein_high_prediction.txt')

# Because excel has a limitation of characters in a cell, so when use execel to split data,
# got very complicated. But bascially I did it manually.

# Now read the file that I manually curated. 
f <- read.table("highPPI.csv.txt",sep="\t",header=F,na.strings =c(""," ","NA") )
ncol(f)
nrow(f)




# Change the transcript level to gene level.
f <- apply(f,2,function(x) sub("FGP","FG",x)) # change genes with "AC"
f <- apply(f,2,function(x) sub("_P\\d+","",x))  # change genes with "GR"
f[1:5,1:5] #check data


# take a look at the data.
new_f <- t(f)
colnames(new_f) <- new_f[1,]
colnames(new_f)[1:5]
new_f[1:5,1:5]
new_f<- new_f[-1,,drop=F]  # new_f is the file that I need with column name the bait protein name,
                           # rest of the column are the proteins that interact with it.
dim(new_f)
save(new_f,file="highConfidentPPData_geneLevel.RBD")
# This new_f file will be the high confident prediction of protein-p interaction.
# the column name is the query gene number, in total 3272.
# rows are proteins that interact with it. Because I need to keep it as a matrix,
# so fill black cell with "NA".

non_na <- apply(new_f,2,function(x) length(which(!is.na(x))))
table(non_na) # most of them just have one protein interaction.
hist(non_na,main="total highPP",xlab="number of interacted proteins")
sum(non_na) #155845
max(non_na)
fivenum(non_na)
d <- density(non_na)
hist(non_na)
plot(d)






################################################################
#Test code here# No use
f[1:5,1:5]
test <- f[1:5,1:5]
test <- apply(test,2,function(x) sub("FGP","FG",x))
test[1:5,1:5]
test <- apply(test,2,function(x) sub("_P\\d+","",x))
test[1:5,1:5]
################################################################
