library(rentrez)
# Just search one year of SRA samples.
entrez_search(db="sra",
              term="(Arabidopsis[ORGN] OR Arabidopsis thaliana[ORGN]) AND 2015[PDAT] AND RNA-Seq[STRA]",
              retmax=0)

search_year_sra <- function(year, term){
  query <- paste(term, "AND (", year, "[PDAT])")
  entrez_search(db="sra", term=query, retmax=0)$count
}

year <- 2008:2016
sra_arabidopsis <- sapply(year, search_year_sra, term="(Arabidopsis[ORGN] OR Arabidopsis thaliana[ORGN]) AND RNA-Seq[STRA] AND illumina[PLAT]",
                    USE.NAMES=FALSE)

entrez_search(db="gds",
              term="GPL198[ACCN] AND 2016[PDAT]",
              retmax=0)

# 1022, 1655, 1283, 1241, 1100, 1273,1216, 1015, 802
# from 2008 to 2016



setwd("D:\\Users\\jhuang\\Documents\\Co-expression\\arabidopsis")
library(edgeR)
#Raw FPKM was downloaded from PODC.
raw_data  <- read.delim("arabi.txt",stringsAsFactors = FALSE,check.names = FALSE,header=T)
colnames(raw_data)[1] <- "gene_name"

colnames(raw_data)[1364]

filter_data <- raw_data[rowSums(raw_data>2)>=1000,]
filter_data[1:5,1:5]

filter_data[,1] <- substr(filter_data[,1],1,9)
noDupli_data <- filter_data[!duplicated(filter_data$gene_name),]

rownames(noDupli_data) <- noDupli_data[,1]
noDupli_data[,1] <- NULL
noDupli_data[1:5,1:5]

save(noDupli_data,file="arabidopsis_expression.RData")

ara_pcc <- cor(t(noDupli_data),method = "p")
save(ara_pcc,file="arabidopsis_pcc.RData")

library (EGAD)

gotable <- read.table("gene_association.tair_parsed.txt",header=T,sep="\t")
gotable[1:5,1:2]

genes <- unique(gotable$ID)
annotationlist <- unique(gotable$GO)

annotataions_arab <- make_annotations(gotable,genes,annotationlist)
GO_arab_pcc <- run_GBA(ara_pcc,annotataions_arab)
save(GO_arab_pcc,file="arabidopsis_GO_pcc.RData")


