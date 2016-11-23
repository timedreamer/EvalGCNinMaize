# this script is to retrieve how many Affy microarrya and illumina RNASeq
# have been produced every year(2008-2015) for maize gene expression analysis.
# Karen wants to include Nimblgen array(GPL12620)

library(rentrez)

# what key words are searchable.
entrez_dbs()
entrez_db_summary("sra")
entrez_db_summary("geoprofiles")

entrez_db_searchable("sra")
entrez_db_searchable("gds")


# Just search one year of SRA samples.
entrez_search(db="sra",
              term="(Zea mays[ORGN] OR maize[ORGN]) AND 2015[PDAT] AND RNA-Seq[STRA])",
              retmax=0)

# Search one year of GEO DataSet(GDS).
# GPL4032(Affy);GPL12620(Nimblegen)
entrez_search(db="gds",
              term="(GPL4032[ACCN] OR GPL12620[ACCN]) AND 2015[PDAT]",
              retmax=0)


# use function and sapply to calculate year 2008-2015
search_year_sra <- function(year, term){
  query <- paste(term, "AND (", year, "[PDAT])")
  entrez_search(db="sra", term=query, retmax=0)$count
}


year <- 2008:2015
sra_maize <- sapply(year, search_year_sra, term="(Zea mays[ORGN] OR maize[ORGN]) AND RNA-Seq[STRA] AND illumina[PLAT]",
                 USE.NAMES=FALSE)


# Final result. year 2008-2015.2007-geo:88;2006-geo:1; previous year 0.[GPL4032 Only]
geo_maize <- c(177,26,95,108,60,34,69,14)
# Both GPL4032 and GPL12620
geo_maize <- c(177,26,95,290,60,34,69,14)
sra_maize <- c(0,6,9,31,460,524,811,2389)


png("seq_array_1.png", units = "px", width=3000, height=3000, res=300)
par(mar=c(6,6,3,1))
plot(sra_maize,type="b",col="black",xaxt="n",ylim=c(0,2500), main="Maize RNA-Seq vs Microarray Samples",
     ylab="number of samples",xlab="year",lwd=2.5,cex.axis =1.5,cex.lab=1.5,cex.main=2)
lines(geo_maize,type="b",col="black",lty=2,lwd=2.5)
axis(1, at=1:8, labels=c(2008:2015),cex.axis=1.5,cex.lab=2)
dev.off()


# Or use Excel to draw graph
geo_maize;sra_maize

geo_total <- 718 # plus GPL12620 is 900

geo_total <- 900

sra_total <- 5056
c_total <- c(geo_total,sra_total)
names(c_total) <- c("Microarray","RNA-Seq")

png("seq_array_2.png", units = "px", width=3000, height=3000, res=300)
par(mar=c(6,10,10,6))
barplot(c_total,col="black",main="Total maize RNA-Seq vs Microarray",
        ylab="number of samples",cex.main=2.5,cex.lab=1.5,cex.axis=2,cex.names = 2)
dev.off()