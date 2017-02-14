# This scipt is to draw correlation matrix plot that show correlation methods are 
# more similar with each other than MI methods.
# the graph is without legend, delete one of "gudie=F" and draw one graph to
# retain the legend.

# correlation map for CPM
cpm_go <- read.table("maize_cpm_1266_allGOMatrix.txt",sep="\t",header=T)
row.names(cpm_go) <- cpm_go[,1];cpm_go<- cpm_go[,-1]
cpm_pppty <- read.table('maize_cpm_1266_allPPPTYMatrix_nullDelete.txt')

vst_go <- read.table("maize_vst_1266_allGOMatrix.txt",sep="\t",header=T)
row.names(vst_go) <- vst_go[,1];vst_go<- vst_go[,-1]
vst_pppty <- read.table('maize_vst_1266_allPPPTYMatrix__nullDelete.txt')

rpkm_go <- read.table("maize_rpkm_1266_allGOMatrix.txt",sep="\t",header=T)
row.names(rpkm_go) <- rpkm_go[,1];rpkm_go<- rpkm_go[,-1]
rpkm_pppty <- read.table('maize_rpkm_1266_allPPPTYMatrix_nullDelete.txt',sep='\t',header=T)

corData <- melt(cor(cpm_go,method = 'p'))
corData2 <- melt(cor(cpm_pppty,method = 'p'))
p1 <- ggplot(corData, aes(x=Var1, y=Var2, fill=value)) +geom_tile() + xlab("") + ylab("") + 
  scale_fill_gradient2(low = "navyblue", high = "red", mid = "white",midpoint = 0.5,
                       guide=F,limit = c(0,1), space = "Lab") + 
  theme(text = element_text(size=20),
        axis.text.x=element_text(angle = 90, hjust = 0)) + ggtitle("CPM_GO") 
p2 <- ggplot(corData2, aes(x=Var1, y=Var2, fill=value)) +geom_tile() + xlab("") + ylab("") + 
  scale_fill_gradient2(low = "navyblue", high = "red", mid = "white",midpoint = 0.5,
                       guide=F,limit = c(0,1), space = "Lab") + 
  theme(text = element_text(size=20),
        axis.text.x=element_text(angle = 90, hjust = 0)) + ggtitle("CPM_PPPTY")

# get_upper_tri<-function(cormat){
#   cormat[upper.tri(cormat)] <- NA
#   return(cormat)
# } 


corData <- melt(cor(vst_go,method = 'p'))
corData2 <- melt(cor(vst_pppty,method = 'p'))
p3 <- ggplot(corData, aes(x=Var1, y=Var2, fill=value)) +geom_tile() + xlab("") + ylab("") + 
  scale_fill_gradient2(low = "navyblue", high = "red", mid = "white",midpoint = 0.5, 
                       guide=F,limit = c(0,1), space = "Lab") + 
  theme(text = element_text(size=20),
        axis.text.x=element_text(angle = 90, hjust = 0)) + ggtitle("VST_GO") 
p4 <- ggplot(corData2, aes(x=Var1, y=Var2, fill=value)) +geom_tile() + xlab("") + ylab("") + 
  scale_fill_gradient2(low = "navyblue", high = "red", mid = "white",midpoint = 0.5, 
                       guide=F,limit = c(0,1), space = "Lab") + 
  theme(text = element_text(size=20),
        axis.text.x=element_text(angle = 90, hjust = 0)) + ggtitle("VST_PPPTY")

corData <- melt(cor(rpkm_go,method = 'p'))
corData2 <- melt(cor(rpkm_pppty,method = 'p'))
p5 <- ggplot(corData, aes(x=Var1, y=Var2, fill=value)) +geom_tile() + xlab("") + ylab("") + 
  scale_fill_gradient2(low = "navyblue", high = "red", mid = "white",midpoint = 0.5, 
                       guide=F,limit = c(0,1), space = "Lab") + 
  theme(text = element_text(size=20),
        axis.text.x=element_text(angle = 90, hjust = 0)) + ggtitle("RPKM_GO") 
p6 <- ggplot(corData2, aes(x=Var1, y=Var2, fill=value)) +geom_tile() + xlab("") + ylab("") + 
  scale_fill_gradient2(low = "navyblue", high = "red", mid = "white",midpoint = 0.5,
                       guide=F,limit = c(0,1), space = "Lab") + 
  theme(text = element_text(size=20),
        axis.text.x=element_text(angle = 90, hjust = 0)) + ggtitle("RPKM_PPPTY") 

multiplot(p1,p2,p3,p4,p5,p6,cols =3 )

