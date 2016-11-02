# This script is to calcualte AUC from Gene ontology
# for individual networks.

library (EGAD)

# read table for maize Gene ontology
gotable <- read.table("ZeamaysAGOv3.30GeneOntology.txt",header=TRUE,sep="\t")
gotable[,1] <- NULL

genes <- unique(gotable$ID)
annotationlist <- unique(gotable$GO)

annotataions_maize <- make_annotations(gotable,genes,annotationlist)

load("15116_geneNames.RData")
annotations_sub <- filter_network(annotataions_maize,flag = 2,min=20,max = 300)
rm(gotable,annotataions_maize)

#Calc. took 261s for vst_pcc
GO_12_pcc <- run_GBA(cpm12_pcc,annotations_sub);GO_12_mrnet <- run_GBA(cpm12_mrnet,annotations_sub)
GO_12_clr <- run_GBA(cpm12_clr,annotations_sub)

GO_36_pcc <- run_GBA(cpm36_pcc,annotations_sub);GO_36_mrnet <- run_GBA(cpm36_mrnet,annotations_sub)
GO_36_clr <- run_GBA(cpm36_clr,annotations_sub)

GO_65_pcc <- run_GBA(cpm65_pcc,annotations_sub);GO_65_mrnet <- run_GBA(cpm65_mrnet,annotations_sub)
GO_65_clr <- run_GBA(cpm65_clr,annotations_sub)

GO_108_pcc <- run_GBA(cpm108_pcc,annotations_sub);GO_108_mrnet <- run_GBA(cpm108_mrnet,annotations_sub)
GO_108_clr <- run_GBA(cpm108_clr,annotations_sub)

GO_270_pcc <- run_GBA(cpm270_pcc,annotations_sub);GO_270_mrnet <- run_GBA(cpm270_mrnet,annotations_sub)
GO_270_clr <- run_GBA(cpm270_clr,annotations_sub)

GO_404_pcc <- run_GBA(cpm404_pcc,annotations_sub);GO_404_mrnet <- run_GBA(cpm404_mrnet,annotations_sub)
GO_404_clr <- run_GBA(cpm404_clr,annotations_sub)

save(GO_12_pcc,GO_12_mrnet,GO_12_clr,
     GO_36_pcc,GO_36_mrnet,GO_36_clr,
     GO_65_pcc,GO_65_mrnet,GO_65_clr,
     GO_108_pcc,GO_108_mrnet,GO_108_clr,
     GO_270_pcc,GO_270_mrnet,GO_270_clr,
     GO_404_pcc,GO_404_mrnet,GO_404_clr,
     file="GO_AUROC_SixAggregate_ntwk.RData")


# save mean value of each method and each size.
GO_12_avg <- c(GO_12_pcc[[3]],GO_12_mrnet[[3]],GO_12_clr[[3]])
GO_36_avg <- c(GO_36_pcc[[3]],GO_36_mrnet[[3]],GO_36_clr[[3]])
GO_65_avg <- c(GO_65_pcc[[3]],GO_65_mrnet[[3]],GO_65_clr[[3]])
GO_108_avg <- c(GO_108_pcc[[3]],GO_108_mrnet[[3]],GO_108_clr[[3]])
GO_270_avg <- c(GO_270_pcc[[3]],GO_270_mrnet[[3]],GO_270_clr[[3]])
GO_404_avg <- c(GO_404_pcc[[3]],GO_404_mrnet[[3]],GO_404_clr[[3]])


save(GO_12_avg,GO_36_avg,GO_65_avg,GO_108_avg,GO_270_avg,GO_404_avg,
     file="GO_AUROC_meanAggregate_ntwk.RData")

rm(cpm12_pcc,cpm12_mrnet,cpm12_clr,cpm36_pcc,cpm36_mrnet,cpm36_clr,
   cpm65_pcc,cpm65_mrnet,cpm65_clr,cpm108_pcc,cpm108_mrnet,cpm108_clr,
   cpm270_pcc,cpm270_mrnet,cpm270_clr,cpm404_pcc,cpm404_clr,cpm404_mrnet)
