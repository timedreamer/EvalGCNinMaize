# this script is to

library (EGAD)

# read table for maize Gene ontology
gotable <- read.table("ZeamaysAGOv3.30GeneOntology.txt",header=TRUE,sep="\t")
gotable[,1] <- NULL

genes <- unique(gotable$ID)
annotationlist <- unique(gotable$GO)

annotataions_maize <- make_annotations(gotable,genes,annotationlist)

load("15116_geneNames.RData")
genelist <- rownames(vst_gcc)
annotations_sub <-filter_network(annotataions_maize,flag = 2,min=20,max = 300)
rm(gotable)


#VST. took 261s for vst_pcc
GO_vst_pcc <- run_GBA(vst_pcc,annotations_sub);GO_vst_scc <- run_GBA(vst_scc,annotations_sub);
GO_vst_gcc <- run_GBA(vst_gcc,annotations_sub);GO_vst_kcc <- run_GBA(vst_kcc,annotations_sub);
GO_vst_bicor <- run_GBA(vst_bicor,annotations_sub);GO_vst_aa <- run_GBA(vst_aaracne,annotations_sub);
GO_vst_ma <- run_GBA(vst_maracne,annotations_sub);GO_vst_mrnet <- run_GBA(vst_mrnet,annotations_sub);
GO_vst_clr <- run_GBA(vst_clr,annotations_sub);
save(GO_vst_pcc,GO_vst_scc,GO_vst_kcc,GO_vst_gcc,GO_vst_bicor,
     GO_vst_aa,GO_vst_ma,GO_vst_mrnet,GO_vst_clr,file="GO_vst.RData")
rm(GO_vst_pcc,GO_vst_scc,GO_vst_kcc,GO_vst_gcc,GO_vst_bicor,
   GO_vst_aa,GO_vst_ma,GO_vst_mrnet,GO_vst_clr,vst_pcc,vst_scc,vst_kcc,vst_gcc,
   vst_bicor,vst_aaracne,vst_maracne,vst_mrnet,vst_clr)

#CPM
GO_cpm_pcc <- run_GBA(cpm_pcc,annotations_sub);GO_cpm_scc <- run_GBA(cpm_scc,annotations_sub);
GO_cpm_gcc <- run_GBA(cpm_gcc,annotations_sub);GO_cpm_kcc <- run_GBA(cpm_kcc,annotations_sub);
GO_cpm_bicor <- run_GBA(cpm_bicor,annotations_sub);GO_cpm_aa <- run_GBA(cpm_aaracne,annotations_sub);
GO_cpm_ma <- run_GBA(cpm_maracne,annotations_sub);GO_cpm_mrnet <- run_GBA(cpm_mrnet,annotations_sub);
GO_cpm_clr <- run_GBA(cpm_clr,annotations_sub);
save(GO_cpm_pcc,GO_cpm_scc,GO_cpm_kcc,GO_cpm_gcc,GO_cpm_bicor,
     GO_cpm_aa,GO_cpm_ma,GO_cpm_mrnet,GO_cpm_clr,file="GO_cpm.RData")
rm(GO_cpm_pcc,GO_cpm_scc,GO_cpm_kcc,GO_cpm_gcc,GO_cpm_bicor,
   GO_cpm_aa,GO_cpm_ma,GO_cpm_mrnet,GO_cpm_clr,cpm_pcc,cpm_scc,cpm_kcc,cpm_gcc,
   cpm_bicor,cpm_aaracne,cpm_maracne,cpm_mrnet,cpm_clr)

#RC
GO_rc_pcc <- run_GBA(rc_pcc,annotations_sub);GO_rc_scc <- run_GBA(rc_scc,annotations_sub);
GO_rc_gcc <- run_GBA(rc_gcc,annotations_sub);GO_rc_kcc <- run_GBA(rc_kcc,annotations_sub);
GO_rc_bicor <- run_GBA(rc_bicor,annotations_sub);GO_rc_aa <- run_GBA(rc_aaracne,annotations_sub);
GO_rc_ma <- run_GBA(rc_maracne,annotations_sub);GO_rc_mrnet <- run_GBA(rc_mrnet,annotations_sub);
GO_rc_clr <- run_GBA(rc_clr,annotations_sub);
save(GO_rc_pcc,GO_rc_scc,GO_rc_kcc,GO_rc_gcc,GO_rc_bicor,
     GO_rc_aa,GO_rc_ma,GO_rc_mrnet,GO_rc_clr,file="GO_rc.RData")
rm(GO_rc_pcc,GO_rc_scc,GO_rc_kcc,GO_rc_gcc,GO_rc_bicor,
   GO_rc_aa,GO_rc_ma,GO_rc_mrnet,GO_rc_clr,rc_pcc,rc_scc,rc_kcc,rc_gcc,
   rc_bicor,rc_aaracne,rc_maracne,rc_mrnet,rc_clr)

#RPKM
GO_rpkm_pcc <- run_GBA(rpkm_pcc,annotations_sub);GO_rpkm_scc <- run_GBA(rpkm_scc,annotations_sub);
GO_rpkm_gcc <- run_GBA(rpkm_gcc,annotations_sub);GO_rpkm_kcc <- run_GBA(rpkm_kcc,annotations_sub);
GO_rpkm_bicor <- run_GBA(rpkm_bicor,annotations_sub);GO_rpkm_aa <- run_GBA(rpkm_aaracne,annotations_sub);
GO_rpkm_ma <- run_GBA(rpkm_maracne,annotations_sub);GO_rpkm_mrnet <- run_GBA(rpkm_mrnet,annotations_sub);
GO_rpkm_clr <- run_GBA(rpkm_clr,annotations_sub);
save(GO_rpkm_pcc,GO_rpkm_scc,GO_rpkm_kcc,GO_rpkm_gcc,GO_rpkm_bicor,
     GO_rpkm_aa,GO_rpkm_ma,GO_rpkm_mrnet,GO_rpkm_clr,file="GO_rpkm.RData")
rm(GO_rpkm_pcc,GO_rpkm_scc,GO_rpkm_kcc,GO_rpkm_gcc,GO_rpkm_bicor,
   GO_rpkm_aa,GO_rpkm_ma,GO_rpkm_mrnet,GO_rpkm_clr,rpkm_pcc,rpkm_scc,rpkm_kcc,rpkm_gcc,
   rpkm_bicor,rpkm_aaracne,rpkm_maracne,rpkm_mrnet,rpkm_clr)


#record average AUROC for different network.
GO_vst_avg <- c(GO_vst_pcc[[3]],GO_vst_scc[[3]],GO_vst_kcc[[3]],GO_vst_gcc[[3]],GO_vst_bicor[[3]],
GO_vst_aa[[3]],GO_vst_ma[[3]],GO_vst_mrnet[[3]],GO_vst_clr[[3]])

GO_cpm_avg <- c(GO_cpm_pcc[[3]],GO_cpm_scc[[3]],GO_cpm_kcc[[3]],GO_cpm_gcc[[3]],GO_cpm_bicor[[3]],
GO_cpm_aa[[3]],GO_cpm_ma[[3]],GO_cpm_mrnet[[3]],GO_cpm_clr[[3]])

GO_rc_avg <- c(GO_rc_pcc[[3]],GO_rc_scc[[3]],GO_rc_kcc[[3]],GO_rc_gcc[[3]],GO_rc_bicor[[3]],
GO_rc_aa[[3]],GO_rc_ma[[3]],GO_rc_mrnet[[3]],GO_rc_clr[[3]])

GO_rpkm_avg <- c(GO_rpkm_pcc[[3]],GO_rpkm_scc[[3]],GO_rpkm_kcc[[3]],GO_rpkm_gcc[[3]],GO_rpkm_bicor[[3]],
GO_rpkm_aa[[3]],GO_rpkm_ma[[3]],GO_rpkm_mrnet[[3]],GO_rpkm_clr[[3]])

save(GO_vst_avg,GO_cpm_avg,GO_rc_avg,GO_rpkm_avg,file="GO_eval_avg_all.RData")
