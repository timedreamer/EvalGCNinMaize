# Method

**RNA-Seq data collection and process**

Maize genome and its annotation were downloaded from Ensembl Plant Release 31. The original 1303 RNA-Seq samples based on illumina HiSeq2000 or Hiseq2500 were downloaded from NCBI Sequence Read Archive(SRA). The downloaded files were converted to fastq format using fastq-dump command in SRA Toolkit(version 2.5.2). The adapters for the fastq files were trimmed by Cutadapt 1.8.1. The adapter-removed files were then quality checked by FastQC v0.11.2. To allign short sequences to maize genome(Refv3), the HISAT2 v2.0.4 was used. Gene-level expression counts were calculated by FeatureCounts 1.5.0 from alligned bam files.

As a result, 26 libraries with less than 5 milliion reads were excluded as well as
11 libraries with less than 70% of total allignment rate reported by HISAT2.

**Gene count normalization**

The obtain gene count expression data was normalization using four different methods before constructing GCN. Counts Per Milliion(CPM) and Reads Per Killobase Per Million(RPKM) were calcualted by edgeR package in R environment and then log2 normalized; Variance Stabilizating Transformation(VST) was calcualted by DESeq2 package. Raw count(RC) for gene expressions were normalized by log2. Low expressed genes which defined as CPM less than 2 in more than 1000 samples were exlcuded from the following analysis.

**Network construction**

Five correlation coefficient methods and four Mutual Information methods were applied to normalized gene expression data to construct GCN. All computing steps were done in R 3.3.1 environment. Pearson Correlation Coefficient(PCC) and Spearman Correlation Coefficient(SCC) was calculated by cor() function. Kendall rank Correlation Coefficient was calculated using cor.fk() function in pcaPP library. Gini Correlation Coefficient was calculated by adjacencymatrix() function in rsgcc library. Biweight midcorrelation was computed by bicor() function in WGCNA pacakge. Mutual information results were all computed by parmigene package.

**Network Performance Evaluation**

We used three four datasets to evaluate different networks' performance. First, maize protein-protein interactons were downloaded from PPIM(version 1.1). Only high-confidence interactions were left for evaluation. Second, maize pathway information was downloaded from MaizeCyc(version 2.2) and same pathway genes were considered as co-expressed together. Third, maize gene ontology data for AGPv3.30 was downloaded from AgriGO. Fourth, ChIP-Seq confirmed targets for HDA101(GRMZM2G172883) was used as poitive co-exressed examples for evaluation.

In terms of evaluation method, we adopted the widely used Receiver operating characteristic (ROC) Curve and Area Under Curve(AUC) for binary classfier problems. R pacakge ROCR and EGAD were used for this purpose.

**Network Clustering**

For each network, top 1 million edges were selected as stringent co-expression networks. The topology of networks were calcaulted using powerRlaw package in R and Cy. To test whether the distribution of nodes were following power-law distribution, the Kolomogorov Smirnoff statistic (KS) was applied with 1000 bootstrap simulations.

Graph clustering was performed using Markov Cluster Algorithm (MCL) in Cytoscape clusterMake with inflatino value set to 1.8. Networks were also visualized in Cytoscape.
