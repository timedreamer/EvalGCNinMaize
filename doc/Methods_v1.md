# Method

**RNA-Seq data collection and process**

Maize genome and its annotation were downloaded from Ensembl Plant Release 31. The original 1303 RNA-Seq samples based on illumina HiSeq2000 or Hiseq2500 were downloaded from NCBI Sequence Read Archive(SRA). The downloaded files were converted to fastq format using fastq-dump command in SRA Toolkit(version 2.5.2). The adapters for the fastq files were trimmed by Cutadapt 1.8.1. The adapter-removed files were then quality checked by FastQC v0.11.2. To allign short sequences to maize genome(Refv3), the HISAT2 v2.0.4 was used. Gene-level expression counts were calculated by FeatureCounts 1.5.0 from alligned bam files.

As a result, 26 libraries with less than 5 milliion reads were excluded as well as
11 libraries with less than 70% of total allignment rate reported by HISAT2.

**Gene count normalization**

The obtain gene count expression data was normalization using four different methods before constructing GCN. Counts Per Milliion(CPM) and Reads Per Killobase Per Million(RPKM) were calcualted by edgeR package in R environment and then log2 normalized; Variance Stabilizating Transformation(VST) was calcualted by DESeq2 package. Raw count(RC) for gene expressions were normalized by log2.

**Network construction**

Five correlation coefficient methods and four Mutual Information methods were applied to normalized gene expression data to construct GCN. Pearson Correlation Coefficient(PCC) and Spearman Correlation Coefficient(SCC) was   
