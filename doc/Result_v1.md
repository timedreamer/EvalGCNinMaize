# Result



## Point one:

### The best method for building Gene Co-Expression Network(GCN) in maize

Within recent years, the use of next-generation sequencing is growing rapidly especially for plant biology. In maize studies, the application of illumina sequencing increased from 6 samples in 2008 to 2389 samples in 2015. On the contrary, microarry,the previous favorable genome-wide gene expression measurement technology, dropped to only14 samples in maize when considering two most used platform (GPL4032 and GPL12620). With the decreasing cost of next-generation sequencing, the improvement of sequencing techniques and the amount of data gained, it is highly likely that next-generation sequencing will continue to be the mainstream genome-wide analysis of gene expression method.  

Except working for the researchers own purposes, with large mount of expression data available, we are also able to build Gene Co-expression Network (GCN).  . One of usage of large amount of gene-expression data is to construct Gene Co-expression Network. Several maize co-expression studies have been reported. But most of them were based on microarray data which is limited in the genes that can be detected as well as probe bias. The large amount of genome-wide RNA-seq data makes it possible to build maize GCN and it will very helpful for the maize community.

In this project we used 1266 high quality RNA-seq samples to build maize GCN. These samples are publicly available from NCBI SRA (Supplemental Table1). We manually annotated their original tissue and genotype from either SRA annotations or their published paper. The downloaded raw files were adapter removed, quality checked and aligned to maize reference genome(V3) by HISAT2. The gene-level count data for all samples were calculated by FeatureCount. The gene-level count matrix was then used in R environment for further steps.

Due to the characteristic of illumina sequencing, several normalization methods have been reported including VST, CPM, RPKM rawcount







Point two: The aggregation network performs better than single network.



Points three: The usage of GCN maize
