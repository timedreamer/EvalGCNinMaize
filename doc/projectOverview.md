
# Maize Gene Co-expression Network Project Overview

The aim of this project is to build maize GCN using several different methods and test which methods is the best from RNA-Seq data.

1.  Choose publicly available RNA-Seq libraries in SRA.

2.  Download SRA files from NCBI using wget command.

3.  Convert .sra files into .fastq files using fastq-dump from SRA toolkit

4.  Cut adapters using Cutadapt. Since all libraries are from standard illumina sequencing, so the same adapter sequence were applied here (AGATCGGAAGAGC). Taken from Trim Galore manual.
5. The processed fastq files were quality checked with FASTQC.

6. HISAT2 was used in the allignment step to get bam files. Then gene counts were calculated from bam files using Feature Count.

7. Also, allignment-free kallisto was used to get transcript level expression. But I don't think I will use that data.
> Above steps were done by  pauper serve or KIN2077 XUbuntu machine. In order to streamline all above steps, *_snakemake_* workflow management was used here.

8. After FeactureCount step, all gene count files were combined togther. Libraries with less than 5M reads were discarded.

9. The genecount file was then imported to R environment for the following steps.

10. We used four different normalization method for the genecount data. 1) log2 normalized raw count. 2) VST normalization using DESEQ2 pacakge. 3) CPM using EDGER. 4) RPKM using EDGER.

11. In order to build GCN, five correlaiton methods and four mutual information methods were applied to four normalization results, including PCC, SCC, GCC, KCC, Bicor, AA, MA, MRNET and CLR.

12. The created GCNs were evaluated by MaizeCyc, maize protein-protein interaction database, Gene Ontology as well as ChIP-Seq data(HDA101). Results were calculated by Area Under the Curve.

13. Then the best method from correlation and mutual information method was chose to test the power of aggreated GCN.

14. Report the result.

15. From Ballouz et al(2015) paper, they reported better performance when aggregating networks from different experiments together. Also they found a positive correlatin between the sample size and the performance of a network when using PCC/SCC. So I also want to test this in our maize samples.
> In order to save time, six experiments in maize varying in sample sizes(12,36,65,,108,270, 404) was chosen.

16. PCC, MRNET and CLR were used to calculate their co-expression matrix. Then each network performance was evaluated by PP_PTY as well as Gene Ontology, reported as Area Under the Curve. The aggregating network which add the absolute value of edge's weight together was also evaluated.

17. The highly connected genes in network are considered as "hub genes". I also extracted hub genes from cpm_total and mrnet_total networks and did Gene Ontology enrichment analysis.

18. How about the performance of the newest proteomics network from Walley et al(2016) paper.
