# MaizeCyc process method

In order to use pathway(**MaizeCyc**) to evaluate different networks, I need to parse the orignal "pathway.col" file.

File was download from **Ensemble**(http://goo.gl/vMTV6J). Then Excel was used to parse the file.

1. transpose the file.
2. For those elements with gene names instead of GRMZM numbers, find the associated GRMZM numbers in MaizeGDB.
3. If associated GRMZM number existed, then change gene name into GRMZM number. If not, leave it blank.
4. Only leave pathways with more than two genes(that's how co-expression can be evaluated).
5. In excel, added a extra row showing how many genes in each pathway. Order columns by the number of genes in the pathway, from small to big.
6. The result file was saved as "pathways.col_final_orginzed.txt".
7. Use R to convert transcript level into gene level. "pathway_process.R".
8. Final readable file is "pathway_namedChanged.txt".
