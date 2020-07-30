# pairwise_genes_comparison
Compare individual CDS from two genomes. 


### Provide path to the filename

filename1 <- paste("data/SARS-Cov-2-Italy.txt")

filename2 <- paste(""data/SARS-CoV-2-Wuhan.txt"")


### Run function 
pairwise_genes_comparison(filename1,filename2)



###Notes
There are two important things that should be ketp in mind.

1) Identifier : This function uses "Gene" as an identifier to exract gene name from the fasta names. But, this is not 
consistent across different files. Other variations like locus tag, gene_id, protein_id are used to identify loci. So, 
there are two ways to get around this. Either you can replace locus_tag or gene_id with gene in your file. Or, you can change those keys within the package. Download the codes and change gene identifier in the line 154 and line 164. I have included a comment within the code, where the change needs to be made.


2)Identity: Sometimes two strains does not have similar number of genes. Then the package will throw an error or it might not even run if the number of genes are not same. You can remove the genes that does not have pair. So, the number of genes are the same and they are paired properly. You can look at the final table to confirm if the pairs was not formed or not by looking the name of the gene. Gene_Name1 and Gene_Name2 should be identical or should refer to the same gene.


