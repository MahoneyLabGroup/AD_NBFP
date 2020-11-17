get.network.genes <- function(Gij, Lij, cur.gene, mart){
	
	
	edge.list <- which(Gij != 0, arr.ind = TRUE)	
	gene.mat <- matrix(NA, nrow = nrow(edge.list), ncol = 3)
	
	for(i in 1:nrow(edge.list)){
		block1 <- rownames(Gij)[edge.list[i,1]]
		gene1 <- cur.gene[block1]
		block2 <- rownames(Gij)[edge.list[i,2]]
		gene2 <- cur.gene[block2]

		cape.val <- Lij[block2, block1]
		if(cape.val == 0){
			cape.val <- Lij[block1, block2]
			}
		gene.mat[i,] <- c(gene1, gene2, cape.val)
		}
	
	u_genes <- unique(c(gene.mat[,1], gene.mat[,2]))
	gene.names <- getBM(c("entrezgene","external_gene_name"), "entrezgene", u_genes, mart)
	
	for(i in 1:nrow(gene.names)){
		gene.locale <- which(gene.mat == gene.names[i,1])
		gene.mat[gene.locale] <- gene.names[i,2]
		}
	
	return(gene.mat)
	
}