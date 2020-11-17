#This function gets the maximum fitness
#network out of a quadprog result object

qp.max.net <- function(qp.result, locus.genes, lib){
	require(biomaRt)

	max.net <- matrix(NA, ncol = 1, nrow = length(locus.genes))
	rownames(max.net) <- names(locus.genes)
	
	all.genes <- unlist(locus.genes)
	
	for(i in 1:length(locus.genes)){
		gene.locale <- match(locus.genes[[i]], all.genes)
		max.gene <- which.max(qp.result$solution[gene.locale])
		max.net[i,1] <- all.genes[gene.locale[max.gene]]
		}

	cat("Looking up gene names...\n")
	gene.names <- getBM(c("external_gene_name", "entrezgene"), "entrezgene", max.net[,1], lib)
	
	gene.name.locale <- match(max.net[,1], gene.names[,2])
	# cbind(max.net[,1], gene.names[gene.name.locale,2])
	max.net <- cbind(max.net, gene.names[gene.name.locale,1])
	return(max.net)
	
}