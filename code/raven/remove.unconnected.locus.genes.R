#This function removes genes from locus genes
#that are not represented in the adjacency matrix
#this happens when locus genes are unconnected
#and therefore removed from the adjacency matrix

remove.unconnected.locus.genes <- function(locus.genes, tissue.mat){
	
	unconnected.genes <- setdiff(unlist(locus.genes), rownames(tissue.mat))
	if(length(unconnected.genes) == 0){
		return(locus.genes)
		}
	
	cat("Removing", length(unconnected.genes), "genes...\n")	
	for(i in 1:length(unconnected.genes)){
		gene.locus <- grep(unconnected.genes[i], locus.genes)
		gene.idx <- which(locus.genes[[gene.locus]] == unconnected.genes[i])
		locus.genes[[gene.locus]] <- locus.genes[[gene.locus]][-gene.idx]
		}
		
	return(locus.genes)
	
	
}