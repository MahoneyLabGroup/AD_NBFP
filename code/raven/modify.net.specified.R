#this function adjusts the network
#by picking a vertex and selecting
#a new gene in the locus
#Lij is the locus network from the cape
#results.
#Gij is the adjacency matrix holding
#the functional connections between 
#randomly selected genes in the locus
#network.
#locus genes is the list holding the possible
#genes in each locus from the cape network.
#adj.mat is the full adjacency matrix from the
#functional network from which new values are 
#drawn to fill Gij.
#cur.gene is the vector holing which genes are
#currently selected for each block.

modify.net.specified <- function(data.obj, locus.idx, gene.idx, verbose = FALSE){
	

	for(i in 1:length(data.obj)){
		assign(names(data.obj)[i], data.obj[[i]])
		}	
		
		
	just.int <- as.matrix(Lij[,1:nrow(Lij)])

	#figure out which edges we can modify based on the cape networks
	edges <- which(just.int != 0, arr.ind = TRUE)

	
	block <- names(locus.genes)[locus.idx]

	the.gene <- locus.genes[[block]][gene.idx,1]
	#change the current gene to reflect the update
	cur.gene[block] <- the.gene
	# print(the.gene)
	#now update any edges that the new gene participates in
	gene.targets <- names(which(Lij[block,] != 0))
	if(length(gene.targets) > 0){
		for(i in 1:length(gene.targets)){
			gene.selection <- cur.gene[gene.targets[i]]
			if(is.na(gene.selection)){
				new.edge <- 0
				}else{
				new.edge <- adj.mat[as.character(the.gene), as.character(gene.selection)]
				}
			if(verbose){
				# cat("\t", Gij[block, gene.targets[i]], "->", new.edge, "\n")	
				}
			Gij[block, gene.targets[i]] <- new.edge
			}
		}
	
	gene.sources <- names(which(Lij[,block] != 0))
	if(length(gene.sources) > 0){
		for(i in 1:length(gene.sources)){
			gene.selection <- cur.gene[gene.sources[i]]
			if(is.na(gene.selection)){
				new.edge <- 0
				}else{
				new.edge <- adj.mat[as.character(the.gene), as.character(gene.selection)]
				}
			if(verbose){
				# cat("\t", Gij[gene.sources[i], block], "->", new.edge, "\n")
				}
			Gij[gene.sources[i], block] <- new.edge
			}
		}
	
	results <- list(Gij, Lij, locus.genes, cur.gene, adj.mat, fit.fun)
	names(results) <- c("Gij", "Lij", "locus.genes", "cur.gene", "adj.mat", "fit.fun")
	return(results)
	
	
}