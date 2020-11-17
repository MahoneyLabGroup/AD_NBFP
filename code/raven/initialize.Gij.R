
initialize.Gij <- function(Lij, locus.genes, adj.mat, verbose = FALSE){
	
	just.int <- as.matrix(Lij[,1:nrow(Lij)])
	
	Gij <- matrix(0, nrow = nrow(just.int), ncol = nrow(just.int))	
	rownames(Gij) <- rownames(just.int)
	colnames(Gij) <- rownames(just.int)
	
	#we also need to make an object that keeps track of the current
	#gene assignment for each block. Only one gene per block can be
	#selected at one time
	cur.gene <- rep(NA, nrow(just.int))
	names(cur.gene) <- rownames(just.int)
		
	for(b in 1:length(locus.genes)){
		block <- names(locus.genes)[b]
		
		if(length(locus.genes[[block]]) > 0){
			
			#pick a random gene for the block if it doesn't have a gene already
			if(is.na(cur.gene[block])){
				possible.genes <- locus.genes[[block]][,1]
				if(length(possible.genes) == 1){
					rnd.gene <- possible.genes
					}
				if(length(possible.genes) > 1){
					rnd.gene <- sample(possible.genes, 1)
					}
				#change the current gene to reflect the update
				cur.gene[block] <- rnd.gene
				}else{
				rnd.gene <- cur.gene[block]	
				}
		
			#now update any edges that the new gene participates in
			gene.targets <- names(which(Lij[block,] != 0))
			if(length(gene.targets) > 0){
				for(i in 1:length(gene.targets)){
					if(is.na(cur.gene[gene.targets[i]])){
						target.genes <- locus.genes[[gene.targets[i]]][,1]
						if(length(target.genes) == 1){
							gene.selection <- target.genes
							}
						if(length(target.genes) > 1){
							gene.selection <- sample(target.genes, 1)
							cur.gene[gene.targets[i]] <- gene.selection
							}
						}else{
						gene.selection <- cur.gene[gene.targets[i]]
						}

					new.edge <- try(adj.mat[as.character(rnd.gene), as.character(gene.selection)], silent = TRUE)
					if(class(new.edge) == "try-error"){new.edge = 0}

					if(verbose){
						cat("\t", Gij[block, gene.targets[i]], "->", new.edge, "\n")	
						}
					Gij[block, gene.targets[i]] <- new.edge
					}
				}
			
			gene.sources <- names(which(Lij[,block] != 0))
			if(length(gene.sources) > 0){
				for(i in 1:length(gene.sources)){
					if(is.na(cur.gene[gene.sources[i]])){
						source.genes <- locus.genes[[gene.sources[i]]][,1]
						if(length(source.genes) == 1){
							gene.selection <- source.genes
							}
						if(length(source.genes) > 1){
							gene.selection <- sample(source.genes, 1)
							cur.gene[gene.sources[i]] <- gene.selection
							}
						}else{
						gene.selection <- cur.gene[gene.sources[i]]
						}
					new.edge <- adj.mat[as.character(rnd.gene), as.character(gene.selection)]
					if(verbose){
						cat("\t", Gij[gene.sources[i], block], "->", new.edge, "\n")
						}
					Gij[gene.sources[i], block] <- new.edge
					}
				}
			}
		}
	
	results <- list(Gij, cur.gene)
	names(results) <- c("Gij", "current.gene.selection")
	return(results)
	
	
}