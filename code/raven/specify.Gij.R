#This function builds a Gij network from specific gene inputs


specify.Gij <- function(Lij, cur.gene, adj.mat, verbose = FALSE){
	
	Gij <- matrix(0, nrow = nrow(Lij), ncol = ncol(Lij))
	rownames(Gij) <- rownames(Lij)
	colnames(Gij) <- colnames(Lij)

	
	for(i in 1:length(cur.gene)){
		if(!is.na(cur.gene[i])){
			block <- names(cur.gene)[i]
			the.gene <- cur.gene[i]
			
			#now update any edges that the new gene participates in
			gene.targets <- names(which(Lij[block,] != 0))
			gene.targets <- gene.targets[which(!is.na(cur.gene[gene.targets]))]
			if(length(gene.targets) > 0){
				for(i in 1:length(gene.targets)){
					gene.selection <- cur.gene[gene.targets[i]]
					new.edge <- adj.mat[as.character(the.gene), as.character(gene.selection)]
					if(verbose){
						cat("\t", Gij[block, gene.targets[i]], "->", new.edge, "\n")	
						}
					Gij[block, gene.targets[i]] <- new.edge
					}
				}
			
			gene.sources <- names(which(Lij[,block] != 0))
			gene.sources <- gene.sources[which(!is.na(cur.gene[gene.sources]))]
			if(length(gene.sources) > 0){
				for(i in 1:length(gene.sources)){
					gene.selection <- cur.gene[gene.sources[i]]
					new.edge <- adj.mat[as.character(the.gene), as.character(gene.selection)]
					if(verbose){
						cat("\t", Gij[gene.sources[i], block], "->", new.edge, "\n")
						}
					Gij[gene.sources[i], block] <- new.edge
					}
				}
			}	
		} #end looping through genes
	# image(Gij)
	return(Gij)
	
}
