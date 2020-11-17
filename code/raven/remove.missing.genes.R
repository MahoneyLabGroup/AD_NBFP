#This function removes genes from the locus list that
#are not in the adjacency matrix

remove.missing.genes <- function(locus.genes, adj.mat){
	
	all.locus.genes <- as.vector(unlist(lapply(locus.genes, function(x) x[,1])))
	missing.genes <- setdiff(all.locus.genes, colnames(adj.mat))

	while(length(missing.genes) > 0){
		missing.locale <- lapply(locus.genes, function(x) which(x[,1] == missing.genes[1]))
		missing.idx <- which(unlist(lapply(missing.locale, length)) > 0)
		if(length(missing.idx) > 0){
			orig.genes <- locus.genes[[names(missing.idx)]]
			missing.pos <- which(orig.genes[,1] == missing.genes[1])
			new.genes <- orig.genes[-missing.pos,,drop=FALSE]
			locus.genes[[names(missing.idx)]] <- new.genes
			}
		all.locus.genes <- as.vector(unlist(lapply(locus.genes, function(x) x[,1])))
		missing.genes <- setdiff(all.locus.genes, colnames(adj.mat))
		}

	return(locus.genes)
	}