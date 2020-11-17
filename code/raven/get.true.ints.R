get.true.ints <- function(Lij, net.genes, ko.ints){
	
	if(!is.null(dim(net.genes))){
		net.names <- colnames(net.genes)
		net.genes <- as.vector(net.genes)
		names(net.genes) <- net.names
		}

	
	#find all gene pairs in the network
	edge.list <- which(Lij == 1, arr.ind = TRUE)
	edge.list[,1] <- rownames(Lij)[as.numeric(edge.list[,1])]
	edge.list[,2] <- rownames(Lij)[as.numeric(edge.list[,2])]	
	gene.pairs <- t(apply(edge.list, 1, function(x) c(net.genes[which(names(net.genes) == x[1])], net.genes[which(names(net.genes) == x[2])])))
	
	gene.pair.data <- vector(mode = "list", length = nrow(gene.pairs))
	for(i in 1:nrow(gene.pairs)){
		gene1.locus <- which(ko.ints[,1:2] == gene.pairs[i,1], arr.ind = TRUE)
		gene2.locus <- which(ko.ints[,1:2] == gene.pairs[i,2], arr.ind = TRUE)
		pair.locus <- intersect(gene1.locus[,1], gene2.locus[,1])
		if(length(pair.locus) > 0){
			gene.pair.data[[i]] <- ko.ints[pair.locus,,drop=FALSE]
			}
		}
	
	int.table <- Reduce("rbind", lapply(gene.pair.data, function(x) x[which.min(as.numeric(x[,3])),]))
	
	return(int.table)


	
}