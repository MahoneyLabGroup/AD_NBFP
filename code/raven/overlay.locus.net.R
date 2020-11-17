#This function builds the adjusted gene-gene matrix
#for quadratic programing. All gene interaction values
#between non-interacting loci are zeroed out

overlay.locus.net <- function(full.adj.mat, Lij, locus.genes){

	require(igraph)	
	
	edge.list <- which(Lij != 0, arr.ind = TRUE)
	
	total.genes <- unlist(locus.genes)
	
	new.mat <- matrix(0, nrow = length(total.genes), ncol = length(total.genes))
	colnames(new.mat) <- rownames(new.mat) <- total.genes
	 # colnames(new.mat) <- rownames(new.mat) <- rownames(full.adj.mat)
	 	
	for(i in 1:nrow(edge.list)){
		locus.pair.genes <- c(locus.genes[[rownames(Lij)[edge.list[i,1]]]], locus.genes[[rownames(Lij)[edge.list[i,2]]]])
		sub.mat <- full.adj.mat[as.character(locus.pair.genes), as.character(locus.pair.genes)]
		new.mat[as.character(locus.pair.genes), as.character(locus.pair.genes)] <- sub.mat
		}
	
	# image(full.adj.mat, main = "Full Matrix")
	# quartz();image(new.mat, main = "Overlayed Matrix")
	# quartz();image(sub.mat, main = "Sub Matrix")	
	# heatmap(full.adj.mat)	
	# quartz();heatmap(new.mat)
	diag(new.mat) <- 0
	return(new.mat)

}