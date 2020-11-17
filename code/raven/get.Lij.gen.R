#this function returns the locus-locus
#adjacency matrix

get.Lij.gen <- function(qtl.edgelist){
	
	require(igraph)
	
	net <- graph_from_edgelist(as.matrix(qtl.edgelist[,1:2]))	
	adj.mat <- as_adjacency_matrix(net)
	return(adj.mat)
	
	}