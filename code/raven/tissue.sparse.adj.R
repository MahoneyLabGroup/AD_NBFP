#This function turns the tissue edge list into a sparse
#adjacency matrix for faster subsetting.

tissue.sparse.adj <- function(tissue.net){

	#sort the indices numerically so that we only fill in 
	#the top triangle of the adjacency matrix
	
	switch.idx = which(tissue.net[ ,1] < tissue.net[ ,2])
	tissue.net[switch.idx , 1:2] = tissue.net[switch.idx , c(2, 1)]
	
	#Sort the labels for the adjacency matrix numerically
	adj.lab = sort(union(unique(tissue.net[ ,1]), unique(tissue.net[ ,2])))
	
	n.row = length(adj.lab)
	adj <- sparseMatrix(i = match(tissue.net[ ,1], adj.lab), j = match(tissue.net[ ,2], adj.lab), 
	x = tissue.net[ ,3], dims = c(n.row, n.row), symmetric=TRUE, dimnames = list(NULL, adj.lab))
	return(adj)	
	
}