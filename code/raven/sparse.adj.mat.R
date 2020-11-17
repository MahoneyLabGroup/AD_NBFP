#This function builds a sparse adacency matrix from the
#FNTM edge list for use in TRiAGE. The result of this
#function is used as an argument in triage.one.interaction()
#and rank.pairs.one.interaction()


sparse.adj.mat <- function(tissue.net){
	
	test <- sparseMatrix()
	
	
	# Reorder entrez ids so that they are in numerical order (triangular form)
	switch.idx = which(tissue.net[ , 1] < tissue.net[ , 2])
	tissue.net[switch.idx , 1:2] = tissue.net[switch.idx , c(2, 1)]

	# Adjacency matrix for network
	adj.lab = sort(union(unique(tissue.net[ , 1]), unique(tissue.net[ , 2])))

	n.row = length(adj.lab)
	adj <- sparseMatrix(i = match(tissue.net[ , 1], adj.lab), j = match(tissue.net[ , 2], adj.lab), x = tissue.net[ , 3], dims = c(nrow, nrow), symmetric=TRUE, dimnames = list(NULL, adj.lab))

	rown
	return(adj)
	
}