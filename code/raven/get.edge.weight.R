#This function spot-checks interaction weights
#for checking adjacency matrix correctness.

get.edge.weight <- function(tissue.net, node1, node2){
	
	node1.locale <- which(tissue.net == node1, arr.ind = TRUE)
	node2.locale <- which(tissue.net == node2, arr.ind = TRUE)
	
	edge.locale <- intersect(node1.locale[,1], node2.locale[,1])
	# print(tissue.net[edge.locale,])

	return(tissue.net[edge.locale,])	
}