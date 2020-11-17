#This function plots a random modular graph

plot.modular.graph <- function(n.vert = 81, p = 0.5, n.clust = 9, cols = brewer.pal(8, "Pastel1")){
	
	require(RColorBrewer)
	
	nodes.per.clust <- round(n.vert/n.clust)
	adj <- matrix(0, nrow = n.vert, ncol = n.vert)
	
	clust.idx <- colV <- vector(mode = "list", length = n.clust)
	#form clusters
	start.idx <- 1
	for(i in 1:n.clust){
		max.node <- min(c(start.idx + nodes.per.clust - 1, n.vert))
		nodes <- start.idx:max.node
		clust.idx[[i]] <- nodes
		colV[[i]] <- rep(cols[i], length(nodes))
		rnd.net <- erdos.renyi.game(length(nodes), p)
		adj[nodes, nodes] <- as.matrix(as_adjacency_matrix(rnd.net))
		start.idx <- start.idx + nodes.per.clust
	}
	
	#set connections between clusters
	top.net <- as.matrix(as_adjacency_matrix(erdos.renyi.game(n.clust, p*1.5)))
	clust.connect <- which(top.net == 1, arr.ind = TRUE)
	for(i in 1:nrow(clust.connect)){
		clust1.id <- clust.connect[i,1]
		clust2.id <- clust.connect[i,2]
		#add some connections between clusters
		rnd.connect <- as.matrix(as_adjacency_matrix(erdos.renyi.game(length(clust.idx[[clust1.id]]), p/20)))
		adj[clust.idx[[clust1.id]], clust.idx[[clust2.id]]] <- rnd.connect
	}
	
	
	final.net <- graph_from_adjacency_matrix(adj)
	final.net <- as.undirected(final.net)
	plot(final.net, vertex.color = unlist(colV), vertex.label = NA, vertex.size = 8)
}
