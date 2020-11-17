#random network
plot.random.net <- function(n, p, vertex.color = "lightblue", vertex.size = 20, edge.width = 4){
	library(igraph)
	net <- erdos.renyi.game(n, p)
	plot(net, vertex.label = NA, vertex.color = vertex.color, vertex.size = vertex.size, edge.width = edge.width)
	}