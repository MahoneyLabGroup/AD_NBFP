#This function finds summary weights between modules 
#in a network. It returns a network in which each node
#is a module, and the edges between them are the summary
#weight between them. The summary function can be any 
#function name that can be applied to a vector of numbers,
#like mean, median, max, min, etc.

module.weights <- function(net, modules, summ.fun = "mean"){
	
	u_modules <- sort(unique(modules))
	module.locale <- lapply(u_modules, function(x) which(modules == x))
	mod.pairs <- pair.matrix(u_modules)

	net.adj <- as_adjacency_matrix(net, type = "both", attr = "weight")

	edge.weights <- rep(NA, nrow(mod.pairs))
	for(i in 1:nrow(mod.pairs)){
		mod1.locale <- module.locale[[mod.pairs[i,1]]]
		mod2.locale <- module.locale[[mod.pairs[i,2]]]
		between.edges <- net.adj[mod1.locale, mod2.locale]
		fun <- call(summ.fun, between.edges)
		edge.weights[i] <- eval(fun)
		}	

	module.net <- graph_from_edgelist(mod.pairs, directed = FALSE)
	E(module.net)$weight <- edge.weights
	return(module.net)

	
	
	
	
}