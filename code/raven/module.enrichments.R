#This function gets module functional enrichments for a fntm network
#and assigns names to modules by selecting the top enrichment term
#for each. The names can also be set manually.
#adj.factor and num.modules only work for hierarchical clustering
#the adj.factor is the power you want to raise the adjacency matrix
#to, and num.modules is the number of modules you want to get out of
#the clustering

module.enrichments <- function(fntm.net, num.terms = 10, set.names.manually = FALSE, module.fun = c("fast_greedy", "hierarchical"), order.by = c("gprofiler", "p.value"), clust.method = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski") ,organism = "mmusculus", adj.factor = 1, num.modules = 2){
	
	clust.method = clust.method[1]
	order.by = order.by[1]
	
	require(igraph)
	if(module.fun == "fast_greedy"){
		comm <- cluster_fast_greedy(as.undirected(fntm.net))$membership
		comm.order <- order(comm)
		}
		
	if(module.fun == "hierarchical"){
		adj.mat <- as.matrix(as_adjacency_matrix(fntm.net, attr = "weight"))
		dist.mat <- dist(adj.mat^adj.factor, method = clust.method)
		clust <- hclust(dist.mat)
		comm <- cutree(clust, k = num.modules)
		comm.order <- clust$order
		}

	u_comm <- sort(unique(comm))
	
	comm.genes <- lapply(u_comm, function(x) V(fntm.net)$name[which(comm == x)])

	cat("looking up network enrichments...\n")
	mod.enrich <- lapply(comm.genes, function(x) gprofiler(x, organism, max_p_value = 0.05))
	
	for(i in 1:length(mod.enrich)){
		write.table(mod.enrich[[i]], paste0("Enrichment.Module.", i, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
		}
	
	mod.names <- unlist(lapply(mod.enrich, function(x) x[1,"term.name"]))
	
	if(set.names.manually){
		for(i in 1:length(u_comm)){
			
			if(order.by == "p.value"){
				ordered.enrich <- mod.enrich[[i]][order(mod.enrich[[i]][,"p.value"]),]
				}else{
				ordered.enrich <- mod.enrich[[i]]
				}
			print(ordered.enrich[1:num.terms,c("term.name", "term.size", "query.size", "overlap.size", "p.value")])
			choice <- readline(prompt = "Type the number of the term desired, or the name of the module.\n")
			is.num <- try(!is.na(as.numeric(choice)), silent = TRUE)
			if(is.num){
				mod.names[i] <- ordered.enrich[as.numeric(choice),"term.name"]
				}else{
				mod.names[i] <- choice
				}
			}
		}
	
	names(comm.genes) <- mod.names
	final.result <- list(comm.genes, comm.order)
	return(final.result)
}