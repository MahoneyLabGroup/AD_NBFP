#This function gets the context of a gene subnetwork
#within a tissue network from IMP, or FNTM, or GIANT
#based on a provided edge list

get.network.context3 <- function(gene.edges, tissue.net, mart, start.thresh = 0.1, step.size = 0.01){
	
	require(igraph)
	all.nodes <- unique(c(gene.edges[,1], gene.edges[,2]))

	#find the subnetwork that includes all nodes		
	all.edges <- unique(unlist(lapply(all.nodes, function(x) c(which(tissue.net[,1] == x), which(tissue.net[,2] == x)))))
	sub.net <- tissue.net[all.edges,]
	# nrow(sub.net)
	#reset: start.thresh = 0.1; step.size = 0.01; chunk.size = 100; nodes.present = TRUE; edges.present = TRUE
		
		
	edges <- sub.net[,1:2]
	edges <- apply(edges, 2, as.character)
	
	net <- graph_from_edgelist(edges, directed = FALSE)
	E(net)$weight <- as.numeric(sub.net[,3])
	
	missing.genes <- setdiff(as.vector(gene.edges), V(net)$name)
	
	if(length(missing.genes) > 0){
		cat("Removing", length(missing.genes), "that are not in the tissue network.\n")
		cat(missing.genes, sep = "\n")
		missing.gene.idx <- Reduce("rbind", lapply(missing.genes, function(x) which(gene.edges == x, arr.ind = TRUE)))
		#take out the rows for these genes
		gene.edges <- gene.edges[-unique(missing.gene.idx[,1]),,drop=FALSE]
		all.nodes <- unique(as.vector(gene.edges))
		}
	
	
	#==================================================================
	#internal functions
	#==================================================================
	check.net <- function(net.now){
		all.paths <- apply(gene.edges, 1, function(x) get.shortest.paths(net.now, from = x[1], to = x[2], mode = "all", weights = 1-E(net.now)$weight))
		path.dist <- unlist(lapply(all.paths, function(x) length(x$vpath)))
		all.edges.present <- all(path.dist > 0)
		all.nodes.present <- all(unlist(lapply(lapply(all.nodes, function(x) grep(x, net.now)), length)) > 0)
		return(c(all.edges.present, all.nodes.present))
		}
	#==================================================================
		
	net.check <- check.net(net)
	edges.present <- net.check[1]
	nodes.present <- net.check[2]		
	new.thresh <- start.thresh
	
	while(edges.present && nodes.present && step.size >= 0.0001){
		cat(new.thresh, "\n")
		low.edges <- which(E(net)$weight < new.thresh)
		if(length(low.edges) > 0){
			cat("Deleting", length(low.edges), "edges\n")	
			test.net <- delete.edges(net, low.edges)
			no.edges <- which(degree(test.net) == 0)
			
			if(length(no.edges) > 0){
				cat("Deleting", length(no.edges), "unconnected vertices\n")
				test.net <- delete.vertices(test.net, no.edges)
				}

			net.check <- check.net(test.net)
			edges.present <- net.check[1]
			nodes.present <- net.check[2]
			if(edges.present && nodes.present){
				net <- test.net
				}else{
				new.thresh <- new.thresh - step.size
				step.size = step.size/10
				edges.present <- TRUE	
				nodes.present <- TRUE
				}
			}
		new.thresh <- new.thresh + step.size
		}
		

	#go through edge by edge trimming unnecessary connections 
	#starting with the lowest weights. Take big chunks and 
	#reduce as necessary
	u_weights <- sort(unique(E(net)$weight))
	chunk.size = 100	
	while(edges.present && nodes.present && chunk.size >= 1){
		u_weights <- sort(unique(E(net)$weight))
		weight.set <- which(E(net)$weight <= u_weights[chunk.size])

		test.net <- delete.edges(net, weight.set)
		no.edges <- which(degree(test.net) == 0)
		
		if(length(no.edges) > 0){
			cat("\tDeleting", length(no.edges), "unconnected vertex\n")
			test.net <- delete.vertices(test.net, no.edges)
			}

		net.check <- check.net(test.net)
		edges.present <- net.check[1]
		nodes.present <- net.check[2]
	
		if(edges.present && nodes.present){
			net <- test.net
			}else{
			chunk.size <- round(chunk.size/2)
			edges.present <- TRUE
			nodes.present <- TRUE
			cat("New chunk size:", chunk.size, "\n")
			}
		}
	
	# vcount(net)
	# ecount(net)
	
	v.names <- V(net)$name

	gene.names <- getBM(c("external_gene_name", "entrezgene"), "entrezgene", v.names, mart)

	name.locale <- match(V(net)$name, as.character(gene.names[,2]))
	node.names <- gene.names[name.locale,1]
	node.names[which(is.na(node.names))] <- v.names[which(is.na(node.names))]
	V(net)$name <- node.names
	V(net)$entrez.name <- v.names
		
	return(net)

}