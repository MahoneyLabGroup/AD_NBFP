#This function gets the context of a gene subnetwork
#within a tissue network from IMP, or FNTM, or GIANT
#It gets the subnetwork that includes all genes one
#step away from the query genes at the given edge 
#threshold


get.network.context2 <- function(net.genes, tissue.net, mart, min.edge.weight = 0.15){
	
	require(igraph)
	all.nodes <- net.genes
		
	all.edges <- unique(unlist(lapply(all.nodes, function(x) c(which(tissue.net[,1] == x), which(tissue.net[,2] == x)))))
	sub.net <- tissue.net[all.edges,]
	nrow(sub.net)
		
	#get the subgraph 	
	#filter the network for only the edges that meet the minimum edge weight standard
	thresh.net <- sub.net[which(as.numeric(sub.net[,3]) >= min.edge.weight),]
	# nrow(thresh.net)
	
	edges <- thresh.net[,1:2]
	edges <- apply(edges, 2, as.character)
	
	net <- graph_from_edgelist(edges, directed = FALSE)
	E(net)$weight <- as.numeric(thresh.net[,3])
	
	adj.mat <- as.matrix(as_adjacency_matrix(net, attr = "weight"))
		
	
	#=========================================================	
	# internal functions
    #=========================================================
	get.direct.connections <- function(gene.name){
		gene.locale <- which(colnames(adj.mat) == gene.name)
		if(length(gene.locale) > 0){
			gene.edges <- adj.mat[gene.locale,]
			non.zero <- which(gene.edges != 0)
			gene.connectors <- colnames(adj.mat)[non.zero]	
			return(gene.connectors)
			}else{
			return(NULL)
			}
		}
    #=========================================================

	all.connectors <- lapply(net.genes, get.direct.connections)
	num.neighbors <- unlist(lapply(all.connectors, length))
		
	all.connectors <- unique(c(unlist(all.connectors), net.genes))
	
	missing.genes <- setdiff(all.connectors, colnames(adj.mat))
	
	if(length(missing.genes) > 0){
		cat("Removing", length(missing.genes), "that are not in the tissue network.\n")
		cat(missing.genes, sep = "\n")
		missing.locale <- match(missing.genes, all.connectors)
		all.connectors <- all.connectors[-missing.locale]
		}
	
	sub.mat <- adj.mat[all.connectors, all.connectors]

	# pheatmap(sub.mat)

	trimmed.net <- graph_from_adjacency_matrix(sub.mat, mode = "undirected", weighted = TRUE)
	# clust <- cluster_fast_greedy(trimmed.net)$membership
	# u_clust <- sort(unique(clust))

	
	v.names <- V(trimmed.net)$name

	
	gene.names <- getBM(c("external_gene_name", "entrezgene"), "entrezgene", v.names, mart)

	name.locale <- match(V(trimmed.net)$name, as.character(gene.names[,2]))
	node.names <- gene.names[name.locale,1]
	node.names[which(is.na(node.names))] <- v.names[which(is.na(node.names))]
	V(trimmed.net)$name <- node.names
	V(trimmed.net)$entrez.name <- v.names
		
	return(trimmed.net)

}