#This function gets the context of a gene subnetwork
#within a tissue network from IMP, or FNTM, or GIANT
#This will be the primary function, but I don't want 
#to lose the original functioning algorithm
#min.edge.weight = 0.15; num.additional.genes = 20; plot.net = FALSE; highlight.nodes = net.genes; layout.call = layout_nicely; layout.matrix = NULL; stretch.layout = 1.2; layout.x.adjust = 0.2; layout.y.adjust = 0.1; vertex.size = 10
get.network.context <- function(net.genes, tissue.net, mart, min.edge.weight = 0.15, add.genes.per.subnet = 10){
	
	require(igraph)
	

	#=========================================================	
	# internal functions
    #=========================================================
	get.subnetwork <- function(gene.names){
		subnet.gene.locale <- intersect(which(tissue.net[,1] %in% gene.names), which(tissue.net[,2] %in% gene.names))
		full.subnet <- tissue.net[subnet.gene.locale,]
		edge.list <- apply(as.matrix(full.subnet[,1:2]), 2, as.character)
		net <- graph_from_edgelist(edge.list, directed = FALSE)
		E(net)$weight <- full.subnet[,3]
		return(net)
		}
		
	gene.neighborhood <- function(gene.name, n.top.edges){
		locale1 <- which(tissue.net[,1] == gene.name)
		locale2 <- which(tissue.net[,2] == gene.name)
		locale <- c(locale1, locale2)
		if(length(locale) == 0){return(NULL)}
		subnet <- tissue.net[locale,]
		edge.min <- sort(subnet[,3], decreasing = TRUE)[n.top.edges]
		above.thresh <- which(subnet[,3] > edge.min)
		above.net <- subnet[above.thresh,]
		subnet.genes <- c(unique(c(above.net[,1], above.net[,2])), gene.name)
		net <- get.subnetwork(subnet.genes)
		return(net)
		}


	
	find.direct.edge <- function(gene1, gene2){
		locale1 <- intersect(which(tissue.net[,1] == gene1), which(tissue.net[,2] == gene2))
		locale2 <- intersect(which(tissue.net[,2] == gene1), which(tissue.net[,1] == gene2))
		locale <- unique(c(locale1, locale2))
		if(length(locale) == 0){
			locale <- "not found"
			}
		return(locale)
		}
    #=========================================================

	#get the neighborhoods around each query gene
	net.list <- lapply(net.genes, function(x) gene.neighborhood(x, add.genes.per.subnet))

	#get the subnetwork for all unique vertices in these networks
	u_verts <- unique(unlist(lapply(net.list, function(x) V(x)$name)))
	expanded.net <- get.subnetwork(u_verts)

	#trim the network to only edges above the given weight
	low.edges <- which(E(expanded.net)$weight < min.edge.weight)
	if(length(low.edges) > 0){
		trimmed.net <- delete.edges(expanded.net, low.edges)
		}else{
		trimmed.net <- expanded.net	
		}

	v.names <- V(trimmed.net)$name
	
	gene.names <- getBM(c("external_gene_name", "entrezgene"), "entrezgene", v.names, mart)

	name.locale <- match(V(trimmed.net)$name, as.character(gene.names[,2]))
	node.names <- gene.names[name.locale,1]
	node.names[which(is.na(node.names))] <- v.names[which(is.na(node.names))]
	V(trimmed.net)$name <- node.names
	V(trimmed.net)$entrez.name <- v.names
		
	return(trimmed.net)

}