#get edge weights from each node in the phenotype ontology (from Mouse Mine)
#to each other gene in the genome (all genes in the downloaded tissue.net)
#the gene list should be in entrez IDs

tissue.adj.mat <- function(tissue.net, gene.list, inc.all.genes = TRUE, 
remove.unconnected = TRUE, verbose = TRUE){

	if(class(gene.list) == "list"){
		gene.list <- as.vector(unlist(lapply(gene.list, function(x) x[,1])))
		}

	test.num <- suppressWarnings(try(as.numeric(gene.list[1]), silent = TRUE))

	require(igraph)
	
	if(!is.na(test.num)){
		gene.list <- as.numeric(gene.list)
		}
	
	if(verbose){cat("Filtering tissue network...\n")}
	gene.list <- unique(gene.list[which(!is.na(gene.list))])
	in.gene.list1 <- which(tissue.net[,1] %in% gene.list)
	in.gene.list2 <- which(tissue.net[,2] %in% gene.list)

	if(inc.all.genes){
		#filter edge list so all edges including
		#our genes are incorporated
		all.edges <- unique(c(in.gene.list1, in.gene.list2))
		}else{
		#filter edge list so all edges are between
		#genes in the gene list
		all.edges <- intersect(in.gene.list1, in.gene.list2)
		}	
	
	net.edges <- tissue.net[all.edges,]

	# test.ind <- 1:10
	test.ind <- 1:nrow(net.edges)

	if(verbose){cat("Converting node names to characters...\n")}
	edges <- apply(as.matrix(net.edges[test.ind,1:2]), 2, as.character)
	if(verbose){cat("Converting edge list to network...\n")}
	net <- graph.edgelist(edges[,1:2], directed = FALSE)
	E(net)$weight <- as.numeric(net.edges[test.ind,3])
	if(verbose){cat("Converting network to adjacency matrix...\n")}
	net.adj <- as_adj(net, attr = "weight")

	#we now have a square adjacency matrix that includes
	#all genes. Filter so that only the gene.list genes
	#are in the rows of the matrix.
	if(inc.all.genes){
		gene.list.locale <- which(row.names(net.adj) %in% gene.list)
		net.adj <- net.adj[gene.list.locale,]
		}
		
	#igraph will return a sparse matrix, which is great, but
	#we need regular matrix behaviors, so convert it back to
	#a regular matrix
	net.adj <- as.matrix(net.adj)

	#now check the matrix to make sure all genes from our original list
	#are in the matrix. If not, add them with all 0 edges
	#these are genes that aren't in the top edges of the tissue network
	#we downloaded.
	missing.genes <- setdiff(gene.list, rownames(net.adj))
	if(length(missing.genes) > 0){
		missing.rows <- matrix(0, nrow = length(missing.genes), ncol = ncol(net.adj))
		rownames(missing.rows) <- missing.genes
		net.adj <- rbind(net.adj, missing.rows)
		missing.cols <- matrix(0, ncol = length(missing.genes), nrow = nrow(net.adj))
		colnames(missing.cols) <- missing.genes
		net.adj <- cbind(net.adj, missing.cols)
		}

	if(remove.unconnected){
		#remove unconnected genes. This function also orders
		#the matrix such that TP genes are in the first block
		net.adj <- remove.unconnected.tp(tissue.mat = net.adj)
		}
	
	return(net.adj)
	
}