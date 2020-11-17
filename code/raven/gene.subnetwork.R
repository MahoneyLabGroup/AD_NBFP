#This function finds a tissue subnetwork associated
#with a given entrez id


gene.subnetwork <- function(entrezid, tissue.adj.mat, num.steps = 1, return.gene.names = TRUE, organism = c("human", "mouse")){
	
	library(biomaRt)
	library(igraph)
	
	gene.locale <- which(colnames(tissue.adj.mat) == entrezid)
	
	if(length(gene.locale) == 0){
		return("That gene is not in this network.")
		}
		
	connectors <- c(entrezid, colnames(tissue.adj.mat)[which(tissue.adj.mat[,gene.locale] > 0)])
	
	if(num.steps > 1){
		for(i in 2:num.steps){
			connector.locale <- match(connectors, colnames(tissue.adj.mat))
			new.connectors <- lapply(connector.locale, function(x) which(tissue.adj.mat[,x] > 0))
			new.connector.names <- unique(unlist(lapply(new.connectors, function(x) names(x))))
			connectors <- c(connectors, new.connector.names)
			}
		}
	
	connector.locale <- match(connectors, colnames(tissue.adj.mat))
	sub.adj <- as.matrix(tissue.adj.mat[,connector.locale])
	row.connections <- rowSums(sub.adj)
	has.connections <- which(row.connections > 0)
	sub.adj <- sub.adj[has.connections,]
	
	if(return.gene.names){
		if(organism[1] == "mouse"){		
		mart = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
		}else{
		mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
		}
	
		row.genes <- getBM(c("entrezgene", "external_gene_name"), "entrezgene", rownames(sub.adj), mart = mart)
		col.genes <- getBM(c("entrezgene", "external_gene_name"), "entrezgene", colnames(sub.adj), mart = mart)
		
		gene.order <- match(rownames(sub.adj), row.genes[,1])
		rownames(sub.adj) <- row.genes[gene.order,2]
		
		gene.order <- match(colnames(sub.adj), col.genes[,1])
		colnames(sub.adj) <- col.genes[gene.order,2]	
	
		#take out any rows/columns that do not have gene names associated with them
		#these have been discontinued
		
		not.na.col <- which(!is.na(colnames(sub.adj)))
		sub.adj <- sub.adj[,not.na.col]
		not.na.row <- which(!is.na(rownames(sub.adj)))
		sub.adj <- sub.adj[not.na.row,]
		
		}
	

	#pare down to a square matrix
	common.names <- intersect(rownames(sub.adj), colnames(sub.adj))
	common.row.locale <- match(common.names, rownames(sub.adj))

	sub.adj <- sub.adj[common.row.locale,,drop=FALSE]
	common.col.locale <- match(common.names, colnames(sub.adj))
	sub.adj <- sub.adj[,common.col.locale,drop=FALSE]

	sub.graph <- graph_from_adjacency_matrix(sub.adj, "undirected", weighted = TRUE)	


	return(sub.graph)
	
	
	
}