#this function finds the true positive interactions
#between pairs of loci. Data should be entered as a 
#list of loci with genes of interest as a vector in
#each entry of the list.
#if locus.gene.list has more than 2 entries, Lij should
#be added as well to define the interactions. Otherwise
#the script will look at all possible pairs
#if ph.genes (a vector of phenotype-related genes) is 
#included, this function also checks to see how many of
#the true positive interactions are between phenotype-
#related genes

get.tp.pairs <- function(locus.gene.list, tp.edges, Lij = NULL, ph.genes = NULL){
	
	require(igraph)
	net <- graph_from_edgelist((as.matrix(tp.edges[,1:2])), directed = FALSE)
	adj.mat <- as.matrix(as_adjacency_matrix(net))
	
	
	is.pair.tp <- function(gene.pair){
		g1 <- which(rownames(adj.mat) == gene.pair[1])
		g2 <- which(colnames(adj.mat) == gene.pair[2])
		is.tp <- adj.mat[g1,g2] > 0
		return(is.tp)
		}
	
	is.pair.ph <- function(gene.pair){
		g1 <- length(which(ph.genes == gene.pair[1]))
		g2 <- length(which(ph.genes == gene.pair[2]))
		is.ph <- g1+g2 > 0
		return(is.ph)
		}
	
	if(!is.null(Lij)){
		int.locale <- which(Lij > 0, arr.ind = TRUE)
		}else{
		int.locale <- pair.matrix(1:length(locus.gene.list))
		}
		
	tp.mat <- matrix(NA, nrow = nrow(int.locale), ncol = 2)
	colnames(tp.mat) <- c("num.TP", "total")
	
	for(i in 1:nrow(int.locale)){
		report.progress(i, nrow(int.locale))
		l1 <- rownames(Lij)[int.locale[i,1]]
		l2 <- rownames(Lij)[int.locale[i,2]]
		l1.locale <- which(names(locus.gene.list) == l1)
		l2.locale <- which(names(locus.gene.list) == l2)
		genes.l1 <- locus.gene.list[[l1.locale]]
		genes.l2 <- locus.gene.list[[l2.locale]]
		all.pairs <- unique(cbind(rep(genes.l1, length(genes.l2)), rep(genes.l2, each = length(genes.l1))))
		tp.pairs <- apply(all.pairs, 1, is.pair.tp)
		ph.pairs <- apply(all.pairs, 1, is.pair.ph)
		tp.mat[i,1] <- length(which(tp.pairs))
		tp.mat[i,2] <- length(tp.pairs)
		}
	
	frac.tp <- tp.mat[,1]/tp.mat[,2]
	# hist(frac.tp)
	final.mat <- cbind(tp.mat, frac.tp)
	return(final.mat)
	
	}