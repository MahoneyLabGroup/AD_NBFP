#This function determines whether each 
#identified gene interaction exists
#in the true network
# true.p = readRDS('~/Documents/Data/yeast/Costanzo/Data File S1. Raw genetic interaction datasets- Pair-wise interaction format/AllIntP.RData')
# true.e = readRDS('~/Documents/Data/yeast/Costanzo/Data File S1. Raw genetic interaction datasets- Pair-wise interaction format/AllIntE.RData')

true.gene.interactions <- function(Lij, all.locus.genes, true.p, true.e){
	require(igraph)
	
	edge.locale <- which(Lij != 0, arr.ind = TRUE)

	gene.int.val <- function(gene1, gene2, p.or.e){
		all.p <- matrix(NA, nrow = 1, ncol = length(true.p))
		colnames(all.p) <- names(true.p)
		for(i in 1:length(true.p)){
			gene1.locale <- which(rownames(true.p[[i]]) == gene1)
			gene2.locale <- which(rownames(true.p[[i]]) == gene2)
			if(length(gene1.locale) > 0 && length(gene2.locale) > 0){
				if(p.or.e == "p"){
					all.p[1,i] <- true.p[[i]][gene1.locale, gene2.locale]
					}else{
					all.p[1,i] <- true.e[[i]][gene1.locale, gene2.locale]	
					}
				}
			}
		return(all.p)
		}

	all.gene.pair.p <- all.gene.pair.e <- vector(mode = "list", length = nrow(edge.locale))	
	int.names <- apply(cbind(rownames(Lij)[edge.locale[,1]],rownames(Lij)[edge.locale[,2]]), 1, function(x) paste(x[1], x[2], sep = "-"))
	names(all.gene.pair.p) <- names(all.gene.pair.e) <- int.names
	
	for(i in 1:nrow(edge.locale)){
		report.progress(i, nrow(edge.locale))
		source.name <- rownames(Lij)[edge.locale[i,1]]
		target.name <- rownames(Lij)[edge.locale[i,2]]
		source.genes <- all.locus.genes[[source.name]]
		target.gene <- all.locus.genes[[target.name]]
		all.gene.pairs <- cbind(rep(source.genes, length(target.genes)), rep(target.genes, each = length(source.genes)))
		
		#get the edge weight for each pair from each experiment set
		all.gene.pair.p[[i]] <- t(apply(all.gene.pairs, 1, function(x) gene.int.val(x[1], x[2], "p")))
		colnames(all.gene.pair.p) <- names(all.true.nets)
		
		all.gene.pair.e[[i]] <- t(apply(all.gene.pairs, 1, function(x) gene.int.val(x[1], x[2], "e")))
		colnames(all.gene.pair.e) <- names(all.true.nets)
		}

	return(list("gene.pair.e" = all.gene.pair.e, "gene.pair.p" = all.gene.pair.p))
	
	}