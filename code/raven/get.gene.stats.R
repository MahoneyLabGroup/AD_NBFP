#For both SA and greedy trials, this function 
#returns the number of genes per locus, a table
#identifying the number of times each of the top
#gene was selected for the locus, ranked tables 
#of the unique genes selected for each locus, the
#number of unique genes in the whole network, and 
#the consensus genes. It also plots
#the comparison of nodes per locus for the different 
#methods.
#mart is only required if entrez.2.name = TRUE

get.gene.stats <- function(Lij, ntrials, entrez.2.name = TRUE, mart = NULL){

	require(vioplot)
	genes.trial.SA <- get.genes.by.trial(Lij, ntrials, SA = TRUE)
	genes.trial.Gr <- get.genes.by.trial(Lij, ntrials, SA = FALSE)

	#the automatic simplify in apply can be an issue here when there are the same number of 
	#genes found in all loci, so I'm am doing something a bit more complicated to avoid that
	#issue
	gene.counts.SA <- lapply(1:nrow(genes.trial.SA), function(x) sort(table(genes.trial.SA[x,]), decreasing =TRUE))
	names(gene.counts.SA) <- rownames(genes.trial.SA)
	
	gene.counts.Gr <- lapply(1:nrow(genes.trial.Gr), function(x) sort(table(genes.trial.SA[x,]), decreasing =TRUE))
	names(gene.counts.Gr) <- rownames(genes.trial.Gr)


	num.genes.per.locus <- unlist(lapply(locus.genes, length))
	num.genes.by.SA <- unlist(lapply(names(locus.genes), function(x) length(gene.counts.SA[[which(names(gene.counts.SA) == x)]])))
	num.genes.by.Gr <- unlist(lapply(names(locus.genes), function(x) length(gene.counts.Gr[[which(names(gene.counts.Gr) == x)]])))
	gene.call.table <- cbind(num.genes.per.locus, num.genes.by.SA, num.genes.by.Gr)

	gene.per.locus.table <- rbind(num.genes.per.locus, num.genes.by.SA, num.genes.by.Gr)

	final.results <- list("gene.per.locus.table" = gene.per.locus.table)

	null.SA <- which(unlist(lapply(gene.counts.SA, length) == 0))
	if(length(null.SA) > 0){
		gene.counts.SA <- gene.counts.SA[-null.SA]
		}

	null.Gr <- which(unlist(lapply(gene.counts.Gr, length) == 0))
	if(length(null.Gr) > 0){
		gene.counts.Gr <- gene.counts.Gr[-null.Gr]
		}

	rank.table.SA <- list2Matrix(gene.counts.SA)
	rank.table.Gr <- list2Matrix(gene.counts.Gr)

	final.results$rank.table.SA <- rank.table.SA
	final.results$rank.table.Gr <- rank.table.Gr	
	
	
	if(entrez.2.name){
		id.table.SA <- entrez2name(list2Matrix(lapply(gene.counts.SA, names)), mart)
		id.table.Gr <- entrez2name(list2Matrix(lapply(gene.counts.Gr, names)), mart)
		}else{
		id.table.SA <- list2Matrix(lapply(gene.counts.SA, names))
		id.table.Gr <- list2Matrix(lapply(gene.counts.Gr, names))
		}

	final.results$id.table.SA <- id.table.SA
	final.results$id.table.Gr <- id.table.Gr
	

	#the number of genes selected by each regime
	#the -1 is for NA
	final.results$unique.genes.SA <- length(unique(as.vector(id.table.SA))) - 1
	final.results$unique.genes.Gr <- length(unique(as.vector(id.table.Gr))) - 1

	if(entrez.2.name){
		final.results$consensus.genes.SA <- entrez2name(unlist(lapply(gene.counts.SA, function(x) names(x)[1])), mart)
		final.results$consensus.genes.Gr <- entrez2name(unlist(lapply(gene.counts.Gr, function(x) names(x)[1])), mart)
		}else{
		final.results$consensus.genes.SA <- unlist(lapply(gene.counts.SA, function(x) names(x)[1]))
		final.results$consensus.genes.Gr <- unlist(lapply(gene.counts.Gr, function(x) names(x)[1]))
		}
	
	return(final.results)
	
	}