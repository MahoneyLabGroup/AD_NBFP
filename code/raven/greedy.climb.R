#This function takes a final data object from an SA run 
#and does a greedy hill climb to find a yet better network

greedy.climb <- function(data.obj, Gij = NULL, cur.gene = NULL){
	
	Lij <- data.obj$Lij
	locus.genes <- data.obj$locus.genes
	adj.mat <- data.obj$adj.mat
	fit.fun <- data.obj$fit.fun
	
	#initialize the gene network randomly
	if(is.null(Gij)){
		Gij.orig <- initialize.Gij(Lij, locus.genes, adj.mat)
		Gij <- Gij.orig$Gij
		cur.gene <- Gij.orig$current.gene.selection
		# imageWithText(Gij, show.text = FALSE)
		}
	
	#set up object for use in simulated annealing
	data.obj <- list(Gij, Lij, locus.genes, cur.gene, full.adj.mat, fit.fun)
	names(data.obj) <- c("Gij", "Lij", "locus.genes", "cur.gene", "adj.mat", "fit.fun")

	orig.fitness <- fitness(data.obj)
	last.fitness <- orig.fitness
	new.data.obj <- data.obj
	final.fitness <- orig.fitness - 1
	locus.genes <- data.obj$locus.genes

	iteration <- 1
	while(final.fitness < orig.fitness){
		cat("round", iteration, "\t")
		orig.fitness <- fitness(new.data.obj)
		for(i in 1:length(locus.genes)){
			report.progress(i, length(locus.genes))
			#update the network to each new gene in the locus
			#and pick the one that givest the maximum fitness
			if(nrow(locus.genes[[i]]) > 0){
				orig.g.fit <- fitness(new.data.obj)
				all.locus.fit <- rep(NA, nrow(locus.genes[[i]]))

				test.nets <- lapply(1:nrow(locus.genes[[i]]), function(x) modify.net.specified(new.data.obj, i, x))
				all.locus.fit <- unlist(lapply(test.nets, fitness))
				gene.weights <- as.numeric(locus.genes[[i]][,2])
				all.locus.fit <- all.locus.fit*gene.weights
				
				min.g <- which.min(all.locus.fit) #keep the new data.obj at the minimum energy level, if it is less than before we came into the loop
				if(all.locus.fit[min.g] < orig.g.fit){
					new.data.obj <- modify.net.specified(new.data.obj, i, min.g)
					}
				}
			} #end looping through loci
		iteration = iteration + 1
		final.fitness <- fitness(new.data.obj)
		cat("\t", final.fitness, "\n")
		}

	return(new.data.obj)
	
}