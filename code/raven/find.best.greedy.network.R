#Gij and cur.gene can be specified. If not,
#they will be initialized randomly

find.best.greedy.network <- function(Lij, locus.genes, adj.mat, Gij = NULL, cur.gene = NULL, trials, fit.fun = "sum"){

	if(is.null(Gij) && !is.null(cur.gene)){
		stop("cur.gene is specified but Gij is not.")
		}

	if(!is.null(Gij) && is.null(cur.gene)){
		stop("Gij is specified but cur.gene is not.")
		}

	init.fit <- rep(NA, trials)
	greedy.fit <- rep(NA, trials)
	for(tr in 1:trials){
		print(tr)
		
		#set up object for use in simulated annealing
		data.obj <- list(Lij, locus.genes, full.adj.mat, fit.fun)
		names(data.obj) <- c("Lij", "locus.genes", "adj.mat", "fit.fun")
		init.fit[tr] <- fitness(data.obj)
		
		final.data.obj <- greedy.climb(data.obj, Gij, cur.gene)
		saveRDS(final.data.obj, paste("Final.Network.Greedy", tr, ".RData", sep = ""))
		
		greedy.fit[tr] <- fitness(final.data.obj)
		
		}
	# boxplot(list(abs(init.fit), abs(greedy.fit)), names = c("Random", "Greedy"))
	return(list("init.fit" = init.fit, "greedy.fit" = greedy.fit))
	
}