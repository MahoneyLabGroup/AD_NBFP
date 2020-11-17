#This function fits a random gene-gene network to the locus-locus network
#it returns the fitness of each network

fit.rnd.net <- function(Lij, locus.genes, adj.mat, n.trials = 100){
	
	rnd.fit <- rep(NA, trials)
for(tr in 1:trials){
	report.progress(tr, trials)
	Gij.orig <- initialize.Gij(Lij, locus.genes, adj.mat = full.adj.mat)
	Gij <- Gij.orig$Gij
	cur.gene <- Gij.orig$current.gene.selection
	# imageWithText(Gij, show.text = FALSE)
	
	#set up object for use in simulated annealing
	data.obj <- list(Gij, Lij, locus.genes, cur.gene, full.adj.mat, "sum")
	names(data.obj) <- c("Gij", "Lij", "locus.genes", "cur.gene", "adj.mat", "fit.fun")
	
	rnd.fit[tr] <- fitness(data.obj)
	}

	return(rnd.fit)
	
	
	
}