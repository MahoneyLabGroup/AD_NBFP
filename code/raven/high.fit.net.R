#This function returns the highest fitness network
#which may be different than the consensus network

high.fit.net <- function(ntrials, SA = TRUE){
	
	if(SA){
		result <- readRDS("Final.Network.Sim.Ann.1.RData")
		n.loci <- length(result[[2]]$locus.genes)
		}else{
		result <- readRDS("Final.Network.Greedy1.RData")
		n.loci <- length(result$locus.genes)
		}
	
	all.fit <- rep(NA, ntrials)
	all.nets <- matrix(NA, nrow = ntrials, ncol = n.loci)

	for(tr in 1:ntrials){
		if(SA){
			result <- readRDS(paste("Final.Network.Sim.Ann.", tr, ".RData", sep = ""))
			final.data.obj <- result[[2]]
			}else{
			final.data.obj <- readRDS(paste("Final.Network.Greedy", tr, ".RData", sep = ""))
			}
		all.fit[tr] <- fitness(final.data.obj)
		all.nets[tr,] <- get.net(final.data.obj)	
		}
	
	max.fit <- min(all.fit)
	fit.locale <- which(all.fit == max.fit)

	max.fit.net <- all.nets[fit.locale,,drop=FALSE]
	u_nets <- unique(max.fit.net)
	colnames(u_nets) <- names(final.data.obj$locus.genes)
	
	final.result <- list("fitness" = max.fit, "networks" = u_nets)
	return(final.result)
	
	
}