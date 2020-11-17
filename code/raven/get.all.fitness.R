get.all.fitness <- function(ntrials, SA = TRUE){
	
	all.fit <- rep(NA, ntrials)
	for(tr in 1:ntrials){
		if(SA){
			result <- readRDS(paste("Final.Network.Sim.Ann.", tr, ".RData", sep = ""))
			final.data.obj <- result[[2]]
			}else{
			final.data.obj <- readRDS(paste("Final.Network.Greedy", tr, ".RData", sep = ""))
			}
		all.fit[tr] <- fitness(final.data.obj)
		}
	
	return(all.fit)
	
	
}