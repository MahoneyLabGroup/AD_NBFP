#This function finds all unique networks


all.nets <- function(ntrials, SA = TRUE){
	all.nets <- vector(mode = "list", length = ntrials)
	
	for(i in 1:ntrials){
		if(SA){
			data.obj <- readRDS(paste0("Final.Network.Sim.Ann.", i, ".RData"))
			all.nets[[i]] <- SA.final.network(data.obj)
			}	
		}
	net.mat <- Reduce("rbind", all.nets)
	
	return(net.mat)
}