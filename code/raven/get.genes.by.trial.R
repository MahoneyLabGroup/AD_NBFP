#This function builds a table of the genes 
#selected in each trial of an algorithm

get.genes.by.trial <- function(Lij, ntrials, SA = TRUE){
	

	block.genes.by.trial <- matrix(NA, nrow = nrow(Lij), ncol = ntrials)
	rownames(block.genes.by.trial) <- rownames(Lij)
	for(tr in 1:ntrials){
		if(SA){
			result <- readRDS(paste("Final.Network.Sim.Ann.", tr, ".RData", sep = ""))
			final.data.obj <- result[[2]]
			results.text <- "Sim.Ann"
			}else{
			final.data.obj <- readRDS(paste("Final.Network.Greedy", tr, ".RData", sep = ""))
			results.text <- "Greedy"
			}
		block.genes.by.trial[,tr] <- final.data.obj$cur.gene
		}
		return(block.genes.by.trial)

}