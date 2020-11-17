#This function uses simulated annealing to find the best gene network
#for the cape locus network based on edges in the FNTM/GIANT network
#trials = 1:100; start.temp = 8; max.iter = 100; max.trials.at.temp = 20000; max.rejects = 15000; max.temp.with.max.reject = 10; verbose = FALSE
#trials is actually the number of each trial you want to run. If you have already run
#50 trials out of 100 and do not want to re-run these, set trials to 51:100.
#you can specify a starting Gij and current gene selection, otherwise
#Gij and cur.gene will be initialized randomly
#weight.power is the power to raise the gene weights to when using weighted gene selection
#This makes the probablity distribution a little less uniform for faster increase in fitness
#the default is 1, which leaves the weights unchanged

find.best.gene.network <- function(Lij, locus.genes, adj.mat, Gij = NULL, cur.gene = NULL, mod.fun = modify.net.weighted, trials = 1:100, start.temp = 8, max.iter = 100, max.trials.at.temp = 20000, max.rejects = 15000, max.temp.with.max.reject = 10, weight.power = 1, verbose = FALSE, n.cores = 4){
	
	if(n.cores > 1){verbose = FALSE}
	
	if(is.null(Gij) && !is.null(cur.gene)){
		stop("cur.gene is specified but Gij is not.")
		}

	if(!is.null(Gij) && is.null(cur.gene)){
		stop("Gij is specified but cur.gene is not.")
		}
	
	
	run.trial.set <- function(trial.num){
		# cat("Trial", trial.num, "\n")
		for(tr in 1:length(trial.num)){
			#initialize the gene network randomly
			#if no Gij is provided
			if(is.null(Gij)){
				Gij.orig <- initialize.Gij(Lij, locus.genes, adj.mat = full.adj.mat)
				Gij <- Gij.orig$Gij
				cur.gene <- Gij.orig$current.gene.selection
				}
				
			# imageWithText(Gij, show.text = FALSE)
			# test.net <- graph_from_adjacency_matrix(abs(Gij), weighted = TRUE)
			# no.edge <- which(degree(test.net) == 0)
			# test.net <- delete.vertices(test.net, no.edge)
			# plot(test.net)
			
			#set up object for use in simulated annealing
			data.obj <- list(Gij, Lij, locus.genes, cur.gene, full.adj.mat, "sum", start.temp, weight.power)
			names(data.obj) <- c("Gij", "Lij", "locus.genes", "cur.gene", "adj.mat", "fit.fun", "temp", "weight.power")
			
			starting.fitness <- fitness(data.obj)
			
			result <- sim.ann(initial.input = data.obj, fn.opt = fitness, fn.step = mod.fun, start.temp = start.temp, max.iter = max.iter, max.trials.at.temp = max.trials.at.temp, max.rejects = max.rejects, max.temp.with.max.reject = max.temp.with.max.reject, verbose = verbose)
			# net <- graph_from_adjacency_matrix(abs(result[[2]]$Gij), weighted = TRUE, mode = "directed")
			# plot(net, vertex.size = 0.1, arrow.length = 0.1)
			
			saveRDS(result, paste("Final.Network.Sim.Ann.", trial.num[tr], ".RData", sep = ""))
			}
		}
		
		
	if(n.cores > 1){
		chunked.p <- chunkV(trials, n.cores)
		cl <- makeCluster(n.cores)
		registerDoParallel(cl)
		foreach(p = 1:length(chunked.p), .export = c("Gij")) %dopar% {
			run.trial.set(trial.num = chunked.p[[p]])
			}				
		stopCluster(cl)
		}else{
		for(i in trials){
			run.trial.set(i)
			}
		
		}


	
}