#This is a weighted version of modify.net
#It uses the temperature from simulated
#annealing to place a probability distribution
#over the choices of gene in the randomly
#selected locus.
#This function requires that the current temp
#be appended to the data.obj

modify.net.weighted <- function(data.obj, verbose = FALSE){
	
	for(i in 1:length(data.obj)){
		assign(names(data.obj)[i], data.obj[[i]])
		}
		
	#pick a locus at random, but make sure the
	#locus has more than one gene available to select
	num.genes <- unlist(lapply(locus.genes, nrow))
	multiple.genes <- which(num.genes > 1)	

	rnd.sample <- as.vector(sample(multiple.genes, 1))
	rnd.block <- names(locus.genes)[rnd.sample]

	orig.fit <- fitness(data.obj)
	
	#select each gene in the locus and calculate the fitness
	#change associated with the modification
	all.genes <- locus.genes[[rnd.block]][,1]
	gene.weights <- as.numeric(locus.genes[[rnd.block]][,2])

	test.nets <- lapply(1:length(all.genes), function(x) modify.net.specified(data.obj, locus.idx = rnd.sample, gene.idx = x, verbose = FALSE))
	all.fit <- unlist(lapply(test.nets, fitness))
	
	if(all(all.fit == 0)){
		all.fit <- rep(1, length(test.nets))
		}	
	
	#now incorporate the weights from the individual genes in the locus
	#trait-related genes should be bumped up in their fitness over
	#non-trait-related genes
	weighted.fit <- all.fit*gene.weights
	
	# plot(weighted.fit, all.fit, xlim = c(-1,0), ylim = c(-1,0))
	# abline(0,1)
		
	#normalize all.fit so we aren't raising e to enormous powers
	max.abs <- max(abs(weighted.fit))
	all.fit <- weighted.fit/max.abs
	orig.fit <- orig.fit/max.abs
	
	#cbind(all.genes, all.fit)
	#place a probability distribution over this fitness 
	#using the current temperature from simulated annealing

	fit.weights <- exp((orig.fit - weighted.fit)/(temp/weight.power)) #the SA version
	# plot(abs(all.fit), fit.weights) # a positive relationship between fitness change and weights
			
			
	# if(any(!is.finite(fit.weights))){
		# cat("returning test nets\n")
		# return(test.nets)
		# }
			
	#in trying to debug I will make a near greedy set of weights
	#here is a near greedy set of weights
	#some weighs are 0, which we don't want in the final version
	# fit.weights.test <- orig.fit - weighted.fit
	# fit.weights.test[which(fit.weights.test < 0)] <- 0

	# plot(fit.weights, fit.weights.test)
	
	#This is the proper normalization
	if(sum(fit.weights) > 0){
		fit.prob <- fit.weights/sum(fit.weights) #normalize the weights 
		}else{
		fit.prob <- rep(1/length(fit.weights), length(fit.weights))	
		}

	# plot(fit.prob, type = "h")
	
	# if(sum(fit.weights.test) > 0){
		# fit.prob.test <- fit.weights.test/sum(fit.weights.test) #normalize the weights 
		# }else{
		# fit.prob.test <- rep(1/length(fit.weights.test), length(fit.weights.test))	
		# }
	
	
	# plot(fit.prob, fit.prob.test)
	
	#This is a greedy normalization
	# fit.prob <- rep(0, length(fit.weights))
	# fit.prob[which.max(fit.weights)] <- 1
	
	# plot(fit.prob, type = "h")
	# plot(fit.weights, fit.prob) # a positive relationship between probability and weight

	#pick a new gene for this vertex based
	#on the probability distribution
		
	#if there is only one gene in the block, pick it
	if(nrow(locus.genes[[rnd.block]]) == 1){
		rnd.gene = locus.genes[[rnd.block]][1,1]
		}
		
	#otherwise, select a gene using the probability 
	#distribution calculated for the 
	if(nrow(locus.genes[[rnd.block]]) > 1){
		rnd.gene <- sample(locus.genes[[rnd.block]][,1], 1, prob = fit.prob)
		#change the current gene to reflect the update
		
		# rnd.gene <- sample(locus.genes[[rnd.block]], 10000, replace = TRUE, prob = fit.prob)
		# times.chosen <- unlist(lapply(all.genes, function(x) length(which(rnd.gene == x))))/length(all.genes)
		# plot(times.chosen, fit.prob) #positive relationship between probability and times chosen
		
		}

	if(length(locus.genes[[rnd.block]]) == 0){
		if(verbose){cat("Region selected does not contain genes. Returning unchanged network.\n")}
		results <- list(Gij, Lij, locus.genes, cur.gene, adj.mat, fit.fun, temp)
		names(results) <- c("Gij", "Lij", "locus.genes", "cur.gene", "adj.mat", "fit.fun", "temp")
		return(results)
		}	
		
	cur.gene[rnd.block] <- rnd.gene
		
	#now update any edges that the new gene participates in
	if(verbose){	cat("Updating edges...\n")}
	other.gene.targets <- names(which(Lij[rnd.block,] != 0))
	if(length(other.gene.targets) > 0){
		for(i in 1:length(other.gene.targets)){
			gene.selection <- cur.gene[other.gene.targets[i]]
			if(is.na(gene.selection)){
				new.edge <- 0
				}else{
				new.edge <- adj.mat[as.character(rnd.gene), as.character(gene.selection)]
				}
			if(verbose){
				cat("\t", Gij[rnd.block, other.gene.targets[i]], "->", new.edge, "\n")	
				}
			Gij[rnd.block, other.gene.targets[i]] <- new.edge
			}
		}
	
	other.gene.sources <- names(which(Lij[,rnd.block] != 0))
	if(length(other.gene.sources) > 0){
		for(i in 1:length(other.gene.sources)){
			gene.selection <- cur.gene[other.gene.sources[i]]
			if(is.na(gene.selection)){
				new.edge <- 0
				}else{
				new.edge <- adj.mat[as.character(rnd.gene), as.character(gene.selection)]
				}
			if(verbose){
				cat("\t", Gij[other.gene.sources[i], rnd.block], "->", new.edge, "\n")
				}
			Gij[other.gene.sources[i], rnd.block] <- new.edge
			}
		}
	
	results <- list(Gij, Lij, locus.genes, cur.gene, adj.mat, fit.fun, temp, weight.power)
	names(results) <- c("Gij", "Lij", "locus.genes", "cur.gene", "adj.mat", "fit.fun", "temp", "weight.power")
	# }
	return(results)
	
	
}