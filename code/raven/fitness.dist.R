#This function initializes a network and modifies it many
#times to look at the distribution of fitness changes given
#a single gene modification


find.start.temp <- function(data.obj, start.temp = 8, cooling.factor = 0.8, slope.min = 0.1){
	
	
	for(i in 1:length(data.obj)){
		assign(names(data.obj)[i], data.obj[[i]])
		}
				
	#pick a locus at random, but make sure the
	#locus has more than one gene available to select
	num.genes <- unlist(lapply(locus.genes, length))
	multiple.genes <- which(num.genes > 1)	
	orig.fit <- fitness(data.obj)
	
	locus.temps <- list()
	locus.ent <- list()
	
	for(m in 1:length(multiple.genes)){
		cat(names(locus.genes)[m], "\n")
		
		block <- names(locus.genes)[m]
		all.genes <- locus.genes[[block]]
		
		test.nets <- lapply(1:length(all.genes), function(x) modify.net.specified(data.obj, locus.idx = m, gene.idx = x, verbose = FALSE))
	
		all.fit <- unlist(lapply(test.nets, fitness))	
	
		temp.ent <- 1
		all.ent <- NULL
		all.temp <- NULL
		temp <- start.temp
		while(is.finite(temp.ent) && round(temp.ent, 2) > 0){
			
			fit.weights <- exp((orig.fit - all.fit)/temp) #the SA version
			fit.prob <- fit.weights/sum(fit.weights) #normalize the weights 

			# plot(fit.prob, type = "h")
			temp.ent <- entropy(fit.prob)
			all.ent <- c(all.ent, temp.ent)
			all.temp <- c(all.temp, temp)
			temp <- temp*cooling.factor
			}
		
		locus.ent[[m]] <- all.ent
		locus.temps[[m]] <- all.temp
		}

	max.ent <- unlist(lapply(locus.ent, function(x) max(x, na.rm = TRUE)))
	num.genes <- lapply(locus.genes, length)
	
	
	par(mfrow = c(1,2))
	plot(num.genes, max.ent, xlab = "Number of Genes in Locus", ylab = "Maximum Entropy")
		
	all.slopes <- vector(mode = "list", length = length(locus.ent))
	for(i in 1:length(locus.ent)){
		all.slopes[[i]] <- instant.slope(locus.temps[[i]], locus.ent[[i]])
		}
	
	#where do the slopes flatten
	turning.point <- unlist(lapply(all.slopes, function(x) min(which(x > slope.min))))
	max.temp <- rep(NA, length(turning.point))
	for(i in 1:length(all.slopes)){
		consec.mat <- consec.pairs(1:length(locus.temps[[i]]))
		max.temp[i] <- locus.temps[[i]][consec.mat[turning.point[i],1]]
		}
	
	opt.start.temp <- max(max.temp)
	
	plot.new()
	plot.window(xlim = c(0, start.temp), ylim = c(0, max(max.ent)))
	for(i in 1:length(locus.ent)){
		points(locus.temps[[i]], locus.ent[[i]], type = "l")
		}	
	axis(1);axis(2)	
	mtext("Temperature", side = 1, line = 2)
	mtext("Entropy", side = 2, line = 2.5)
	mtext("Entropy changes with temperature", line = 2)

	
	abline(v = opt.start.temp, col = "red")
	text(x = opt.start.temp + (start.temp*0.01), y = max(max.ent)*0.1, labels = signif(opt.start.temp, 2), adj = 0, col = "red")

	return(opt.start.temp)
	
}