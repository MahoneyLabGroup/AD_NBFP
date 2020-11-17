plot.fitness.climb <- function(ntrials, ymax = NULL){
	
	all.fit <- vector(mode = "list", length = ntrials)
	all.temp <- vector(mode = "list", length = ntrials)
	for(tr in 1:ntrials){
		result <- readRDS(paste("Final.Network.Sim.Ann.", tr, ".RData", sep = ""))
		all.temp[[tr]] <- result$Temp_over_time
		all.fit[[tr]] <- result$E_over_time
		}

	max.fit <- max(unlist(lapply(all.fit, function(x) max(abs(x)))))
	max.runs <- max(unlist(lapply(all.fit, function(x) length(x)))) 
	
	if(is.null(ymax)){
		# ymax <- max.fit + 10 - (max.fit%%10) #round up to the nearest 10
		ymax <- max.fit
		}
	
	temp.mat <- list2Matrix(all.temp)
	mean.temps <- apply(temp.mat, 2, function(x) mean(x, na.rm = TRUE))

	plot.new()
	plot.window(xlim = c(0, max.runs), ylim = c(0, ymax))
	for(i in 1:length(all.fit)){
		points(abs(all.fit[[i]]), type = "l")
		}
	axis(1, at = 1:length(mean.temps), labels = signif(mean.temps, 2))
	axis(2)
	mtext("Temperature", side = 1, line = 2.5)
	mtext("Fitness", side = 2, line = 2.5)		
	return(all.fit)	
}