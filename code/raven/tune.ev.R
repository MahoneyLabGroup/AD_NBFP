#This function tunes the number of eigenvectors
#to take from a decomposed tissue matrix.


tune.ev <- function(full.mat, decomp.mat, ev.seq = c(seq(10, 40, 10), seq(50, 300, 50)), starting.cost.list = 10^seq(-5, 2, 1)){
	
	if(max(ev.seq) > nrow(decomp.mat)){stop("There are more eigenvectors specified than exist in the decomposed matrix. Pleast edit ev.seq.\n")}
	
	all.cost.list <- vector(mode = "list", length = length(ev.seq))
	fitnesses <- vector(mode = "list", length = length(ev.seq))
	names(fitnesses) <- names(all.cost.list) <- ev.seq

	get.fit.cost <- function(fitness.list){
		#first find the final trial
		trial.locale <- grep("Trial 1", fitness.list)
		final.list <- fitness.list[max(trial.locale):length(fitness.list)]
		split.fit <- strsplit(final.list,  ": ")
		
		acc <- unlist(lapply(split.fit, function(x) as.numeric(x[2])))
		max.acc <- which.max(acc)

		cost.part <- split.fit[[max.acc]]
		cost <- as.numeric(multi.strsplit(cost.part, c("C = ", " :"))[1])
		
		return(c("Cost" = cost, "Accuracy" = acc[max.acc]))
		}

	fit.cost.mat <- matrix(NA, ncol = 2, nrow = length(ev.seq))
	rownames(fit.cost.mat) <- ev.seq
	best.acc <- 0
	num.worse <- 0
	counter <- 1
	
	while(num.worse < 2 && ev.seq[counter] <= max(ev.seq)){
	
		cat("Number of EV:", ev.seq[counter], "\n")
		final.mat <- t(decomp.mat[,1:ev.seq[counter]]) #tp.mat has genes in columns and eigenvectors in rows
		rownames(final.mat) <- paste("EV_", 1:length(1:ev.seq[counter]))
		colnames(final.mat) <- colnames(full.mat)
	
		fitnesses[[counter]] <- capture.output(all.cost.list[[counter]] <- tune.cost.parameter(tissue.adj = final.mat, num.tp = nrow(full.mat), n.trials = 1, C.list = starting.cost.list, verbose = TRUE))
		
		fit.cost.mat[counter,] <- get.fit.cost(fitnesses[[counter]])
		cur.acc <- fit.cost.mat[counter,2]
		best.acc <- max(fit.cost.mat[,2], na.rm = TRUE)
		
		if(cur.acc < best.acc){
			num.worse <- num.worse + 1
			}
		if(cur.acc == best.acc){ #if we are at the maximum accuracy so far, reset the counter
			num.worse <- 0
			}

		# print(fit.cost.mat[1:counter,])		
		counter <- counter + 1
		}


		# plot(ev.seq, fit.cost.mat[,2], type = "l", lwd = 3, xlab = "Number of Eiegenvectors", ylab = "Maximum Accuracy")
		
		max.fit.idx <- which.max(fit.cost.mat[,2])
		best.ev <- ev.seq[max.fit.idx]
		tuned.cost <- all.cost.list[[max.fit.idx]]
		

		final.results <- list("best.ev.num" = best.ev, "tuned.cost.list" = tuned.cost, "cost.fit.table" = fit.cost.mat)
		return(final.results)
		}