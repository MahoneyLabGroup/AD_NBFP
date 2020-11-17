#max.temp.with.max.reject is the number of temperatures overall at which we have
#hit the maximum number of rejections
#once we have hit max rejections in this many temperatures, we stop

sim.ann <- function(initial.input, fn.opt, fn.step, start.temp = 10, max.iter = 10000, max.trials.at.temp = 100, max.rejects = 10, max.temp.with.max.reject = 3, verbose = TRUE){

	if(verbose){	
		plot.new()
		plot.window(xlim = c(0, max.trials.at.temp), ylim = c(sum(initial.input$Lij)*-1, 0))
		axis(1);axis(2)
		}
	
	opt.fun <- match.fun(fn.opt) 
	step.fun <- match.fun(fn.step)
	
	E <- fn.opt(initial.input) #find the starting point for the function we are optimizing
	
	iter <- 1	#initialize the number of iterations
	temp <- start.temp	#initialize the temperature at the start temp
	num.temp.with.max.rejects <- 0 #initialize the number of temps with the maximum number of rejects

	E.tracker <- E #start off our E tracker with the initial E
	temp.tracker <- temp #start off the temp tracker with the initial temp
	
	E.opt <- E #our starting E is our optimal E so far
	obj.opt <- initial.input #our initial object is our optimal object so far
	
	while(iter < max.iter && num.temp.with.max.rejects < max.temp.with.max.reject){

		if(verbose){
			cat("\n\nIteration: ", iter, "\n")
			}
			
			
		trials.at.temp <- 1 #initialize the number of tials we have had at this temperature
		rejects <- 0 #initialize the number of rejections at this temperature
	
		if(verbose){
			cat("\tTemp: ", temp, "\n") #report the current temperature
			 }
		
		#while we haven't hit any of our stopping criteria for the current temperature
		while(trials.at.temp < max.trials.at.temp & rejects < max.rejects){ 
			# if(verbose){		
				# report.progress(trials.at.temp, max.trials.at.temp) #print out our progress
				# }

			stepped.obj <- step.fun(data.obj = initial.input, verbose = TRUE)	#step the object according to the given function
			# stepped.obj <- modify.net.weighted(data.obj = initial.input, verbose = TRUE) #for debugging
			
			E.prime <- opt.fun(data.obj = stepped.obj)	#find the new value of the function you are trying to optimize
			
			if(E.prime < E.opt){		#if E.prime is lower than the current optimum,
				E.opt <- E.prime		#keep the E 
				obj.opt <- stepped.obj	#and the current object state
				if(verbose){points(trials.at.temp, E.opt, col = iter)}
				
				}else{ #otherwise, calculate a probability for accepting the step based on the temperature

				trans.prob <- exp((E-E.prime)/temp)	#find the transition probability
				rnd <- runif(1)	#flip a coin
			
				if(trans.prob >= rnd){	#if the transition probability is greater than your coin flip,

					initial.input <- stepped.obj	#accept the transition, by setting our initial input to the new value
					E <- E.prime					#and E to the new value
					if(verbose){points(trials.at.temp, E)}
					
					}else{ #otherwise reject the transition, and do not record E or the stepped object

					rejects <- rejects + 1
					# if(verbose){
						# cat(".", rejects, ".") #write out the number of rejections we have made
						# }
										
					} #end case for reject
				
		
				trials.at.temp <- trials.at.temp + 1 #increment the number of trials we have done at this temperature
				}#end case for if E.prime is worse (greater) E 
				
			} #end loop for before we have hit any of our stopping criteria


		initial.input <- obj.opt #set the input for the next temperature set to be the optimum object we found in this set
		E <- E.opt #set E to be the optimum E we found this round

		E.tracker <- c(E.tracker, E) #keep track of E at the end of each temperature set
		temp.tracker <- c(temp.tracker, temp) #keep track of the temp

		temp <- temp*0.8 	#adjust the temperature
		iter <- iter + 1	#increment the iterations
		
		cat("Current Fitness:", E, "\n")
		if(verbose){points(trials.at.temp, E, col = "red", pch = 16)}
		
		#if we have hit our maximum for rejections at a temperature
		#increment the number of temperatures at which we have hit
		#that maximum.
		if(rejects == max.rejects){
			num.temp.with.max.rejects <- num.temp.with.max.rejects + 1
			# if(verbose){
				# cat("\n", num.temp.with.max.rejects, "steps with max rejects")
				# }
			}
					
		
		
		} #end loop for before we've hit our max iterations
		
	# cat("iter: ", iter, "\n")
	# cat("trials.at.temp: ", trials.at.temp, "\n")
	# cat("rejects: ", rejects, "\n")
	# cat("temps.with.max.rejects: ", num.temp.with.max.rejects, "\n")
	cat("Final fitness:", fitness(initial.input), "\n")
	results <- list(E, initial.input, temp.tracker, E.tracker)
	names(results) <- c("E", "optimized.output", "Temp_over_time", "E_over_time")
	
	return(results)
	
	
	
} #end sim.ann function